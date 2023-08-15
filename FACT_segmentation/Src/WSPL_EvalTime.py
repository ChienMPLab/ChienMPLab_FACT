import WSPL_func
from WSPL_ImageProcessing import *
import os
from os import path, sep, mkdir
import cv2
from cv2 import imread, imwrite
from glob import glob
import joblib
import json
import pandas as pd
import time


from evaluation import *
from skimage.measure import label, regionprops
from skimage.segmentation import relabel_sequential


def encodingMasks( labels):
    labels_flatten =labels.flatten()
    loc = np.flatnonzero( labels_flatten)
    value = labels_flatten[loc]

    result = np.zeros( (loc.size+1,2), loc.dtype)
    result[1:,0]= loc
    result[1:,1]= value
    result[0,:] = labels.shape
    return result

def decodinMasks( encodingLabels):
    imgshape = encodingLabels[0,:]
    imgtmp = np.zeros(imgshape, dtype=np.int64).flatten()
    imgtmp[encodingLabels[1:,0]]=  encodingLabels[1:,1]
    return imgtmp.reshape(imgshape)

def load_model( model_path):
    modelInfo_tmp_path = model_path.replace('.model','.json')
    cuml_model = joblib.load(model_path)

    with open(modelInfo_tmp_path, 'r') as openfile:
        modelInfo_tmp = json.load(openfile)

    if 'select' in modelInfo_tmp:
        select = np.array(modelInfo_tmp['select'])
    else:
        select = []
    sigmaarray = np.array(modelInfo_tmp['sigmaarray'])

    return cuml_model, sigmaarray, select

def segmentation_( cuml_model, sigmaarray, select, path_img):
    st = time.time()
    img = imread( path_img, -1).astype(np.float32)
    imgshape = img.shape[:2]
    if len(select)<1:
        # whole features
        imgarray = WSPL_func.ggenFS( path_img, sigmaarray)
    else:
        # subfeatures
        method_list, parameters_dic = subFeatures_temeplate( sigmaarray)
        imgarray = imgprocessing_subprocessingV2( img, 0, parameters_dic, select)
    imgarray = imgarray.reshape( (-1, imgarray.shape[-1]))
    y = WSPL_func.segmenation_GPU(imgarray,  cuml_model, 0)
    y = y.reshape(imgshape)
    label_img = WSPL_func.gpu_label( y, 0)
    props = WSPL_func.regionprops( label_img, img)
    # remove the small segmentation
    props_array = np.zeros( (len(props), 1), dtype=int) # area, x, y
    for i in range( len( props)):
        props_array[ i, 0] = props[i].area

    for i_label in np.flatnonzero( props_array[:,0] <= 10):
        x, y = props[i_label].coords[:,1], props[i_label].coords[:,0]
        label_img[ y,x ] = 0
    end = time.time()
    eval_time = end- st
    return label_img, eval_time

def eval_timming( path_model, folder_images, path_result):

    running = 1
    if not path.isfile( path_model):
        running =0
        print('Can not find the model!')
    else:
        cuml_model, sigmaarray, select = load_model( path_model)
    
    if not path.isdir( folder_images):
        running =0
        print('Can not find the folder of images!')
    else:
        imglist = sorted(
            glob(f'{folder_images}{sep}*.tif')
            )
        if len(imglist)< 1:
            running =0
            print('Can not find any images!')
    
    if running>0:
        result_segmentation = []
        # fileName, encodingMasks, eval_time
        for path_image in imglist:
            label_img, eval_time = segmentation_(
                cuml_model, sigmaarray, select, path_image  
                )
            img_masks_encoding = encodingMasks( label_img)
            result_segmentation.append(
                [
                    path.basename(path_image), img_masks_encoding, eval_time
                ]
            )
        df = pd.DataFrame(
            result_segmentation,
            columns=['filename', 'encodingMasks', 'eval_time']
            )
        df.to_pickle(path_result)# .pickle
    else:
        print('Did not eval the model!')


def eval_df_df( file_GT, file_seg):

    df_gt = pd.read_pickle(file_GT)
    df_seg = pd.read_pickle(file_seg)

    results = pd.DataFrame(
            columns=["Image", "Threshold", "F1", "Jaccard", "TP",\
                    "FP", "FN", "Official_Score", "Precision", "Recall"]
            )
    list_filename = df_seg['filename'].to_list()
    for i in range( len(df_gt)):
        filename = df_gt['filename'][i]+".tif"
        w,d,x,y = df_gt['wdxy'][i]
        img = decodinMasks( df_gt['encodingMasks'][i])

        index_seg = list_filename.index(filename)
        img_pred = decodinMasks( df_seg['encodingMasks'][index_seg])[y:y+d,x:x+w]
        ground_truth = label(img)
        label_img = label(img_pred)
        # Relabel objects
        ground_truth = relabel_sequential(ground_truth)[0] 
        prediction = relabel_sequential(label_img)[0]


        results = compute_af1_results(
                    ground_truth, 
                    prediction, 
                    results, 
                    filename
                )
    return results



# def eval_timming( folder_dataset, folder_GT, list_model):

#     folder_dataset = f'{folder_dataset}{sep}Train_data'
#     running = 1
#     if not path.exists(folder_dataset):
#         print('Can not find the folder!')
#         print(f'{folder_dataset}')
#         running = 0
    
#     imglist = sorted( glob(f'{folder_GT}{sep}GT{sep}raw_whole{sep}*.tif') )
#     if len(imglist)<1:
#         print('There is no image at the GT folder!')
#         print(f'{folder_GT}')
#         running = 0

#     if running >0:

# load the model
# list of images

# each images
# masks image
