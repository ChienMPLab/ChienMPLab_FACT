import os
from os import sep, path, mkdir
import numpy as np
from glob import glob
import cv2
from cv2 import imread, imwrite
import shutil
import pickle
import re
from cupy.cuda import Device
import cupy as cp
from evaluation import *
# from WSPL_ImgFeatures import *
from WSPL_ImageProcessing import *


from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier

from cuml.ensemble import RandomForestClassifier as curfc
from cucim.skimage.measure import label as culabel #,regionprops_table

from skimage.morphology import label
from skimage.color import label2rgb
from skimage.measure import regionprops #,regionprops_table, label
from skimage.segmentation import relabel_sequential

from IPython.display import clear_output
from scipy.stats import spearmanr

# from joblib import dump, load
import joblib
import json

def generate_TrainFolder( folder_path, imglist):
    folder_train = f'{folder_path}{sep}Train_data'
    if not path.exists( folder_train): mkdir( folder_train)
    folder_train_model = f'{folder_train}{sep}model'
    if not path.exists( folder_train_model): mkdir( folder_train_model)
    folder_train_images = f'{folder_train}{sep}images'
    if not path.exists( folder_train_images): mkdir( folder_train_images)
    folder_tmp = f'{folder_path}{sep}GT_Opt'
    if not path.exists( folder_tmp): mkdir( folder_tmp)
    folder_tmp_sub = f'{folder_tmp}{sep}Raw'
    if not path.exists( folder_tmp_sub): mkdir( folder_tmp_sub)
    folder_tmp_sub = f'{folder_tmp}{sep}label'
    if not path.exists( folder_tmp_sub): mkdir( folder_tmp_sub)
    folder_tmp = f'{folder_path}{sep}GT'
    if not path.exists( folder_tmp): mkdir( folder_tmp)


    for ipath in imglist:
        tmp_folder = f'{folder_train_images}{sep}{path.basename( ipath).split(".tif")[0]}'
        if not path.exists( tmp_folder): mkdir( tmp_folder)

        raw_img_folder = f'{tmp_folder}{sep}Raw'
        if not path.exists( raw_img_folder): mkdir( raw_img_folder)
        shutil.move(
            ipath,
            f'{raw_img_folder}{sep}{path.basename(ipath)}'
        )
        ann_folder = f'{tmp_folder}{sep}Annotation'
        if not path.exists( ann_folder): mkdir( ann_folder)

def saveNumpyArray( filename, array):
    with open(f'{filename}.pkl','wb+') as f: pickle.dump(array, f, protocol=4)

def convert_Annoation( folder_train):
    folder_train = f'{folder_train}{sep}images{sep}'
    train_folder_list =  sorted( os.listdir( folder_train) )
    for im, iforlder in enumerate( train_folder_list):
        sub_folder = f'{folder_train}{sep}{iforlder}'
        if path.isfile( sub_folder):continue
        annotation_folder = f'{sub_folder}{sep}Annotation'
        label_imgs = sorted( glob( f"{annotation_folder}{sep}*tif") )
        print( sub_folder)
        for i, i_imgs in enumerate(label_imgs):
            img = imread( i_imgs, -1)
            label_index = re.findall(
                r'_(\d+)'
                ,path.basename(i_imgs).split('.tif')[0]
                )[0]
            # print( label_index)
            label_loc = np.flatnonzero( img)
            # print( label_loc.size)
            if i < 1:
                imgtemp = np.zeros_like( img, dtype=int)
            imgtemp[ img> 0] = int( label_index)

        loc = np.flatnonzero( imgtemp)
        result = np.zeros( (loc.size,2), dtype=int)
        result[:, 0] = loc
        result[:, 1] = imgtemp.ravel()[loc]

        imwrite(
            f'{sub_folder}{sep}Label_img.tif',
            imgtemp
        )
        filename = f'{sub_folder}{sep}Annotation'
        saveNumpyArray( filename, result)

def getTrainingData( folder_train, sigmaarray):
    train_folder_list =  sorted( os.listdir( folder_train) )
    ii = 0
    for im, iforlder in enumerate( train_folder_list):
        sub_folder = f'{folder_train}{sep}{iforlder}'
        if path.isfile( sub_folder):continue

        annotation_label = np.load(
            f'{sub_folder}{sep}Annotation.pkl', allow_pickle=True
            )
        ytemp = np.zeros( (annotation_label.shape[0], 3), int)
        ytemp[:,0] = im
        ytemp[:,1:3] = annotation_label
        img_path = glob(f'{sub_folder}{sep}Raw{sep}*.tif')[0]
        imgarray= ggenFS( img_path, sigmaarray)
        imgarray = imgarray.reshape( (-1, imgarray.shape[-1]))

        if ii < 1:
            imgFS = imgarray[annotation_label[:,0],:]
            y = ytemp
        else:
            imgFS = np.concatenate(
                (imgFS,imgarray[annotation_label[:,0],:])
            )
            y = np.concatenate((y,ytemp))
        ii += 1
    return imgFS, y


def ggenFS( imgpath, sigmaarray):
    # st = time()
    img = imread( imgpath, -1).astype(np.float32)
    # img ratio
    img_ref = imgprocessing_gaussian( img, 300, 1)
    img_ratio = img/ img_ref
    # img features
    imgarray = imgprocessing( img_ratio, sigmaarray, 1)
    # end = time()
    # print( end-st)
    return imgarray

def getGTopt( GT_opt, sigmaarray):
    
    folder_raw = f'{GT_opt}{sep}Raw'
    img_list = sorted( glob( f'{folder_raw}{sep}*.tif') )
    folder_label = f'{GT_opt}{sep}label'
    xy_total = np.zeros( (len(img_list), 4), int)
    xytt = 0
    for im, img_path in enumerate( img_list):
        # img = imread( img_path, -1).astype(np.float32)
        imgarray= ggenFS( img_path, sigmaarray)
        xy_total[im, :2] = imgarray.shape[:2]
        xy_total[im, 2:4] = xytt, xytt+int(imgarray.shape[0]*imgarray.shape[1])
        xytt = xy_total[im, 3]
        # print(xy_total )
        imgarray = imgarray.reshape( (-1, imgarray.shape[-1]))

        imgpath = f'{folder_label}{sep}{path.basename( img_path)}'
        imgtmp = imread( imgpath, -1).flatten()

        if im < 1:
            imgFS = imgarray
            imgGT = imgtmp
        else:
            imgFS = np.concatenate(
                (imgFS,imgarray)
            )
            imgGT = np.concatenate(
                (imgGT,imgtmp)
            )
        
    return imgFS, xy_total, imgGT

def combin_paraterArray( test_n_setimators, test_max_features,\
    test_max_depth, test_smaples_leaf):
    parameters_array = np.zeros( 
        (int(test_n_setimators.size* len( test_max_features)\
        * test_max_depth.size* test_smaples_leaf.size), 4),
        dtype=int
        )
    iv = 0
    for i1 in test_n_setimators:
        for i2 in range(len(test_max_features)):
            for i3 in test_max_depth:
                for i4 in test_smaples_leaf:
                    parameters_array[iv, :] = i1, i2, i3, i4
                    iv +=1
    return parameters_array

def get_kFold_train_test_index( k, label_data_all):
    kf = KFold(n_splits=k, random_state=None, shuffle=False)
    # train_label =1 , test_label = 0
    train_test_index_result = np.zeros( 
        (label_data_all.shape[0], k),
        dtype=int
        )
    # get from each image
    image_ID = np.unique( label_data_all[:,0])
    for i, iD in enumerate( image_ID):
        ipage_index = np.flatnonzero( label_data_all[:,0] == iD )
        temp = label_data_all[ ipage_index,:]
        # print( temp.shape)
        ii_split = 0
        for train_index , test_index in kf.split(temp):
            train_test_index_result[ ipage_index[train_index], ii_split] = 1
            ii_split+=1
    return train_test_index_result



def opt_traing_parameters( select, parameters_array,test_max_features,\
    kf_num, imgFS, y, imgFS_GT, xy_total, imgGT):

    # model_sub_folder = f'{folder_model}{sep}Opt_{model_name}'
    # if not path.exists( model_sub_folder): mkdir( model_sub_folder)

    train_test_index_result = get_kFold_train_test_index( kf_num, y)
    train_accuracy = np.zeros( parameters_array.shape[0], dtype=float)

    for i_result, parameters in enumerate( parameters_array):
        nf = test_max_features[parameters[1]]
        cuml_model = curfc(
                n_estimators = parameters[0],
                max_features = nf,
                max_depth = parameters[2],
                min_samples_leaf = parameters[3]
        )
        kfold_tmp = np.zeros( kf_num, dtype = float)
        iik = 0
        for ik in range( train_test_index_result.shape[1]):
            train_index = np.flatnonzero( train_test_index_result[:,ik])
            test_index = np.flatnonzero( train_test_index_result[:,ik] < 1)
            cuml_model.fit(imgFS[train_index,:][:, select], y[train_index,:][:,2]-1)
            clear_output()
            kfold_tmp[iik] = model_accuracy_GTopt(imgFS_GT, xy_total, imgGT, cuml_model, select)
            iik +=1
        
        train_accuracy[i_result] = kfold_tmp.mean()

    idmax0 = train_accuracy.argmax()
    txt_output = f"Max: {train_accuracy.max()}\n"+\
        f"n_estimators:{ parameters_array[idmax0,0]}\n"+\
        f"max_features:{ test_max_features[parameters_array[idmax0,1]]}\n"+\
        f"max_depth:{ parameters_array[idmax0,2]}\n"+\
        f"min_samples_leaf:{ parameters_array[idmax0,3]}"
    result_list1 = ['n_estimators', 'max_features',\
        'max_depth', 'min_samples_leaf']
    result_list2 = [parameters_array[idmax0,0],test_max_features[parameters_array[idmax0,1]],\
        parameters_array[idmax0,2], parameters_array[idmax0,3]]
    
    return train_accuracy.max(), [result_list1, result_list2]
    # file_path = f"{model_sub_folder}{sep}Opt_{model_name}.txt"
    # with open( file_path, 'w+') as f: f.write( txt_output)
    # print( txt_output)


def getGT( folder_GT, sigmaarray):
    
    img_list = sorted( glob( f'{folder_GT}{sep}*{sep}GT{sep}Roi*.tif') )
    xy_total = np.zeros( (len(img_list), 4), int)
    xytt = 0
    for im, img_path in enumerate( img_list):
        imgtmp = imread( img_path, -1).flatten()

        imgpathImg = glob( f'{path.dirname(img_path)}{sep}TCC*.tif')[0]

        imgarray= ggenFS( imgpathImg, sigmaarray)
        xy_total[im, :2] = imgarray.shape[:2]
        xy_total[im, 2:4] = xytt, xytt+int(imgarray.shape[0]*imgarray.shape[1])
        xytt = xy_total[im, 3]
        # print(xy_total )
        imgarray = imgarray.reshape( (-1, imgarray.shape[-1]))


        if im < 1:
            imgFS = imgarray
            imgGT = imgtmp
        else:
            imgFS = np.concatenate(
                (imgFS,imgarray)
            )
            imgGT = np.concatenate(
                (imgGT,imgtmp)
            )
        
    return imgFS, xy_total, imgGT

def forwardSelection_opt(sub_folder_model, Top,\
    kf_num, imgFS, y, imgFS_GT, xy_total, imgGT):
    select = np.zeros( Top, int)
    if not path.exists( sub_folder_model): mkdir( sub_folder_model)
    numb_features = imgFS.shape[1]
    train_test_index_result = get_kFold_train_test_index( kf_num, y)
    cuml_model = curfc(
                n_estimators = 100
            )
    for i in range( Top):
        features_acc = np.zeros( numb_features, float)

        for i_f in range( numb_features):
            select[i] = i_f
            select_test = select[:i+1]
            uniqueSize = np.unique( select_test).size
            if select_test.size> uniqueSize: continue
            # with KFold
            for ik in range( train_test_index_result.shape[1]):
                train_index = np.flatnonzero( train_test_index_result[:,ik])
                cuml_model.fit(
                    imgFS[train_index,:][:, select_test],
                    y[train_index,:][:,2]-1
                    )
                clear_output()
                features_acc[i_f] += model_accuracy_GTopt(
                    imgFS_GT, xy_total, imgGT, cuml_model, select_test
                    )
            # wo kFold
            # cuml_model.fit(
            #         imgFS[:, select_test],
            #         y[:,2]-1
            #         )
            # clear_output()
            # features_acc[i_f] += model_accuracy_GTopt(
            #         imgFS_GT, xy_total, imgGT, cuml_model, select_test
            #         )
        select[i] = features_acc.argmax()
    np.savetxt(
        f'{sub_folder_model}{sep}select.txt',
        select
    )
    return select

def train_model(imgFS , y, modelInfo_tmp_path):
    with open(modelInfo_tmp_path, 'r') as openfile:
            modelInfo_tmp = json.load(openfile)
    if not'select' in modelInfo_tmp: select = np.arange( imgFS.shape[1])
    else: select = np.array(modelInfo_tmp['select'])
    
    pa_value = modelInfo_tmp['model_parameters_value']
    cuml_model = curfc(
        n_estimators= int(pa_value[0]),
        max_features=  pa_value[1],
        max_depth= int(pa_value[2]),
        min_samples_leaf = int(pa_value[3]),
        )
    cuml_model.fit(
            imgFS[:,select],
            y[:,2]-1
            )
    clear_output()
    model_path = modelInfo_tmp_path.replace('.json', '.model')
    joblib.dump( cuml_model, model_path)



def gpu_label( imgarray_test_y, gpuid):
    with Device(gpuid):
        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()
        
        imgarray_gpu = cp.asarray( imgarray_test_y)
        label_img = culabel(imgarray_gpu==0, connectivity=1)

        result_host = cp.asnumpy( label_img)

        del imgarray_gpu
        del label_img

        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()

        return result_host

def segmenation_GPU(X,  cuml_model, gpuid):
    with Device(gpuid):
        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()
        y = cuml_model.predict( X)
        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()
        return y


def model_accuracy_GTopt( imgFS, xy_total, imgGT, cuml_model, select):
    gt_predict = segmenation_GPU( imgFS[:, select],  cuml_model, 0)
    # gt_predict = imgGT
    accuracy_result = np.zeros( xy_total.shape[0])
    for i in range( xy_total.shape[0]):
        gt_predict_tmp = gt_predict[xy_total[i,2]:xy_total[i,3]].reshape(xy_total[i,:2])
        GT_data = imgGT[xy_total[i,2]:xy_total[i,3]].reshape(xy_total[i,:2])
        # generate the labeled Image
        label_img = gpu_label( gt_predict_tmp, 0)
        props = regionprops( label_img)
        # remove the small segmentation
        props_array = np.zeros( (len(props), 1), dtype=int) # area, x, y
        for ipp in range( len( props)): props_array[ ipp, 0] = props[ipp].area
        for i_label in np.flatnonzero( props_array[:,0] < 11):
            x, y = props[i_label].coords[:,1], props[i_label].coords[:,0]
            label_img[ y,x ] = 0

        ground_truth = label(GT_data)
        # Relabel objects
        ground_truth = relabel_sequential(ground_truth)[0] 
        prediction = relabel_sequential(label_img)[0]

        IOU = intersection_over_union(ground_truth, prediction)
        for t in np.arange(0.5, 1.0, 0.05):
            f1, tp, fp, fn, os, prec, rec = measures_at(t, IOU)
            accuracy_result[i] += os
        accuracy_result[i] = accuracy_result[i]/10
        # # # f1, tp, fp, fn, os1, prec, rec = measures_at(0.7, IOU)
        # accuracy_result[i] = os
        # for i_IOU in np.arange( 0.1,1,0.05):
        #     f1, tp, fp, fn, os, prec, rec = measures_at(i_IOU, IOU)
        #     accuracy_result[i] += os
    return accuracy_result.mean()

def model_eval( imgFS_select, xy_total, imgGT, cuml_model):
    results = pd.DataFrame(
        columns=["Image", "Threshold", "F1", "Jaccard", "TP",\
                "FP", "FN", "Official_Score", "Precision", "Recall"]
    )

    gt_predict = segmenation_GPU( imgFS_select,  cuml_model, 0)
    # gt_predict = imgGT
    accuracy_result = np.zeros( xy_total.shape[0])
    for i in range( xy_total.shape[0]):
        gt_predict_tmp = gt_predict[xy_total[i,2]:xy_total[i,3]].reshape(xy_total[i,:2])
        GT_data = imgGT[xy_total[i,2]:xy_total[i,3]].reshape(xy_total[i,:2])
        # generate the labeled Image
        label_img = gpu_label( gt_predict_tmp, 0)
        props = regionprops( label_img)
        # remove the small segmentation
        props_array = np.zeros( (len(props), 1), dtype=int) # area, x, y
        for ipp in range( len( props)): props_array[ ipp, 0] = props[ipp].area
        for i_label in np.flatnonzero( props_array[:,0] < 11):
            x, y = props[i_label].coords[:,1], props[i_label].coords[:,0]
            label_img[ y,x ] = 0

        ground_truth = label(GT_data)
        # Relabel objects
        ground_truth = relabel_sequential(ground_truth)[0] 
        prediction = relabel_sequential(label_img)[0]

        results = compute_af1_results(
                ground_truth, 
                prediction, 
                results, 
                i
            )
    return results

def features_importance( folder_model, imgFS, y, clf, kf_num, Top):
    if not path.exists( folder_model): mkdir( folder_model)
    train_test_index_result = get_kFold_train_test_index( kf_num, y)
    kfold_tmp = np.zeros( (imgFS.shape[1],kf_num), dtype = float)

    ii = 0
    for ik in range( train_test_index_result.shape[1]):
        train_index = np.flatnonzero( train_test_index_result[:,ik])
        test_index = np.flatnonzero( train_test_index_result[:,ik] < 1)
        # train the model
        clf.fit(imgFS[train_index,:], y[train_index,2]-1)
        kfold_tmp[:, ii] = clf.feature_importances_
    np.savetxt(
            f"{folder_model}{sep}feature_importances.txt",
            kfold_tmp,
            delimiter = ',', fmt = '%.16f'
        )

    corr = spearmanr(imgFS).correlation
    corr = np.nan_to_num(corr)
    # Ensure the correlation matrix is symmetric
    corr = (corr + corr.T) / 2
    np.fill_diagonal(corr, 1)
    distance_matrix = 1 - np.abs(corr)

    kfold_tmp_mean = kfold_tmp.mean( axis=1)
    rank = kfold_tmp_mean.argsort()[::-1]

    select = np.zeros( Top, int)
    i_select = 0
    select[i_select] = rank[0]
    i_select += 1
    for i, ir in enumerate( rank):
        test_cor = distance_matrix[ select[:i_select], ir ]
        if np.any( test_cor< 0.5): continue
        select[i_select] = ir
        i_select += 1
        if i_select >= 3: break
    for i, ir in enumerate( rank):
        test_cor = distance_matrix[ select[:i_select], ir ]
        if np.any( test_cor< 0.4): continue
        select[i_select] = ir
        i_select += 1
        if i_select >= Top: break
    
    print( select)
    np.savetxt(
        f'{folder_model}{sep}select.txt',
        select
    )
    return select