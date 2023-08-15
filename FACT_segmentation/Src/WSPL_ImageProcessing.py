from xmlrpc.client import INVALID_METHOD_PARAMS
import numpy as np
import cupy as cp
from cupy.cuda import Device
from numpy import exp


from cucim.skimage.filters import gaussian, sobel, difference_of_gaussians, median
from cucim.skimage.feature import hessian_matrix, hessian_matrix_eigvals
from cucim.skimage.morphology import disk
from cupyx.scipy.ndimage import convolve
from  scipy.ndimage import rotate


def mbj_kernel_host():
    ttarray = np.zeros( (19,19,30), dtype = np.float32)
    ttarray[:, 9, 0] = 1
    ir = 1
    for i in np.arange( 6, 180, 6):
        ttarray[:,:, ir] = rotate(ttarray[:,:, 0], i, reshape=False)
        ir +=1
    return ttarray

def hessian_matrix_gpu_input( img, sigma_value, result, mempool, pinned_mempool, gpuid):
    with Device(gpuid):
        H_elems = hessian_matrix(img, sigma=sigma_value, mode='reflect', cval=0, order='rc')
        eigenvalues = hessian_matrix_eigvals(H_elems)
        for i in range(2):
            result[:,:,i] = eigenvalues[i]

        ##Determinant
        result[:,:,2] = H_elems[0]*H_elems[2] - H_elems[1]**2
        ##Module
        result[:,:,3] =cp.sqrt(H_elems[0]**2+ H_elems[1]**2 + H_elems[2]**2 )
        ##Trace
        result[:,:,4] = H_elems[0]+H_elems[2]
        ##Square of Gamma-normalized eigenvalue difference
        result[:,:,5] = ((H_elems[0]-H_elems[2])**2+4*(H_elems[1])**2)*(H_elems[0]-H_elems[2])**2
        ##Orientation
        result[:,:,6]= 0.5*cp.arccos( cp.clip(result[:,:,5], -1, 1) )
        ##Gamma-normalized square eigenvalue difference
        result[:,:,7] = result[:,:,5]*(H_elems[0]-H_elems[2])**2


        del H_elems
        del eigenvalues
        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()

def membrane_projections( img, result, temp, mbj_kernel, mempool, pinned_mempool, gpuid):
    with Device(gpuid):
        for i in range( mbj_kernel.shape[2]):
            temp[:,:, i] = convolve(img, mbj_kernel[:,:,i], mode='reflect')

        # sum
        result[:,:, 0] = temp.sum( axis = 2)
        # mean
        result[:,:, 1] = temp.mean( axis = 2)
        # standard deviation
        result[:,:, 2] = temp.std( axis = 2)
        # median
        result[:,:, 3] = cp.median( temp, axis=2)
        # maximum
        result[:,:, 4] = temp.max( axis = 2)
        # minimum
        result[:,:, 5] = temp.min( axis = 2)

        del mbj_kernel
        del temp

        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()

def imgprocessing( img, sigmaarray, gpuid):
    sigaSize =  sigmaarray.size
    datatype = img.dtype
    method_ = np.zeros(6,np.int8)
    method_[0] = 1+sigaSize
    method_[1] = method_[0] + sigaSize
    method_[2] = method_[1] + int((sigaSize)*(sigaSize-1)/2)
    method_[3] = method_[2] + sigaSize* 8
    method_[4] = method_[3] +  6
    method_[5] = method_[4] + sigaSize
    imgyx = img.shape[:2]
    
    with Device(gpuid):
        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()

        result = cp.zeros( (imgyx[0],imgyx[1], method_[-1]), dtype=datatype)
        result[ :, :, 0] = cp.array( img)
        
        # gaussian
        ii = 0
        for i in range(  1,method_[0]):
            result[ :, :, i] = gaussian(
                result[ :, :, 0], sigma=sigmaarray[ii], mode='reflect'
                )
            ii += 1
        # sobel
        ii = 0
        for i in range(  method_[0],method_[1]):
            result[ :, :, i] = sobel(
                result[ :, :, i-sigaSize],
                mode='reflect'
                )
            ii += 1
        # difference_of_gaussians
        ii = method_[1]
        for i in range(sigaSize):
            low_sigma = sigmaarray[i]
            for j in range(i+1, sigaSize):
                high_sigma = sigmaarray[j]
                result[:, :, ii] = difference_of_gaussians(
                    result[ :, :, 0], low_sigma, high_sigma, mode='reflect'
                    )
                ii += 1
        # Hessian matrix
        ii = 0
        for i in range(  method_[2],method_[3], 8):
            hessian_matrix_gpu_input( 
                result[ :, :, 0],
                sigmaarray[ii],
                result[ :, :, i:i+8],
                mempool, pinned_mempool,
                gpuid
                )
            ii += 1
        # Membrane projections
        ii = method_[3]
        mbj_kernel = cp.array( mbj_kernel_host())
        temp = cp.zeros( (imgyx[0],imgyx[1], 30), dtype=datatype)
        membrane_projections( 
            result[ :, :, 0],
            result[ :, :, ii:ii+6], temp,
            mbj_kernel,
            mempool, pinned_mempool,
            gpuid
            )
        # Median filter
        ii = 0
        for i in range( method_[4], method_[5]):
            kernel = disk( sigmaarray[ii])
            result[ :, :, i] = median(
                result[ :, :, 0], kernel, mode='reflect'
                )
            ii += 1

        result_host = cp.asnumpy( result)

        del mbj_kernel
        del temp
        del kernel
        del result

        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()
        return result_host


def img_processing_method_template( sigmaarray):
    sigaSize =  sigmaarray.size
    method_ = np.zeros(6,np.int8)
    method_[0] = 1+sigaSize
    method_[1] = method_[0] + sigaSize
    method_[2] = method_[1] + int((sigaSize)*(sigaSize-1)/2)
    method_[3] = method_[2] + sigaSize* 8
    method_[4] = method_[3] +  6
    method_[5] = method_[4] + sigaSize

    result = np.zeros( method_[-1], dtype='U301')
    result[0] = 'Image ratio'
    ii = 0
    for i in range(  1,method_[0]):
        result[ i] = f'Gaussain {sigmaarray[ii]}'
        ii += 1
    ii = 0
    for i in range(  method_[0],method_[1]):
        result[ i] = result[ i-sigaSize]+ f'Sobel {sigmaarray[ii]}'
        ii += 1

    ik = method_[1]
    for i in range(sigaSize):
        low_sigma = sigmaarray[i]
        for j in range( i+1, sigaSize):
            high_sigma = sigmaarray[j]
            result[ ik] = f'Difference Gaussian {low_sigma} {high_sigma}'
            ik += 1

    hm_= [
            "ei0", "ei1", "Determinant", "Module", "Trace",
            "Square of Gamma-normalized eigenvalue difference",
            "Orientation", "Gamma-normalized square eigenvalue difference"
            ]
    isigma = 0
    for i in range(  method_[2],method_[3], 8):
        ii = 0
        for i_hm in hm_:
            result[i+ii] = f'hessian martix {sigmaarray[isigma]} {i_hm} '
            ii += 1
        isigma += 1

    mp = [
            "sum", "mean", "SD", "median", "max", "min"
        ]

    i = method_[3]
    ii= 0
    for imp in mp:
        result[i+ii] = f'Membrane_projections {imp} '
        ii += 1

    ii = 0
    for i in range(method_[4], method_[5]):
        result[i] = f'MedianFilter {sigmaarray[ii]} '
        ii += 1

    return result


def imgprocessing_sub(sigmaarray, featuresID):
    sigaSize =  sigmaarray.size
    method_ = np.zeros(6,np.int8)
    method_[0] = 1+sigaSize
    method_[1] = method_[0] + sigaSize
    method_[2] = method_[1] + int((sigaSize)*(sigaSize-1)/2)
    method_[3] = method_[2] + sigaSize* 8
    method_[4] = method_[3] +  6
    method_[5] = method_[4] + sigaSize

    parameters_DifferenceGaussian = np.zeros(
            (int(sigaSize*(sigaSize-1)/2),2),\
            dtype=int
            )
    parameters_HessianMatrix = np.zeros(
            (sigaSize*8,2),\
            dtype=int
            )
    ii =0 
    for i in range(0,sigaSize*8,8):
        parameters_HessianMatrix[i:i+8,0] = ii
        parameters_HessianMatrix[i:i+8,1] = np.arange(8)
        ii += 1
    ik = 0
    for i in range(sigaSize):
        low_sigma = sigmaarray[i]
        for j in range( i+1,sigaSize):
            high_sigma = sigmaarray[j]
            parameters_DifferenceGaussian[ik,:2] = np.array([low_sigma, high_sigma])
            ik += 1

    sub_processing = np.zeros((featuresID.size,2), int)
    for ii, i_ID in enumerate(featuresID): 
        
        if i_ID <1:
            # image ratio [img/img_ref]
            i_method = 0
            i_sigma = 0

        elif np.logical_and(1<= i_ID, i_ID< method_[0]):
            # Gaussain filter
            i_method = 1
            i_sigma = sigmaarray[i_ID-1]

        elif np.logical_and(method_[0]<= i_ID, i_ID< method_[1]):
            # Sobel filter
            i_method = 2
            i_sigma = sigmaarray[ i_ID-method_[0]]
        
        elif np.logical_and(method_[1]<= i_ID, i_ID< method_[2]):
            # Difference of gaussians
            i_method = 3
            i_sigma = i_ID-method_[1]

        elif np.logical_and(method_[2]<= i_ID, i_ID< method_[3]):
            # Hessian Matrix
            i_method = 4
            i_sigma = i_ID-method_[2]

        elif np.logical_and(method_[3]<= i_ID, i_ID< method_[4]):
            # MBJ
            i_method = 5
            i_sigma = i_ID-method_[3]

        elif np.logical_and(method_[4]<= i_ID, i_ID< method_[5]):
            # Median filter
            i_method = 6
            i_sigma = sigmaarray[i_ID-method_[4]]
        
        sub_processing[ii,0] = i_method
        sub_processing[ii,1] = i_sigma

    return sub_processing, parameters_DifferenceGaussian, parameters_HessianMatrix


def imgprocessing_subprocessing( img, gpuid, imgprocessing_sub_Info):
    sub_processing, parameters_DifferenceGaussian, parameters_HessianMatrix=  imgprocessing_sub_Info
    
    imgyx = img.shape[:2]
    datatype = img.dtype
    methods = np.unique( sub_processing[:,0])
    
    with Device(gpuid):
        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()

        result = cp.zeros( (imgyx[0],imgyx[1], sub_processing.shape[0]), dtype=datatype)
        img_device = cp.array( img)
        for i in range( sub_processing.shape[0]):
            i_method  = sub_processing[i,0]
            i_sigma = sub_processing[i,1]
            if i_method < 1:
                result[:,:,i] = img_device
            elif i_method  == 1:
                result[:,:,i] = gaussian(
                    img_device,
                    sigma=i_sigma, mode='reflect'
                    )
            elif i_method  == 2:
                result[:,:,i] = gaussian(
                    img_device,
                    sigma=i_sigma, mode='reflect'
                    )
                result[:,:,i] = sobel(
                    result[:,:,i],
                    mode='reflect'
                    )
            elif i_method == 3:
                low_sigma = parameters_DifferenceGaussian[i_sigma,0]
                high_sigma = parameters_DifferenceGaussian[i_sigma,1]
                print( f'{low_sigma}\t {high_sigma}')
                result[:,:,i] = difference_of_gaussians(
                    img_device, low_sigma, high_sigma, mode='reflect'
                    )
                # result[:,:,i] = img_device
                # result[:,:,i] = difference_of_gaussians(
                #     img_device, 2, 4, mode='reflect'
                #     )
            elif i_method == 4:
                result_Hessian = cp.zeros( (imgyx[0],imgyx[1],8), dtype=datatype)
                sigma_hessianMatrix = parameters_HessianMatrix[i_sigma,0]
                which_hessianMatrix = parameters_HessianMatrix[i_sigma,1]
                hessian_matrix_gpu_input( 
                    img_device,
                    sigma_hessianMatrix,
                    result_Hessian,
                    mempool, pinned_mempool,
                    gpuid
                    )
                result[:,:,i] = result_Hessian[:,:,which_hessianMatrix]
            elif i_method == 5:
                mbj_kernel = cp.array( mbj_kernel_host())
                temp = cp.zeros( (imgyx[0],imgyx[1], 30), dtype=datatype)
                result_mbp = cp.zeros( (imgyx[0],imgyx[1], 6), dtype=datatype)
                membrane_projections( 
                    img_device,
                    result_mbp, temp,
                    mbj_kernel,
                    mempool, pinned_mempool,
                    gpuid
                    )
                result[:,:,i] = result_mbp[:,:,i_sigma]
            
            elif i_method == 6:
                kernel = disk( i_sigma)
                result[ :, :, i] = median(
                    img_device, kernel, mode='reflect'
                    )
        
        result_host = cp.asnumpy( result)

        del img_device
        if np.any(methods==4): del result_Hessian
        if np.any(methods==5):
            del mbj_kernel
            del temp
            del result_mbp
        if np.any(methods==6): del kernel


        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()
        return result_host

        

def imgprocessing_gaussian( img, sigma, gpuid):
    with Device(gpuid):
        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()


        datatype = img.dtype
        size_result = 1
        x = img.shape[1]
        y = img.shape[0]
        # gaussian
        result = cp.asarray( img)
        result =  gaussian(
                result, sigma=sigma, mode='reflect'
                )
        result_host = cp.asnumpy( result)
        del result

        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()
        
    
        return result_host

def subFeatures_temeplate( sigmaarray):
    sigaSize =  sigmaarray.size
    # method_ = np.zeros( 7,np.int8)
    ##0 image
    ##1 gaussian
    ##2 sobel
    ##3 difference of gaussians
    ##4 Hessian matrix
    ##5 Membrane projections
    ##6 Median filter
    # method_[0] = 1+sigaSize # img_ratio + gaussian with various sigma
    # method_[1] = method_[0] + sigaSize
    # method_[2] = method_[1] + int((sigaSize)*(sigaSize-1)/2)
    # method_[3] = method_[2] + sigaSize* 8
    # method_[4] = method_[3] +  6
    # method_[5] = method_[4] + sigaSize

    # (gaussian, sobel, median)+ Hessian matrix+ differenceOfGaussians+ image+ membrane projections
    # amount_features = sigaSize* 3+ sigaSize* 8+ int((sigaSize)*(sigaSize-1)/2)+ 1+ 6

    # print( method_)
    method_ = [
        'Image ratio', 'Gaussian', 'Sobel', 'Difference of Gaussain', 'Hessian Matrix',
        'Membrane projections', 'Median filter'
    ]
    method_size = np.array([ 1, sigaSize, sigaSize, int((sigaSize)*(sigaSize-1)/2), sigaSize* 8, 6, sigaSize])
    submethod_HM = [
        "ei0", "ei1", "Determinant", "Module", "Trace",
        "Square of Gamma-normalized eigenvalue difference",
        "Orientation", "Gamma-normalized square eigenvalue difference"
    ]
    submethod_MP = ["sum", "mean", "SD", "median", "max", "min"]
    result = np.zeros( method_size.sum(), dtype='U301')
    parameters_ = {} # [method, [parameteres]]
    # print( result.size)
    ii = 0
    for i, imethod in enumerate(method_):
        if i <1:
            # image ratio
            result[ii] = imethod
            parameters_[ii] = [0, [0]]
            ii += 1
        elif i in {1,2,6} :
            # gaussian, sobel, median
            for i_sigma in sigmaarray:
                result[ii] = f'{imethod} {i_sigma}'
                parameters_[ii] = [i, [i_sigma]]
                ii += 1
        elif i == 3:
            # difference of gaussain
            for j in range(sigaSize):
                low_sigma = sigmaarray[j]
                for jj in range( j+1, sigaSize):
                    high_sigma = sigmaarray[jj]
                    result[ii] = f'{imethod} {low_sigma} {high_sigma}'
                    parameters_[ii] = [i, [low_sigma, high_sigma]]
                    ii += 1
        elif i == 4:
            # Hessian Matrix
            for i_sigma in sigmaarray:
                for ii_, i_hm in enumerate(submethod_HM):
                    result[ii] = f'{imethod}  {i_sigma} {i_hm} '
                    parameters_[ii] = [i, [i_sigma, ii_]]
                    ii += 1
        elif i == 5:
            # Membrane projections
            for ii_, i_mp in enumerate(submethod_MP):
                result[ii] = f'{imethod}  {i_mp}'
                parameters_[ii] = [i, [ii_]]
                ii += 1

    return result, parameters_

def imgprocessing_subprocessingV2( img, gpuid, parameters_dic, id_subFeatures):
    imgyx = img.shape[:2]
    datatype = img.dtype
    
    with Device(gpuid):
        mempool = cp.get_default_memory_pool()
        pinned_mempool = cp.get_default_pinned_memory_pool()

        result = cp.zeros( (imgyx[0],imgyx[1], id_subFeatures.shape[0]), dtype=datatype)
        img_device = cp.array( img)
        total_method = np.zeros(id_subFeatures.size, int)
        for i, i_id_subFeatures in enumerate( id_subFeatures):
            
            i_method, i_parameters = parameters_dic[i_id_subFeatures]
            total_method[i] = i_method

            if i_method < 1:
                # image ratio
                result[:,:,i] = img_device
            elif i_method == 1:
                # gaussian
                i_sigma = i_parameters[0]
                result[:,:,i] = gaussian(
                    img_device,
                    sigma=i_sigma, mode='reflect'
                )
            elif i_method == 2:
                # sobel
                i_sigma = i_parameters[0]
                result[:,:,i] = gaussian(
                    img_device,
                    sigma=i_sigma, mode='reflect'
                )
                result[:,:,i] = sobel(
                    result[:,:,i],
                    mode='reflect'
                )
            elif i_method == 6:
                # median
                i_sigma = i_parameters[0]
                kernel = disk( i_sigma)
                result[ :, :, i] = median(
                    img_device, kernel, mode='reflect'
                    )
            elif i_method == 3:
                # difference of gaussians
                low_sigma, high_sigma = i_parameters
                result[:,:,i] = difference_of_gaussians(
                    img_device, low_sigma, high_sigma, mode='reflect'
                )
            elif i_method == 4:
                # hessian matrix
                result_Hessian = cp.zeros( (imgyx[0],imgyx[1],8), dtype=datatype)
                sigma_hessianMatrix, which_hessianMatrix = i_parameters
                hessian_matrix_gpu_input( 
                    img_device,
                    sigma_hessianMatrix,
                    result_Hessian,
                    mempool, pinned_mempool,
                    gpuid
                )
                result[:,:,i] = result_Hessian[:,:,which_hessianMatrix]
            elif i_method == 5:
                # membrane projections
                i_sigma = i_parameters[0]
                mbj_kernel = cp.array( mbj_kernel_host())
                temp = cp.zeros( (imgyx[0],imgyx[1], 30), dtype=datatype)
                result_mbp = cp.zeros( (imgyx[0],imgyx[1], 6), dtype=datatype)
                membrane_projections( 
                    img_device,
                    result_mbp, temp,
                    mbj_kernel,
                    mempool, pinned_mempool,
                    gpuid
                    )
                result[:,:,i] = result_mbp[:,:,i_sigma]

        result_host = cp.asnumpy( result)
        # empty the memory on gpu
        del img_device
        del result

        if np.any(total_method==4): del result_Hessian
        if np.any(total_method==5):
            del mbj_kernel
            del temp
            del result_mbp
        if np.any(total_method==6): del kernel

        mempool.free_all_blocks()
        pinned_mempool.free_all_blocks()
        return result_host