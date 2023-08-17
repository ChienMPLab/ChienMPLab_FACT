# FACT: a real-time cell segmentation and tracking algorithm to instantly export cellular characteristics from big image data

This repository contains the source codes for FACT segmentation and tracking, a real-time pipeline to instantly analyze large-scale of imaging data. Details can be found in:

Ting-Chun Chou*, Li You*, Cecile Beerens, Kate J. Feller, Jelle Storteboom, Miao-Ping Chien.
Fast-and-Accurate Cell Tracking (FACT): a real-time cell segmentation and tracking algorithm to instantly export cellular characteristics from big image data (*authors contributed equally)

Information on publication will be updated on time. 

## 1. Overview
Instant analysis of large imaging data is the advantage of using FACT: Once an image (e.g., ~10,000 cells) is generated, it is segmented (<5 secs) and tracked (< 1 min) right after. Quantification of cellular characteristics (such as cell movement, and lineage over 24 hrs) can be obtained soon after image acquisition. FACT consists of two interactive parts: GTWeka segmentation and tracking. 

FACT segmentation is called _GTWeka_: it is built based on [_Trainable Weka_](https://github.com/fiji/Trainable_Segmentation) or [_Illastik_](https://github.com/ilastik/ilastik) with ground-truth annotation as a 'golden standard' to optimize the Weka model (i.e., a random forest) to select necessary image features for the best segmentation (measured by IoU values). An optimized random forest is referred to as a GTWeka model, which can be used as a pre-trained cell segmentation model. 
![image](https://github.com/ChienMPLab/ChienMPLab_FACT/assets/42544588/52cb9954-eb3b-4561-8e75-27a9e4c973cd)

We release our pre-trained GTWeka in this repository (see below 'use-pretrained-GTWeka'). For the dataset that requires a different GTWeka model, we guide you on how to train your GTWeka (see below 'train-your-GTWeka'). 

FACT tracking combines [_Graussian Mixture Model Tracking_](https://sourceforge.net/projects/funseq/) with real-time cell track correction. It helps to reduce contamination caused by cell overlap. We release FACT tracking as a Matlab package in this repository 'FACT_Tracking'. Its Python version will be uploaded by the end of 2023.  

![image](https://github.com/ChienMPLab/ChienMPLab_FACT/assets/42544588/43e92ab9-1c8d-462c-86e6-11732ba02020)

## 2 Usage: FACT Segmentation - GTWeka
## 3 Usage: FACT Tracking
Prepare the training dataset to train your own model. See the steps shown in the figure below.

![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/PrepareTrainingDataset.png)

### Train the segmentation model
After having the training dataset. Open train.ipynb from FACT_segmentation, and set up the path of the training dataset folder. Then, the user can start to train their own model. The workflow of training is shown in the figure.

![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/TrainingFlow.png)


### Cell tracking
After starting image acquirement and segmentation, set up the folder path for raw image and segmentation. Use main.m from FACT_Tracking for real-time cell tracking. The real-time cell tracking is shown in the following figure.

![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/Tracking.png)
