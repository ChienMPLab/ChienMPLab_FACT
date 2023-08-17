# FACT: a real-time cell segmentation and tracking algorithm to instantly export cellular characteristics from big image data

This repository contains the source codes for FACT segmentation and tracking, a real-time pipeline to instantly analyze large-scale of imaging data. Details can be found in:

Ting-Chun Chou*, Li You*, Cecile Beerens, Kate J. Feller, Jelle Storteboom, Miao-Ping Chien.
Fast-and-Accurate Cell Tracking (FACT): a real-time cell segmentation and tracking algorithm to instantly export cellular characteristics from big image data (*authors contributed equally)

**Information on publication will be updated on time.**

## 1. Overview
Instant analysis of large imaging data is the advantage of using FACT: Once an image (e.g., ~10,000 cells) is generated, it is segmented (<5 secs) and tracked (<1 min) right after. Quantification of cellular characteristics (such as cell movement and lineage over 24 hrs) can be obtained soon after image acquisition. FACT consists of two interactive parts: GTWeka segmentation and tracking. 

FACT segmentation is called GTWeka: it is built based on [Trainable Weka](https://github.com/fiji/Trainable_Segmentation) or [Illastik](https://github.com/ilastik/ilastik) with ground-truth annotation as a 'golden standard' to optimize the Weka model (i.e., a random forest) to select necessary image features for the best segmentation (measured by IoU values). An optimized random forest is referred to as a GTWeka model, which can be used as a pre-trained cell segmentation model. 
![image](https://github.com/ChienMPLab/ChienMPLab_FACT/assets/42544588/52cb9954-eb3b-4561-8e75-27a9e4c973cd)

We release our pre-trained GTWeka in this repository (see below 'Use pre-trained GTWeka'). For the dataset that requires a different GTWeka model, we guide you on how to train your GTWeka (see below 'Train your GTWeka'). 

FACT tracking combines [Gaussian Mixture Model Tracking](https://sourceforge.net/projects/funseq/) with real-time cell track correction. It helps to reduce contamination caused by cell overlap. We release FACT tracking as a Matlab package in this repository 'FACT_Tracking'. Its Python version will be uploaded by the end of 2023.  

![image](https://github.com/ChienMPLab/ChienMPLab_FACT/assets/42544588/43e92ab9-1c8d-462c-86e6-11732ba02020)

## 2. Usage: FACT Segmentation - GTWeka
### 2.1 Use pre-trained GTWeka
**@Jason **
### 2.2 Train your GTWeka
### Annotation 
We suggest annotation with Labkit:
1. Install [Fiji](https://fiji.sc/) and [Labkit](https://imagej.net/plugins/labkit/).
2. Load images with Labkit: Plugins > Labkit > Open Current Image With Labkit.
3. Add training data for each labeling: foreground, background, and edge, then click 'Pixel Classification'.
   
   ![image](https://github.com/ChienMPLab/ChienMPLab_FACT/assets/42544588/e9258eda-ce25-4f16-9ad3-47655867e18a)

4. Generate foreground (i.e., pixels that are not background) segmentation from Labkit: Segmentation > Create label from segmentation > foreground. Here the segmentation is only binary. 
  
   ![image](https://github.com/ChienMPLab/ChienMPLab_FACT/assets/42544588/b5111081-53c9-42ad-8854-aac81d053d09)

5. Convert the outcome of 4 to cell instances: Plugin > BIOP > Image Analysis > ROIs > Label image to ROIs.
6. Maunual check on cell instances. For example, if Cell15 (in blue) is not nicely separated, one can re-do the cell contour and correct the segmentation.
   
   ![image](https://github.com/ChienMPLab/ChienMPLab_FACT/assets/42544588/83cca740-fda8-41b2-aa99-3b2f19b839dd)

7. Generate final segmentation and save it to a 16-bit image: Plugins > LOCI > ROI Map
   
   ![image](https://github.com/ChienMPLab/ChienMPLab_FACT/assets/42544588/404318c7-57a9-46bb-ae0b-6f5405be2abf)


### Training
When the training dataset is ready:
1. load 'train.ipynb' from 'FACT_Segmentation', and set up the path of the training dataset folder.
2. train your model, while successively adding annotations (via Labkit, see above) if higher performance is expected. 

## 3. Usage: FACT Tracking
'FACT_Tracking' is a complete package ready to run. It is scripted in Matlab. **We will upload a Python version in the near future**. In this package: 
1. 'FACT_Tracking/Data/Raw' contains all raw images.
2. 'FACT_Tracking/Data/GTWeka' contains all cell masks.
3. 'FACT_Tracking/Data/Input' contains input images for tracking, where Input = Raw .* (GTWeka > 0).
4. 'set_paths.m' let you edit folder paths of raw images, input images etc. It sets also default parameters such as gap-closing-distance, gap-closing-window, etc. 
5. 'main.m' is the function to run tracking. 
6. The rest .m files are helper functions, most of the time, you do not need to edit them.
7. The tracking outcome will be saved as XML files. One can visualize tracks using Fiji and [Mamut](https://imagej.net/plugins/mamut/) Plugin. 

## 4. Tracking Visualization
