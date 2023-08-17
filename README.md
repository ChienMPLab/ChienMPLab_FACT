# FACT - Fast-and-Accurate Cell Tracking: a real-time cell segmentation and tracking algorithm to instantly export cellular characteristics from big image data

This repository contains the FACT segmentation and FACT tracking, an algorithm for real-time live-cell image segmentation and tracking, as described in the papers: 

Ting-Chun Chou*, Li You*, Cecile Beerens, Kate J. Feller, Jelle Storteboom, Miao-Ping Chien.
Fast-and-Accurate Cell Tracking (FACT): a real-time cell segmentation and tracking algorithm to instantly export cellular characteristics from big image data (*authors contributed equally)





## Workflow
### Prepare the training data for segmentation
Prepare the training dataset to train your own model. See the steps shown in the figure below.

![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/PrepareTrainingDataset.png)

### Train the segmentation model
After having the training dataset. Open train.ipynb from FACT_segmentation, and set up the path of the training dataset folder. Then, the user can start to train their own model. The workflow of training is shown in the figure.

![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/TrainingFlow.png)


### Cell tracking
After starting image acquirement and segmentation, set up the folder path for raw image and segmentation. Use main.m from FACT_Tracking for real-time cell tracking. The real-time cell tracking is shown in the following figure.

![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/Tracking.png)
