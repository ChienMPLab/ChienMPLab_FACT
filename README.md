# FACT - Fast-and-Accurate Cell Tracking: a real-time segmentation and tracking algorithm to instantly export cellular characteristics from large scale imaging

**Ting-Chun Chou<sup>1,2,3,4</sup>, Li You<sup>1,2,3,4</sup>, Cecile Beerens<sup>1,2,3</sup>, Kate J. Feller<sup>1,2,3</sup>, Jelle Storteboom<sup>1,2,3</sup>, Miao-Ping Chien<sup>1,2,3,5,6</sup>**

<sup><sup>1</sup>Department of Molecular Genetics, Erasmus University Medical Center, Rotterdam, The Netherlands

<sup><sup>2</sup>Erasmus MC Cancer Institute, The Netherlands.

<sup><sup>3</sup>Oncode Institute, Utrecht, The Netherlands.

<sup><sup>4</sup>These authors contributed equally

<sup><sup>5</sup>Senior author

<sup><sup>6</sup>Lead contact




## Workflow
### Prepare the training data for segmentation
Prepare the training dataset to train your own model. See the steps shown in the figure below.
![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/PrepareTrainingDataset.png)

### Train the segmentation model
After having the training dataset. Open train.ipynb from FACT_segmentation, and set up the path of the training dataset folder. Then, the user can start to train their own model. The workflow of training shown in the figure:
![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/TrainingFlow.png)


### Cell tracking
After starting image acquirement and segmentation, set up the folder path for raw image and segmentation. Use main.m from FACT_Tracking for real-time cell tracking. The real-time cell tracking shown in the following figure
![](https://github.com/ChienMPLab/ChienMPLab_FACT/blob/main/images/Tracking.png)
