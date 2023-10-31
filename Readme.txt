RFSW_RP_DLM MATLAB Software Package, Beta version 1.0
by Yi-Zhong Wang
Retina Foundation of the Southwest
Copyright 2023. All Rights Reserved.

RFSW_RP_DLM is a software package in MATLAB used to obtain photoreceptor outer segment (OS) metrics measurements by deep learning model (DLM) reported in a longitudinal study of X-linked retinitis pigmentosa (RP) (IOVS, 2023). This beta version of the software is based on the work published previously:

Wang Y-Z, Galles D, Klein M, Locke, KG, Birch DG (2020) TVST 9.2.15
Wang Y-Z, Wu W, Birch DG. (2021) TVST 10.13.9
Wang Y-Z and Birch DG (2022) Front. Med. 9:932498

and consists of the following parts:

01_OCTScanFileGeneration
02_DLMTrainingDatasetGeneration
03_DLMTraining
04_OCTScanClassification
05_OCTScanMeasurement
06_SharedFunctions
07_Utilities

Download the software package and add it to the search path of MATLAB. The software package was tested in MATLAB 2023b on both Apple silicon and Intel processor Macs running macOS Sonoma 14.1 or macOS Big Sur 11.7.10. In addition, accompanying this software package are an Excel file of OS metrics measurement results ("OS Metrics Measurement Results - XLRP Longitudinal Study - IOVS 2023.xlsx") and four examples of XML exports of two deidentified OCT volume scans, obtained using Heidelberg Spectralis SD-OCT, in the folder "Volume_Scan_Examples/XML_RP_Auto_HE_Manu_YZW_ID_P22". These XML exports are included for demonstrating the use of the software package. Each volume scan has two XLM exports: one with the layer boundary segmentation results by the Spectralis' built-in automatic segmentation software, and the other with manual correction upon the automatic segmentation by Spectralis for five boundary lines: inner limiting membrane (ILM), distal inner nuclear layer (dINL), ellipsoid zone (EZ), proximal retinal pigment epithelium (pRPE), and Bruch's membrane (BM). The following are the instructions on how to run each part of this software package.


##
01_OCTScanFileGeneration

Run the following script to convert XLM files exported from Spectralis to the OCTScan data structure used by the RFSW_RP_DLM software package:
	script_sort_HE_XML_Export_to_OCTScan_multi_graders.m
Note that converting XML exports to OCTScan data structure requires a MATLAB Central File Exchange function "xml2struct" by Wouter Falkena (2023). Need to download this file to run the above script. (https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct).

Select a folder containing XML exports, e.g., "XML_RP_Auto_HE_Manu_YZW_ID_P22", the script will generate OCTScan files in the designated output folder or in a new folder named "OCTScans_RP_Auto_HE_Manu_YZW_ID_P22" under the parent folder of XML export folder if no output folder is specified.

For details of XML export folder naming convention as well as OCTScan structure, please refer to the notes in "script_sort_HE_XML_Export_to_OCTScan_multi_graders.m"


##
02_DLMTraningDatasetGeneration

Run the following script to generate datasets for the training and validation of the sliding-window (SW) CNN models:
	script_generate_training_datasets_SW.m
and run the following script to generate datasets for the training and validation of U-Net CNN models:
	script_generate_training_datasets_UNet.m
The comments in the scripts contain details on datasets generation. OCTScan files obtained in "01_OCTScanFileGeneration" can be used to generate example datasets. 


##
03_DLMTraining

Run the following scripts to train SW models or UNet models:
	train_RP_model_SW_batch.m
	train_RP_model_UNet_batch.m
The comments in these two scripts contain details on the SW model and UNet model. Datasets generated with "02_DLMTrainingDatasetGeneration" can be used as examples for CNN training.


##
04_OCTScanClassification

Run the following script:
	script_classify_OCTScans_single_folder_paired_models.m
to classify OCTScans in a folder by specified trained DLMs (in pairs of UNet and the SW models). Refer to the comments in this script for details.


##
05_OCTScanMeasurement

Run the following script:
	script_get_layer_metrics_from_OCTScan_all_segmentations
to obtain retinal layer metrics. The default layer is PR1(EZ)-pRPE, i.e., OS layer. Refer to the comments in this script for details.

##
06_SharedFunctions

This part of the software package contains shared functions, as well as a local connected area searching algorithm (LCASA) for the SW model (Wang et al., 2020 TVST), as well as a modified version (LCASAP) for the hybrid model (Wang et al., 2021 TVST).


##
07_Utilities

Two utility scripts are included:

(i) script_check_OCTBScan_segmentation.m
		which displays layer boundary segmentations on a B-scan image from OCTScan.

(ii) script_OCTScan_Display_Layer_Metrics.m
		which plots OS thickness, EZ area, and OS volume from a volume OCTScan.

	