RFSW_RP_DLM MATLAB Software Package, Beta version 1.0
by Yi-Zhong Wang
Retina Foundation of the Southwest
Copyright 2023. All Rights Reserved.

RFSW_RP_DLM is a software package used to generate the photoreceptor outer segment (OS) metrics measurement by deep learning model (DLM) reported in a longitudinal study of X-linked retinitis pigmentosa (RP) (ref). This beta version of the software is based on the work published previously:

Wang Y-Z, Galles D, Klein M, Locke, KG, Birch DG (2020) TVST 9.2.15
Wang Y-Z, Wu W, Birch DG. (2021) TVST 10.13.9
Wang Y-Z and Birch DG (2022) Front. Med. 9:932498

and consists of the following five components:

01_OCTScanFileGeneration
02_TrainingDatasetGeneration
03_DLMTraining
04_OCTScanClassification
05_OCTScanMeasurement

Accompanying this software package are an Excel file of OS metrics measurement results ("OS Metrics Measurement Results - XLRP Longitudinal Study - IOVS 2023.xlsx") and four examples of XML exports of two OCT volume scans, obtained using Spectralis SD-OCT from the left eye of a patient on 2 separate visits in the folder "Volume_Scan_Examples/XML_RP_Auto_HE_Manu_YZW_ID_P22". These XML exports are included for demonstrating the use of the software package. Each volume scan has two XLM exports: one with the layer boundary segmentation results by the Spectralis' built-in automatic segmentation software, and the other with manual correction or segmentation upon the automatic segmentation by Spectralis. The following are the instructions on how to run this software package.


01_OCTScanFileGeneration

Run the following script to convert XLM files to the OCTScan data structure used by the RFSW_RP_DLM software package:

script_sort_HE_XML_Export_to_OCTScan_multi_graders.m

Note that converting XML exports to OCTScan data structure requires a MATLAB Central File Exchange function "xml2struct" by Wouter Falkena (2023). Need to download this file to run the above script. (https://www.mathworks.com/matlabcentral/fileexchange/28518-xml2struct).

Select a folder containing XML exports, e.g., "XML_RP_Auto_HE_Manu_YZW_ID_P22", the script will generate OCTScan files in the specified output folder or in a new folder named "OCTScans_RP_Auto_HE_Manu_YZW_ID_P22" under the parent folder of XML export folder, if no output folder is specified.

For details of XML export folder naming convention as well as OCTScan structure, please refer to the notes in script_sort_HE_XML_Export_to_OCTScan_multi_graders.m


02_TrainingDatasetGeneration

