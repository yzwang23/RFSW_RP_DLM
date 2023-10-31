% script_get_layer_metrics_summary_all_segmentations.m
%   script to obtain mean layer metrics of the measurements for various
%   segmentation methods and construct a summary array for output. The data
%   can be obtained from the patients with longitudinal follow-up visits. 
%   Different patients' OCT scans were stored in separate folders and layer
%   metrics by various segment methods are already obtained and saved by
%   using "script_get_layer_metrics_from_OCTScan_all_segmentations.m"
%
% Yi-Zhong Wang, 10/13/2022
% Retina Foundation of the Southwest
% Copyright. All Rights Reserved.
%

clearvars; clc; close all;

% make sure the current directory is accessible so that uigetdir works
cd(fileparts(mfilename('fullpath')));

clearvars; clc; close all;
warning('off');     % to suppress "Dot indexing is not supported ..."

% make sure the current directory is accessible so that uigetdir works
cd(fileparts(mfilename('fullpath')));


%**************************************************************************
disp('specify a directory containing OCTScan*.mat files or folders ...');
%**************************************************************************
DirFullName = uigetdir;
if DirFullName==0
  disp('OCTScan folder not specified ...');
  return;
end
DirPathName = [DirFullName, filesep];

% get rid of entries with leading '.' in the name
% note that 'dir' command lists names starting with '.' first
D_all = dir(DirFullName);
count = 0;
for i=1:length(D_all)
  if strcmp(D_all(i).name(1), '.')
    count = count + 1;
  end
end
D = D_all(count+1:length(D_all));       % exclude '.' & '..' directory

%**************************************************************************
% init output variables for all patients
%**************************************************************************
PtIDAll                      = [];
SummaryLayerMetricsAutoAll   = [];
SummaryLayerMetricsManuAll   = [];
SummaryLayerMetricsUNetAll   = [];
SummaryLayerMetricsHybridAll = [];
% SummaryThickAutoAll          = [];
% SummaryAreaAutoAll           = [];
% SummaryVolAutoAll            = [];
SummaryThickManuAll          = [];
SummaryAreaManuAll           = [];
SummaryVolManuAll            = [];
SummaryThickUNetAll          = [];
SummaryAreaUNetAll           = [];
SummaryVolUNetAll            = [];
SummaryThickHybridAll        = [];
SummaryAreaHybridAll         = [];
SummaryVolHybridAll          = [];

% ADDED on 07/24/2023
SummaryCTimeUNetAll          = [];
SummaryCTimeHybridAll        = [];


%**************************************************************************
% perform analyses, need to determine a single folder or not
%**************************************************************************
if D(1).isdir==0 && contains(D(1).name, 'OCTScan_') && contains(D(1).name, '.mat')
  % 1st item in D is not a folder & is an OCTScan*.mat file -> a single folder with a list of OCTScan*.mat files
  disp('OCTScan*.mat files are in the selected folder. Get layer measurements from Results folder ...');
  fprintf('\n');

  % get summary results of layer metrics measurement
  DirListToCheck = D;
  SummaryR = GetLayerMetricsSummaryResult(DirListToCheck);

  % summary results
  PtIDAll                      = [PtIDAll; SummaryR.PtID];
  SummaryLayerMetricsAutoAll   = [SummaryLayerMetricsAutoAll; SummaryR.SummaryLayerMetricsAuto];
  SummaryLayerMetricsManuAll   = [SummaryLayerMetricsManuAll; SummaryR.SummaryLayerMetricsManu];
  SummaryLayerMetricsUNetAll   = [SummaryLayerMetricsUNetAll; SummaryR.SummaryLayerMetricsUNet];
  SummaryLayerMetricsHybridAll = [SummaryLayerMetricsHybridAll; SummaryR.SummaryLayerMetricsHybrid];

  % MOD on 07/06/2023 to get separate summary for individual metrics
  % with mean, std, as well as individual measurements
%   SummaryThickAutoAll   = [SummaryThickAutoAll; SummaryR.SummaryThickAuto];
%   SummaryAreaAutoAll    = [SummaryAreaAutoAll; SummaryR.SummaryAreaAuto];
%   SummaryVolAutoAll     = [SummaryVolAutoAll; SummaryR.SummaryVolAuto];

  SummaryThickManuAll   = [SummaryThickManuAll; SummaryR.SummaryThickManu];
  SummaryAreaManuAll    = [SummaryAreaManuAll; SummaryR.SummaryAreaManu];
  SummaryVolManuAll     = [SummaryVolManuAll; SummaryR.SummaryVolManu];

  SummaryThickUNetAll   = [SummaryThickUNetAll; SummaryR.SummaryThickUNet];
  SummaryAreaUNetAll    = [SummaryAreaUNetAll; SummaryR.SummaryAreaUNet];
  SummaryVolUNetAll     = [SummaryVolUNetAll; SummaryR.SummaryVolUNet];

  SummaryThickHybridAll = [SummaryThickHybridAll; SummaryR.SummaryThickHybrid];
  SummaryAreaHybridAll  = [SummaryAreaHybridAll; SummaryR.SummaryAreaHybrid];
  SummaryVolHybridAll   = [SummaryVolHybridAll; SummaryR.SummaryVolHybrid];
  % end of MOD

  % ADDED on 07/24/2023
  SummaryCTimeUNetAll   = [SummaryCTimeUNetAll; SummaryR.SummaryCTimeUNet];
  SummaryCTimeHybridAll = [SummaryCTimeHybridAll; SummaryR.SummaryCTimeHybrid];
  % end
    
else
  disp('OCTScan*.mat files are in the subfolders under the selected folder ...');
  disp('  loop through all subfolders to get layer measurements ...');
  fprintf('\n');

  % remove non-folder items from D
  DirFlg = zeros(length(D),1);
  for j=1:length(DirFlg)
    if D(j).isdir
      DirFlg(j) = 1;
    end
  end
  D_temp = D((DirFlg>0));
  D = D_temp;

  numSubfolder = length(D);
  for l=1:numSubfolder
    disp(['***Sub-folder name: ', D(l).name]);

    DirListToCheck = dir(fullfile(D(l).folder, D(l).name));

    % get rid of entries with leading '.' in the name
    count = 0;
    for i=1:length(DirListToCheck)
      if strcmp(DirListToCheck(i).name(1), '.')
        count = count + 1;
      end
    end
    DirListToCheck = DirListToCheck(count+1:length(DirListToCheck));       % exclude '.' & '..' directory

    % get summary results of layer metrics measurement - the following
    % function was hard-coded to handle 3 DL models -- need to be revised
    SummaryR = GetLayerMetricsSummaryResult(DirListToCheck);

    % append all summary results
    PtIDAll                      = [PtIDAll; SummaryR.PtID];
    SummaryLayerMetricsAutoAll   = [SummaryLayerMetricsAutoAll; SummaryR.SummaryLayerMetricsAuto];
    SummaryLayerMetricsManuAll   = [SummaryLayerMetricsManuAll; SummaryR.SummaryLayerMetricsManu];
    SummaryLayerMetricsUNetAll   = [SummaryLayerMetricsUNetAll; SummaryR.SummaryLayerMetricsUNet];
    SummaryLayerMetricsHybridAll = [SummaryLayerMetricsHybridAll; SummaryR.SummaryLayerMetricsHybrid];

    % MOD on 07/06/2023 to get separate summary for individual metrics
    % with mean, std, as well as individual measurements
%     SummaryThickAutoAll   = [SummaryThickAutoAll; SummaryR.SummaryThickAuto];
%     SummaryAreaAutoAll    = [SummaryAreaAutoAll; SummaryR.SummaryAreaAuto];
%     SummaryVolAutoAll     = [SummaryVolAutoAll; SummaryR.SummaryVolAuto];

    SummaryThickManuAll   = [SummaryThickManuAll; SummaryR.SummaryThickManu];
    SummaryAreaManuAll    = [SummaryAreaManuAll; SummaryR.SummaryAreaManu];
    SummaryVolManuAll     = [SummaryVolManuAll; SummaryR.SummaryVolManu];

    SummaryThickUNetAll   = [SummaryThickUNetAll; SummaryR.SummaryThickUNet];
    SummaryAreaUNetAll    = [SummaryAreaUNetAll; SummaryR.SummaryAreaUNet];
    SummaryVolUNetAll     = [SummaryVolUNetAll; SummaryR.SummaryVolUNet];

    SummaryThickHybridAll = [SummaryThickHybridAll; SummaryR.SummaryThickHybrid];
    SummaryAreaHybridAll  = [SummaryAreaHybridAll; SummaryR.SummaryAreaHybrid];
    SummaryVolHybridAll   = [SummaryVolHybridAll; SummaryR.SummaryVolHybrid];
    % end of MOD

    % ADDED on 07/24/2023
    SummaryCTimeUNetAll   = [SummaryCTimeUNetAll; SummaryR.SummaryCTimeUNet];
    SummaryCTimeHybridAll = [SummaryCTimeHybridAll; SummaryR.SummaryCTimeHybrid];
    % end

  end
end

% get parent folder, then save summary results of all patients to .csv files in the parent folder
ParentFolder = D(1).folder;

% save all patient IDs
OutputFileName = 'Summary_Results_All_Patient_IDs.csv';
writetable(cell2table(PtIDAll), fullfile(ParentFolder, OutputFileName), 'writevariablenames', false, 'quotestrings', true);

% save column names of all patient summary results array
OutputFileName = 'Summary_Results_All_Patients_Array_Column_Names.csv';
SummaryResultsArrayColNames = {'numBscans', 'volScanW_mm', 'volScanH_mm', 'bScanW_pix', 'bScanH_pix', 'TestEyeNum', 'TestTime', 'FollowUpTime', 'FileModDateNumber', 'MeanThick', 'StdThick', 'MeanArea', 'StdArea', 'MeanVol', 'StdVol'};
writetable(cell2table(SummaryResultsArrayColNames), fullfile(ParentFolder, OutputFileName), 'writevariablenames', false, 'quotestrings', true);

% save column names of all patient summary results array
OutputFileName = 'Summary_Results_All_Patients_Array_Column_Names_Individual_Metric.csv';
SummaryResultsArrayColNames = {'numBscans', 'volScanW_mm', 'volScanH_mm', 'bScanW_pix', 'bScanH_pix', 'TestEyeNum', 'TestTime', 'FollowUpTime', 'FileModDateNumber', 'Mean', 'Std', 'Individual Measures'};
writetable(cell2table(SummaryResultsArrayColNames), fullfile(ParentFolder, OutputFileName), 'writevariablenames', false, 'quotestrings', true);

% save test time, follow-up time, mean and std of OS metrics measurements - Auto-segmentation
OutputFileName = 'Summary_Results_All_Patients_Average_Layer_Metrics_Auto.csv';
writematrix(SummaryLayerMetricsAutoAll, fullfile(ParentFolder, OutputFileName));

% save test time, follow-up time, mean and std of OS metrics measurements - Manual-segmentation
OutputFileName = 'Summary_Results_All_Patients_Average_Layer_Metrics_Manual.csv';
writematrix(SummaryLayerMetricsManuAll, fullfile(ParentFolder, OutputFileName));

% save test time, follow-up time, mean and std of OS metrics measurements - UNet
OutputFileName = 'Summary_Results_All_Patients_Average_Layer_Metrics_UNet.csv';
writematrix(SummaryLayerMetricsUNetAll, fullfile(ParentFolder, OutputFileName));

% save test time, follow-up time, mean and std of OS metrics measurements - Hybrid
OutputFileName = 'Summary_Results_All_Patients_Average_Layer_Metrics_Hybrid.csv';
writematrix(SummaryLayerMetricsHybridAll, fullfile(ParentFolder, OutputFileName));

% MOD on 07/06/2023 to save the summary for individual metrics mean, std,
% and individual measurements used to calculate mean & std

% save test time, follow-up time, mean, std, and individual OS metric measurements - Auto-segmentation
% OutputFileName = 'Summary_Results_All_Patients_Average_Thick_Auto.csv';
% writematrix(SummaryThickAutoAll, fullfile(ParentFolder, OutputFileName));
% OutputFileName = 'Summary_Results_All_Patients_Average_Area_Auto.csv';
% writematrix(SummaryAreaAutoAll, fullfile(ParentFolder, OutputFileName));
% OutputFileName = 'Summary_Results_All_Patients_Average_Vol_Auto.csv';
% writematrix(SummaryVolAutoAll, fullfile(ParentFolder, OutputFileName));

% save test time, follow-up time, mean and std of OS metrics measurements - Manual-segmentation
OutputFileName = 'Summary_Results_All_Patients_Average_Thick_Manu.csv';
writematrix(SummaryThickManuAll, fullfile(ParentFolder, OutputFileName));
OutputFileName = 'Summary_Results_All_Patients_Average_Area_Manu.csv';
writematrix(SummaryAreaManuAll, fullfile(ParentFolder, OutputFileName));
OutputFileName = 'Summary_Results_All_Patients_Average_Vol_Manu.csv';
writematrix(SummaryVolManuAll, fullfile(ParentFolder, OutputFileName));

% save test time, follow-up time, mean and std of OS metrics measurements - UNet
OutputFileName = 'Summary_Results_All_Patients_Average_Thick_UNet.csv';
writematrix(SummaryThickUNetAll, fullfile(ParentFolder, OutputFileName));
OutputFileName = 'Summary_Results_All_Patients_Average_Area_UNet.csv';
writematrix(SummaryAreaUNetAll, fullfile(ParentFolder, OutputFileName));
OutputFileName = 'Summary_Results_All_Patients_Average_Vol_UNet.csv';
writematrix(SummaryVolUNetAll, fullfile(ParentFolder, OutputFileName));

% save test time, follow-up time, mean and std of OS metrics measurements - Hybrid
OutputFileName = 'Summary_Results_All_Patients_Average_Thick_Hybrid.csv';
writematrix(SummaryThickHybridAll, fullfile(ParentFolder, OutputFileName));
OutputFileName = 'Summary_Results_All_Patients_Average_Area_Hybrid.csv';
writematrix(SummaryAreaHybridAll, fullfile(ParentFolder, OutputFileName));
OutputFileName = 'Summary_Results_All_Patients_Average_Vol_Hybrid.csv';
writematrix(SummaryVolHybridAll, fullfile(ParentFolder, OutputFileName));
% end of MOD

% ADDED on 07/24/2023
OutputFileName = 'Summary_Results_All_Patients_Classificatoin_Time_UNet.csv';
writematrix(SummaryCTimeUNetAll, fullfile(ParentFolder, OutputFileName));
OutputFileName = 'Summary_Results_All_Patients_Classificatoin_Time_Hybrid.csv';
writematrix(SummaryCTimeHybridAll, fullfile(ParentFolder, OutputFileName));
% end

warning('on');    % turn warning back on

% end of script_get_layer_metrics_summary_all_segmentations.m
