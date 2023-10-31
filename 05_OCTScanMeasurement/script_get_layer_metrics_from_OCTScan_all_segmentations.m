% script_get_layer_metrics_from_OCTScan_all_segmentations.m
%   script to obtain layer metrics measurements from a list of OCTScan*.mat
%   files that have been segmented by various approaches, including auto,
%   manual, and/or DL models. OCTScan*.m files can be in a single selected 
%   folder or in multiple subfolders under a specified directory (e.g., for
%   longitudinal follow-up, each patient's results are in one folder). This
%   script will generate a measurement 'Results_' folder in the same folder
%   that contains OCTScan*.mat files.
%
% Yi-Zhong Wang, 10/05/2022
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

% scan width
d = inf;       % =inf -> full scan

%**************************************************************************
% specify pairs of segmentation lines for area and volume estimates, as
% well as type of segmentation to be evaluated
%**************************************************************************
LinePair = {'PR1', 'RPE'};      % photoreceptor outer segment [top, bottom]

%**************************************************************************
% trained DL models to be evaluated - full model name
%**************************************************************************
% DLModelNames = {'RP240_80_20_UNET'}';
DLModelNames = [];          % evaluate all DL models if empty


%**************************************************************************
disp('specify a directory containing OCTScan*.mat files or folders ...');
%**************************************************************************
DirFullName = uigetdir;
if DirFullName==0
  disp('OCTScan folder not specified ...');
  return;
end
DirPathName = [DirFullName, filesep];
disp(['Directory: ' DirFullName, ' ...']);

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
% perform classification, need to determine a single folder or not
%**************************************************************************p
if D(1).isdir==0        % not a sub-directory, a list of OCTScan*.mat files
  % only keep files starting with 'OCTScan_'
  Indx = zeros(length(D), 1);
  for i=1:length(D)
    if strcmpi(D(i).name(1:8), 'OCTScan_')
      Indx(i) = 1;
    end
  end
  D = D(Indx==1);

  disp('OCTScan*.mat files are in the selected folder. Get layer measurements from these OCTScan files ...');
  fprintf('\n');
  OCTScanList = D;

  % call MeasureOCTScans.m to get layer metrics measurements
  Results = MeasureOCTScans(OCTScanList, DLModelNames, LinePair, d, 0);
else
  disp('OCTScan*.mat files are in the subfolders under the selected folder ...');
  disp('  loop through all subfolders to get layer measurements ...');
  fprintf('\n');
  numSubfolder = length(D);
  for l=1:numSubfolder
    if D(l).isdir~=0
      disp(['***Sub-folder name: ', D(l).name]);
      OCTScanSubFullName = [D(l).folder, filesep, D(l).name];
      OCTScanList = dir([OCTScanSubFullName, filesep, 'OCTScan*.mat']);
  
      % call MeasureOCTScans.m to get layer metrics measurements
      Results = MeasureOCTScans(OCTScanList, DLModelNames, LinePair, d, 1);
    end
  end
end

warning('on');    % turn warning back on

% end of script_get_layer_metrics_from_OCTScan_all_segmentations.m
