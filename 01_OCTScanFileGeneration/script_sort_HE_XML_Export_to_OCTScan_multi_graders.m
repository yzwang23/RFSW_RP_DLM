% script_sort_HE_XML_Export_to_OCTScan_multi_graders.m
% 	a script to process HE OCT scan XML export files stored in a specified
%   folder to generate OCTScan data structure, then save sorted structure
%   OCTScan in an OCTScan*.mat file in the specified destination folder.
%   If a destination folder is not specified,  the sorted OCTScan will be
%   saved in a default folder under the parent folder of the XML files. 
%   If the destination folder contains any OCTScan*.mat files, this script
%   will check if any of the existing files is from the same OCT XML export
%   by checking the patient ID & the date & time stamps. If no same export
%   exsits, a new file will be created. If yes, but from a different
%   grader, new grader's result will be added to OCTScan, then the updated 
%   OCTScan will be saved under the same .mat file name, so that the new 
%   OCTScan*.mat file may contain various types of segmentation, including 
%   auto segmentation, manual segmentation from differen graders, as well 
%   as the segmentation by deep machine learing models added later.
%
% Yi-Zhong Wang, 12/11/2018; 07/14/2021; new comments added on 10/26/2023
% Retina Foundation of the Southwest
% Copyright. All Rights Reserved.

% Note - the following are the naming conventions of the folder that stores
%        XML exports:
%   XML_Diagnosis_[SegmentType_Grader]_[ID_xxxxx]_[others]......
%     Diagnosis:    adRP / arRP / xlRP / isoRP / AMD ......
%     SegmentType:  Auto / Manu / ML ...... (fixed tags)
%     Grader:       HE / GraderName / ModelName ......
%     ID:           ID tag (fixed)
%     xxxxx:        ID number
%     [othres]:     additional fields which are not processed by the script
%   Note that the folder name could contain multiple segment types, e.g.
%     XML_xlRP_Auto_HE_Manu_GraderName_......
%	since the script will sort XML files into auto and manual
%
% Note - file naming for OCTScan*.mat file, to be generated automatically
%     OCTScan_ScanType_Diagnosis_ID_Eye_Date_Time
%   where
%     OCTScan:     'OCTScan' as a fixed field
%     ScanType:    Line ('LnH'/'LnV') or volume ('VoH'/'VoV')
%     Diagnosis:   'Normal', 'AMD', 'DR', 'RP', 'adRP', or others
%     ID:          RFSW/Clinical Trial ID shown on HE OCT
%     Eye:         'OD' or 'OS'
%     Date:        in format 'yyyymmdd'
%     Time:        in format 'hhmmssfs' (fs: fraction of second)
% 
% Note - combining automatic segmentation and manual segmentation together
%   Automatic segmentation is stored as default for each B-scan image under
%       OCTScan.Images(i).Segment
%   Manually corrected segmentation is stored in a 'Grader_XXX' subfield
%       OCTScan.Images(i).Segment.Grader_XXX
%   Note that the "XML_..." folder may contain XML exports (.xml files) of
%   both automatic and manually-corrected segmentation. This script will
%   first identify the type of .xml files (auto or manu) based on the field
%   of Segment.Manual
%
% MOD
%   01/24/2019 - make sure B-scan images have black background
%   02/21/2019 - add handling of sub-folders of XML files
%                change folder naming order: move Diagnosis to after Grader
%                change OCTScan file naming order: move ScanType to 2nd
%                check for what type of scan based on XML output and get
%                rid of ScanType in the folder name
%                make a default output folder to save OCTScan*.mat files
%                when no output directory is specified
%   11/09/2020 - add .FileName field to variable OCTScan for saving
%   07/14/2021 - add the capability to handle auto and manual segmentation
%                XML exports in the same folder

clearvars; clc; close all;

% make sure the current directory is accessible so that uigetdir works
cd(fileparts(mfilename('fullpath')));

% flag to save sorted OCTScan as a combination file or an individual file
combFile = 1;

%**************************************************************************
disp('specify a directory containing OCT XML exports ...');
%**************************************************************************
DirFullName = uigetdir;
if DirFullName==0
  disp('OCT XML exports folder not specified ...');
  return;
end
DirPathName = [DirFullName, filesep];

% find all .xml files in the DirPathName folder
FileNamesAll = dir([DirPathName, '*.xml']);

% get rid of files with leading '.'
count = 0;
for i=1:length(FileNamesAll)
  if ~strcmp(FileNamesAll(i).name(1), '.')
    count = count + 1;
    FileNames(count) = FileNamesAll(i);
  end
end

numXMLFiles = length(FileNames);

%**************************************************************************
disp('specify a directory to save sorted OCTScan*.mat files ...');
% then get output file path. If no output folder is specified, each sorted
% file will be saved in a default folder created in the parent folder of
% the XML folder. the name of the default folder is created based on XML
%	folder name as well as the name of the parenet folder of the XML folder
%**************************************************************************
OutputDirName = uigetdir;
if OutputDirName==0
  Indx = find(DirFullName==filesep);
  XMLFolderName = DirFullName(Indx(end)+1:length(DirFullName));
  XMLParentFolderName = DirFullName(Indx(end-1)+1:Indx(end)-1);
  XMLParentFolderPath = [DirFullName(1:Indx(end)-1), filesep];
  Indx2 = find(XMLFolderName=='_');
  OutputFolderName = ['OCTScans', XMLFolderName(Indx2(1):end)];
  RFSWPsychMakeDirectory([XMLParentFolderPath, OutputFolderName]);	% create the default folder
  OutputPathName = [XMLParentFolderPath, OutputFolderName, filesep]; % output path
else
  OutputPathName = [OutputDirName, filesep];
end

% change current folder to output folder for file saving
OldDir = cd;            % get current dir
cd(OutputPathName);

% analysize XML export files
% AnalysisResults = AnalyzeXMLExportFiles(DirPathName);
AnalysisResults = AnalyzeXMLExportFilesDetails(DirPathName);
Diagnosis       = AnalysisResults.Diagnosis;
ManuFiles       = AnalysisResults.ManuFiles;    % flags for manual files
AutoFiles       = abs(1 - ManuFiles);           % flags for auto files

fprintf('\n');
numAutoFiles = length(find(AutoFiles==1));
%**************************************************************************
disp(['processing automatic segmentation (n = ', num2str(numAutoFiles), ') ...']);
%**************************************************************************
for j=1:length(AutoFiles)
  if AutoFiles(j)==1
    FileName = FileNames(j).name;
  
    % convert XML export to OCTScan
    [OCTScan, OCTScanFileName] = XMLExportToOCTScan(DirPathName, FileName, Diagnosis);
  
    % check if OCTScanFileName exists in the destination folder
    if ~exist([OutputPathName, OCTScanFileName, '.mat'], 'file')	% not exist
      disp('    OCTScan file not exist ...')
      disp('    save OCTScan to a .mat file in the output folder ...');
      save(OCTScanFileName, 'OCTScan', '-v7.3');

    else        % if OCTScanFileName already existed
      disp('    OCTScan file exists ...')

      % load previously saved OCTScan file into OCTScan_pre
      OCTScan_temp = load(OCTScanFileName, 'OCTScan');
      OCTScan_pre  = OCTScan_temp.OCTScan;

      % check to see if OCTScan_pre already has auto-segment results
      % if not exist, copy auto-segment results to OCTScan_pre
      % if exists, replace the previous one
      for i=1:length(OCTScan.Images)
        if ~isempty(OCTScan.Images(i).Segment)
          OCTScan_pre.Images(i).Segment.NumLines   = OCTScan.Images(i).Segment.NumLines;
          OCTScan_pre.Images(i).Segment.LineName   = OCTScan.Images(i).Segment.LineName;
          OCTScan_pre.Images(i).Segment.Enabled    = OCTScan.Images(i).Segment.Enabled;
          OCTScan_pre.Images(i).Segment.Manual     = OCTScan.Images(i).Segment.Manual;
          OCTScan_pre.Images(i).Segment.DataPoints = OCTScan.Images(i).Segment.DataPoints;
          OCTScan_pre.Images(i).Segment.DataType   = OCTScan.Images(i).Segment.DataType;
          OCTScan_pre.Images(i).Segment.Data       = OCTScan.Images(i).Segment.Data;                  
        end
      end
      
      disp('    save updated OCTScan to replace previous OCTScan ...');
      OCTScan = OCTScan_pre;      % OCTScan_pre is now updated
      save(OCTScanFileName, 'OCTScan', '-v7.3');

    end
  end
end

fprintf('\n');
numManuFiles = length(find(ManuFiles==1));
%**************************************************************************
disp(['processing manual segmentation (n = ', num2str(numManuFiles), ') ...']);
%**************************************************************************
SegType = 'Manu';
Grader  = AnalysisResults.ManuGrader;
for j=1:length(ManuFiles)
  if ManuFiles(j)==1
    FileName = FileNames(j).name;
  
    % convert XML export to OCTScan
    [OCTScan, OCTScanFileName] = XMLExportToOCTScan(DirPathName, FileName, Diagnosis);
  
    % check if OCTScanFileName exists in the destination folder
        if ~exist([OutputPathName, OCTScanFileName, '.mat'], 'file')	% not exist
          disp('    OCTScan file not exist ...')
          
          % get grader initial and grader_init subfield name
          Grader_Init = ['Grader_', Grader];

          % move Segment field to Segment.Grader_Init
          OCTScan_temp = OCTScan;         % store OCTScan to a temp var
          for i=1:length(OCTScan.Images)
            OCTScan.Images(i).Segment = [];
            OCTScan.Images(i).Segment.(Grader_Init) = OCTScan_temp.Images(i).Segment;
          end

          disp('    save OCTScan to a .mat file in the folder ...');
          save(OCTScanFileName, 'OCTScan', '-v7.3');
            
        else        % if OCTScanFileName already existed
          disp('    OCTScan file exists ...')

          % load previously saved OCTScan file into OCTScan_pre
          OCTScan_temp = load(OCTScanFileName, 'OCTScan');
          OCTScan_pre  = OCTScan_temp.OCTScan;

          % get grader initial and grader_init subfield name
          Grader_Init = ['Grader_', Grader];

          % replace previous Grader_Init, if any, with the current one,
          % no matter that previous Grader_Init exists or not
          % copy current Segment field to Segment.Grader_Init
          for i=1:length(OCTScan.Images)
            OCTScan_pre.Images(i).Segment.(Grader_Init) = OCTScan.Images(i).Segment;
          end
          
          disp('    save updated OCTScan to replace previous OCTScan ...');
          OCTScan = OCTScan_pre;      % OCTScan_pre is now updated
          save(OCTScanFileName, 'OCTScan', '-v7.3');

        end
  end
end

% restore to previous dir
cd(OldDir);

% end of script_sort_HE_XML_Export_to_OCTScan_multi_graders.m
