% script_check_OCTBScan_segmentation.m
%   to show the results of segmentation on OCTScan by different methods, 
%   including auto-, manual-, and DML model segmentations
%
% Yi-Zhong Wang, 06/11/2020; MOD on 10/31/2023 to check a single B-scan
% Retina Foundation of the Southwest
% Copyright. All Rights Reserved.
%

clearvars; clc; close all;

warning('off','all');       % turn all warnings off

% start file number
startFileNum = 1;

% B-scan to be checked
bScanNum = inf;             % inf -> check mid-line B-scan

% display B-scan image or not
showBscanImage = 1;

% segmentation lines to be checked
LineName = {'PR1', 'RPE'};

%**************************************************************************
% get OCTScan.mat files
%**************************************************************************
[dataDir, OCTScanList] = GetFileDirList('OCTScan*.mat', 'to check');
if dataDir==0
  return;
end
disp(['dataDir = ', dataDir]);
fprintf('\n');

%**************************************************************************
disp('loop through all OCTScan*.mat files to display B-Scan segmentation ...')
%**************************************************************************
for j=startFileNum:length(OCTScanList)
  % load OCTScan.mat to check segmentation
  FileName = OCTScanList(j).name;
  disp(['  Check segmentations of "', FileName, '"']);
  load(fullfile(dataDir, FileName), 'OCTScan');

  % display the first image (infrared)
  Im_1st = OCTScan.Images(1).ImageArray;
  figure; imshow(Im_1st);
  pause(0.1);
  % close all;

  % display the specified B-scan
  numBscans = OCTScan.numImages - 1;    % first image is infrared
  if bScanNum > numBscans
    i = round(numBscans/2);             % actual B-scan # to check
  end
  disp(['    Check segmentation of B-Scan # ', num2str(i)]);
  BScanImage = OCTScan.Images(i+1);       % B-scan image structure
  Im         = BScanImage.ImageArray;     % B-scan image array
  Segments   = HEOCTSortOCTScanSegmentation(BScanImage, FileName);
  
  % check auto segmentation if it exists
  indxEZ = find(contains(LineName, 'PR1'));
    
  % check auto segmentation if it exists
  if ~isempty(Segments) && ~isempty(Segments.Auto)
    % check all auto segmentation
    for k=1:length(Segments.Auto)
      GraderName = Segments.Auto(k).Grader;
      Segment = Segments.Auto(k).Segment;
      [PixCL, PixCA] = HEOCTBScanImagePixClass(Im, Segment, LineName);
      LayerM = GetLayerMeasure(PixCL);
      h = OCTScanCheckSegment(Im, Segment, LineName);
      GraderName(GraderName=='_') = '-';    % replace '_' with '-' for plot title display
      title(['Auto Segmentation by ', GraderName '; EZ width = ', num2str(LayerM.Width(indxEZ)), ' pixels']);
      set(h, 'Name', FileName);
      % set(h, 'NumberTitle', 'off');
    end
  end

  % check manual segmentation if it exists
  if ~isempty(Segments) && ~isempty(Segments.Manual)
    % check all manual segmentation
    for k=1:length(Segments.Manual)
      GraderName = Segments.Manual(k).Grader;
      Segment = Segments.Manual(k).Segment;
      [PixCL, PixCA] = HEOCTBScanImagePixClass(Im, Segment, LineName);
      LayerM = GetLayerMeasure(PixCL);
      h = OCTScanCheckSegment(Im, Segment, LineName);
      GraderName(GraderName=='_') = '-';    % replace '_' with '-' for plot title display
      title(['Manual Segmentation by ', GraderName '; EZ width = ', num2str(LayerM.Width(indxEZ)), ' pixels']);
      set(h, 'Name', FileName);
      % set(h, 'NumberTitle', 'off');
    end
  end

  % check DML model segmentation if it exisits
  if ~isempty(Segments) && ~isempty(Segments.ML)

    % check all DML model segmentation
    for k=1:length(Segments.ML)
      ModelName = Segments.ML(k).Model;
      ModelName(ModelName=='_') = '-';    % replace '_' with '-' for correct title text display
      
      if contains(ModelName, 'MatConvNet', 'IgnoreCase', true) || contains(ModelName, 'SW', 'IgnoreCase', true)
        Im_c    = Segments.ML(k).Im_c;
        ProbMap = Segments.ML(k).ProbMap;
        Param   = Segments.ML(k).Param;
        
        % indxEZ may be different
        indxEZ  = find(contains(Param.PatchParam.LineName, 'PR1'));

        % get results of segmentation
        if isfield(Segments.ML(k), 'SegmentData_SW') && ~isempty(Segments.ML(k).SegmentData_SW)
          LCASAResult = Segments.ML(k).SegmentData_SW;
        else
          LCASAResult = SWPostProcessPixelClasses(Im_c, ProbMap, Param);
        end

        LayerM = GetLayerMeasure(LCASAResult.Im_c_p);     % here is for plotting convenience, since Im_c_p is pixelized, not real type 
        h = OCTScanCheckSegment(Im, LCASAResult.Segment, LineName);
        title(['DML Segmentation by ', ModelName '; EZ width = ', num2str(LayerM.Width(indxEZ)), ' pixels']);
        set(h, 'Name', FileName);
        % set(h, 'NumberTitle', 'off');

      elseif contains(ModelName, 'UNET', 'IgnoreCase', true)

        % post-processing the classifcation results of the model
        Im_c    = Segments.ML(k).Im_c;
        ProbMap = Segments.ML(k).ProbMap;
        Param   = Segments.ML(k).Param;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % modify other parameters for post-processing
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Param.PostP.minPixels = 10;     % minimum # of pixels in a legit local area

        % indxEZ may be different
        indxEZ  = find(contains(Param.PatchParam.LineName, 'PR1'));
                      
        % UNet segmentation
        if isfield(Segments.ML(k), 'SegmentData') && ~isempty(Segments.ML(k).SegmentData)
          SegmentData = Segments.ML(k).SegmentData;
        elseif isfield(Segments.ML(k), 'SegmentData_UNet') && ~isempty(Segments.ML(k).SegmentData_UNet)
          SegmentData = Segments.ML(k).SegmentData_UNet;
        else
          SegmentData = UNetGetSegmentDataFromClass(Im_c, Param);
        end
        LayerM = GetLayerMeasure(SegmentData.PixCL);
        h = OCTScanCheckSegment(Im, SegmentData.Segment, LineName);
        title(['DML Segmentation by ', ModelName '; EZ width = ', num2str(LayerM.Width(indxEZ)), ' pixels']);
        set(h, 'Name', FileName);
        % set(h, 'NumberTitle', 'off');
        
        % the hybrid model segmentation, if exists
        if isfield(Segments.ML(k), 'SegmentData_Hybrid') && ~isempty(Segments.ML(k).SegmentData_Hybrid)
          SegmentData = Segments.ML(k).SegmentData_Hybrid;
          LayerM = GetLayerMeasure(SegmentData.PixCL);
          h = OCTScanCheckSegment(Im, SegmentData.Segment, LineName);
          title(['DML Hybrid Segmentation by ', ModelName '; EZ width = ', num2str(LayerM.Width(indxEZ)), ' pixels']);
          set(h, 'Name', FileName);
          % set(h, 'NumberTitle', 'off');
        end
      end
    end
         
  end

  if showBscanImage~=0
    % display B-scan image
    figure;
    imshow(Im);
  end
  
  pause;
  close all;

end

warning('on', 'all');     % turn all warnings back on
  
% end of script_check_OCTBScan_segmentation.m
