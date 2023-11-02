% script_OCTScan_Display_Layer_Metrics.m
%   Display layer metrics in a volume scan, including 1-D layer thickness 
%   of mid-line Bscan, 2-D layer area and 3-D layer volume for a selected
%   OCTScan*.mat file
%
% Yi-Zhong Wang, 09/12/2023
% Retina Foundation of the Southwest
% Copyright. All Rights Reserved.
% 

clearvars; clc; close all;

% specify the type of segmenation to plot
% SegType = 'Auto';       % HE auto segmentation
% SegType = 'Manu';       % manual segmentation
% SegType = 'DMLM';       % UNet segmentation
SegType = 'HDMLM';        % hybrid model segmentation

% specify two segmentation lines for thickness analysis
LN1 = 'PR1';
LN2 = 'RPE';
LineName = {LN1, LN2};

% other parameters
scanWidth      = inf;   % scan width in mm to be processed
scanPlotXRange = 8;     % X-axis plot range in mm for 3-D layer map
scanPlotYRange = 8;     % Y-axis plot range in mm for 3-D layer map
scanPlotZLimit = 50;    % Z-axis limit in micrometer for 3-D layer map

% post-processing to remove isolated EZ segments or not for DLM segment
elimIsoEZSegment = 1;

% which segment in SegType to display
segNum = 1;                 % which segment in SegType

% plot parameters
PlotParam = DefaultOCTVolScanPlotParam;
PlotParam.limitX      = scanPlotXRange/2;
PlotParam.limitY      = scanPlotYRange/2;
PlotParam.limitZ      = scanPlotZLimit;  

%**************************************************************************
disp('select an OCTScan file to process ...');
%**************************************************************************
[FileName, PathName] = uigetfile('*.mat', 'Select an OCTScan file to open');
if FileName==0
  disp('no OCTScan file selected ...');
  return;
end
load([PathName, filesep, FileName]);    % load OCTScan*.mat

% get infrared image
IrIm = OCTScan.Images(1).ImageArray;
if ~isempty(IrIm)
  [m, n] = size(OCTScan.Images(1).ImageArray(:,:,1));
else
  % for UHR OCT
  m = 1800;
  n = 1800;
end

% get mid-line B-scan image
numBscans = length(OCTScan.Images) - 1;
midLineNum = floor(numBscans/2) + 1;
BScanImage = OCTScan.Images(midLineNum+1);
% BScanImage = OCTScan.Images(midLineNum);        % one line below
Im = BScanImage.ImageArray(:,:,1);      % get one-channel for gray scale

% get segments
Segments = HEOCTSortOCTScanSegmentation(BScanImage);
switch SegType
  case 'Auto'
    if ~isempty(Segments.Auto)
      kA = min(segNum, length(Segments.Auto));
      AutoSegName = Segments.AutoSegNames{kA};
      Segment = Segments.Auto(kA).Segment;
      disp(['OS metrics measurements from ', AutoSegName, ' auto segmentation']);
    else
      disp('no HE auto segmentation. exit ...');
      return;
    end
  case 'Manu'
    if ~isempty(Segments.Manual)
      kM = min(segNum, length(Segments.Manual));
      ManuSegName = Segments.ManualGraders{kM};
      Segment = Segments.Manual(kM).Segment;
      disp(['OS metrics measurements from ', ManuSegName, ' manual segmentation']);
    else
      disp('no manual segmentation. exit ...');
      return;
    end
  case 'DMLM'
    if ~isempty(Segments.ML)
      kD = min(segNum, length(Segments.ML));
      DLMSegName = Segments.MLModelNames{kD};
      Segment = Segments.ML(kD).SegmentData.Segment;
      disp(['OS metrics measurements from ', DLMSegName, ' UNet segmentation']);
      if elimIsoEZSegment~=0
        % remove isolated EZ segments
        OCTScan = UNetPPDeleteIsolatedLayerSegment(OCTScan, LN1, LN2, 'DMLM', kD);
      end
    else
      disp('no UNet segmentation. exit ...');
      return;
    end
  case 'HDMLM'
    if ~isempty(Segments.ML) && isfield(Segments.ML, 'SegmentData_Hybrid')
      kD = min(segNum, length(Segments.ML));
      DLMSegName = Segments.MLModelNames{kD};
      Segment = Segments.ML(kD).SegmentData_Hybrid.Segment;
      disp(['OS metrics measurements from ', DLMSegName, ' hybrid model segmentation']);
      if elimIsoEZSegment~=0
        % remove isolated EZ segments
        OCTScan = UNetPPDeleteIsolatedLayerSegment(OCTScan, LN1, LN2, 'HDMLM', kD);
      end
    else
      disp('no hybrid model segmentation. exit ...');
      return;
    end
  otherwise
    disp('no segmentation type defined. exit ...');
    return;
end

disp(['get thickness results of ', SegType, ' segmentation ...']);
Thickness = OCTScanGetLayerThickness(OCTScan, LN1, LN2, SegType);

% extract [Xs_mm, Ys_mm, Zs_mu]
Xs_mm = Thickness.Xs_mm;
Ys_mm = Thickness.Ys_mm;
Zs_mu = Thickness.Zs_mu;

% extract central area of the layer [Xc, Yc, Zc] specified by scanWidth
[Xc, Yc, Zc] = GetCentralPart(Xs_mm, Ys_mm, Zs_mu, scanWidth);

% make Yc and Xc the same dimension as the infrared image size [m, n] by
% appending rows to Xc and columns to Yc to generate figures for DLM_MS04
Xct = zeros(m, n);
Yct = zeros(m, n);
Xc1d = Xc(1, :);
Yc1d = Xc1d';
for i=1:m
  Xct(i, :) = Xc1d;
end
for i=1:n
  Yct(:, i) = Yc1d;
end

% expand Zc by padding zeros -> assume EZ zone is limited
Zct = AdjustSize(Zc, m, n, 0);
maxZct = max(max(Zct));
Zct_3c = cat(3, 255*Zct/maxZct, 255*Zct/maxZct, 255*Zct/maxZct);


%**************************************************************************
% generate OS thickness, EZ area and OS volume plots
%**************************************************************************
% OS thickness on a B-scan image
h = OCTScanCheckAreaSegmentation(Im, Segment, LineName);

% show EZ area with intensity normalized to maximum OS thickness
figure; imshow(Zct, [0, max(max(Zct))]);
% figure; imshow(Zc, [0, scanPlotZLimit]);

% infuse infrared image with EZ area
if ~isempty(IrIm)
  EZ_EnFace = imfuse(IrIm, Zct_3c, 'blend', 'Scaling', 'joint');
  figure; imshow(EZ_EnFace);
end

% display 2-D area and 3-D surface - thickness normalized to scanPlotZLimit
Hc = PlotThicknessData(Xct, Yct, Zct, PlotParam);

% end of script_OCTScan_Display_Layer_Metrics.m
