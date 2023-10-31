% script_classify_OCTScans_single_folder_paired_models.m
%   script to automatically classify all pixels of B-scan images of a list
%   of OCTScan.mat files in a single directory using a list of paired
%   UNet and SW models: UNet for initial segmentation and the SW model for
%   hybrid model classification.
%
% Yi-Zhong Wang, 03/28/2022
% Retina Foundation of the Southwest
% Copyright. All Rights Reserved.
%

clc; clearvars; close all;

% make sure the current directory is accessible so that uigetdir works
cd(fileparts(mfilename('fullpath')));


%**************************************************************************
% specify a list of paired UNET & SW models
%**************************************************************************
% model pairs for three DLM measurements
ModelPairs = {'RP340_T080_UNET256x32_O28_DLM01', 'RP340_T080_SW33x33_O32_DLM01'; ...
              'RP340_T080_UNET256x32_O28_DLM02', 'RP340_T080_SW33x33_O32_DLM02'; ...
              'RP340_T080_UNET256x32_O28_DLM03', 'RP340_T080_SW33x33_O32_DLM03'};


%**************************************************************************
% modify parameters for OCT B-scan image classification
%**************************************************************************
Param             = DefaultParamClassification;
Param.resultLog   = 0;    % display segmentation results or not
Param.tTop    	  = inf;  % # pixels above ILM; inf->all pixels above
Param.tBot    	  = inf;  % # pixels below BM; inf->all pixels below
Param.testWidth   = inf;  % central area (mm) to test; inf->full width
Param.addHybrid   = 1;    % ~=0->add hybrid model segmentation
Param.saveProbMap = 0;    % = 0 -> not save ProbMap to reduce file size


%**************************************************************************
% get trained U-Net models
%**************************************************************************
[netDir, ModelListAll] = GetFileDirList('*.mat', 'of trained DML models');
if netDir==0
  return;
end
disp(['netDir = ', netDir]);
fprintf('\n');

% only check U-Net models in ModelList
Indx = zeros(size(ModelListAll));
for i=1:length(Indx)
  if contains(ModelListAll(i).name, 'UNET', 'IgnoreCase', true)
    Indx(i) = 1;
  end
end
ModelList = ModelListAll(Indx==1);


%**************************************************************************
% get OCTScan.mat files for B-scan image classification
%**************************************************************************
[dataDir, OCTScanList] = GetFileDirList('OCTScan*.mat', 'for classification');
if dataDir==0
  return;
end
disp(['dataDir = ', dataDir]);
fprintf('\n');


%**********************************************************************
disp('loop through all OCTScan files for classification by DL models ...');
%**********************************************************************
fprintf('\n');
for i=1:length(OCTScanList)
  % load OCTScan.mat for classification
  FileName = OCTScanList(i).name;
  disp(['load OCTScan file ', FileName, ' ...']);
  load([dataDir, filesep, OCTScanList(i).name]);
  disp(['Classifying "', FileName, '"']);

  changeMade = 0;     % flag to indicate if change made to OCTScan or not

  %**************************************************************************
  disp('  loop through all U-Net models to perform classification ...');
  %**************************************************************************
  for k=1:length(ModelList)
    MLmodel = ModelList(k).name;
    MLmodelName = MLmodel(1:end-4);       % get rid of .mat extension
  
    % determine if MLmodel is in the specified model pairs
    % DLModelNames(:,1) -> access all 1st-column U-Net models as an array
    % DLModelNames{:,1} -> access all 1st-column models as indivisual cell
    indxM = find(ismember(ModelPairs(:,1), MLmodelName), 1);
    useModel = ~isempty(indxM);
    if useModel==0
      disp(['    ', MLmodelName, ' is not used ...']);
      fprintf('\n');
    else
      disp(['    Use DL model ', MLmodelName, ' for semantic classification ...']);
    end
    
    if useModel
      %**********************************************************************
      % load UNet then prepare it for classification
      %**********************************************************************
      disp(['    Load UNet model "', MLmodel, '" for semantic segmentation ...']);
      load([ModelList(k).folder, filesep, MLmodel], 'Net');
      
      % get net; note that loaded Net has two fields: .net & .meta
      net = Net.net;
  
      % update PatchParam and ModelParam
      Param.PatchParam = Net.meta.PatchParam;
      Param.ModelParam = Net.meta.ModelParam;
  
      % modify parameters for patch generation for U-Net classification
      Param.PatchParam.patchType = 3;   % 3 -> classification patches
      Param.PatchParam.overlap   = 0;   % no patch overlap for classification
      
      %**********************************************************************
      % get the SW model that is paired with UNet
      %**********************************************************************
      SWModelName = ModelPairs{indxM, 2};
      if any(size(dir([netDir, filesep, SWModelName, '.mat']),1))
        % SWModelName exists, load it
        disp(['    Load paired SW model "', SWModelName, '" for hybrid model refinement ...']);
        SWNet = load([netDir, filesep, SWModelName], 'Net');
        Param.SWNet = SWNet.Net;        % Net has two fields: .net & .meta
        Param.SWNet.SWModelName = SWModelName;        % add SW model name
      else
        % SWModelName doesn't exist, using default SW model
        disp('    default SW model will be used for hybrid model refinement ...');
        SWNet = load('SWNetDefault', 'Net');          % make sure assigned to SWNet, so it will not replace UNet Net
        Param.SWNet = SWNet.Net;
        Param.SWNet.SWModelName = 'SWNetDefault';     % add SW model name
      end
      fprintf('\n');
  

      %******************************************************************
      disp('      loop through all B-scans in OCTScan for classification ...');
      %******************************************************************
      numBscanImage = OCTScan.numImages - 1;
      changeMadeModel = 0;    % flag to indicate change made to B-scans for the current model - consider multiple-models
      for j=1:numBscanImage
        BScanImage = OCTScan.Images(j+1);
        Im_orig = BScanImage.ImageArray;

        %******************************************************************
        % MOD on 01/13/2022 by YZW to handle no signal area in Im, i.e., 
        % those areas filled with 255 (no signal), e.g., RUSH2A data.
        % After pre-processing, need to assign Im back to ImageArray for DL
        % model classification
        Im = PrePHandleNoScanSignalAreas(Im_orig);    % pre-prcessing Im
        BScanImage.ImageArray = Im;                   % assign back
        % end of MOD
        %******************************************************************

        % update Param with B-scan parameters
        Param.BScanParam.ScaleXY      = BScanImage.ScaleXY;
        if isfield(BScanImage, 'StarCoordXY')
          Param.BScanParam.StartCoordXY = BScanImage.StartCoordXY;
        end
        if isfield(BScanImage, 'EndCoordXY')
          Param.BScanParam.EndCoordXY   = BScanImage.EndCoordXY;
        end

        % check if BScanImage already contains the classification results
        % by the same model
        if isfield(BScanImage, 'Segment') && ~isempty(BScanImage.Segment)   % considering no auto segment
          Fields = fieldnames(BScanImage.Segment);
          Indx_ML = find(contains(Fields, ['ML_', MLmodelName]));
        else
          Indx_ML = [];
        end

        % classification by the model if not been done before
        if isempty(Indx_ML)
          % display classification for each B-scan only for a signle DL model pair
          if ~isempty(ModelPairs) && length(ModelPairs)==1
            disp(['      Classify B-scan image #', num2str(j), ' ...']);
          end

          ML_seg = ClassifyBScanUNetHE(net, BScanImage, Param, FileName);

          % using dynamic field name to add a new field of ML
          OCTScan.Images(j+1).Segment.(['ML_', MLmodelName]) = ML_seg;

          if Param.resultLog
            % display the pixel classification
            disp(['      time used for semantic segmentation = ', num2str(Sec)]);
            disp(['      time used for refinement classification = ', num2str(Sec2)]);
            figure;
            imshow(Im_c, [1, max(max(Im_c))]);
            title(FileName);
            colormap(hot);
          end

          % set flag
          changeMade = 1;         % change-made to OCTScan
          changeMadeModel = 1;		% change-made for the current model; note that changeMadeModel reset for each model
        end
      end
    end

    if changeMadeModel==0
      disp(['      Classification by the model "', MLmodelName, '" already done']);
    else
      disp(['      Classification by the model "', MLmodelName, '" completed']);
    end
    fprintf('\n');
    
  end

  % save updated OCTScan back to replace the previous OCTScan*.mat file
  if changeMade==1
    disp(['Change made to ', FileName]);
    disp(['Save modified ', FileName, ' ...']);
    save([dataDir, filesep, FileName], 'OCTScan', '-v7.3');
    fprintf('\n');
  else
    disp(['No change made to ', FileName]);
    fprintf('\n');
  end
  fprintf('\n');
end

% end of script_classify_OCTScans_single_folder_paired_models.m
