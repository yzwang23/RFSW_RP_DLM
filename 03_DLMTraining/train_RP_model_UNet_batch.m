% train_RP_model_UNet_batch.m
%   script to train UNet model with image patch datasets generated from 
%   OCTScan*.mat files using script_generate_datasets_UNet
%
% by YZW, 12/16/2020; additional comments on 10/30/2023
% Retina Foundation of the Southwest
% Copyright. All Rights Reserved.
%

clc; clearvars; close all;

% get running script name for record
RunningScriptName = mfilename;

% make sure the current directory is accessible so that uigetdir works
cd(fileparts(mfilename('fullpath')));

% get computer name for record
[~, ComputerName] = system('hostname');

% shuffle random seed generator for different weights initialization
rng shuffle;

% using testing data in a separate directory or not
sepTestDir = 0;	% =0 -> testing data in the same directory as training data
                % if sepTestDir==0 and no testing data in training data
                % directory, then no testing/validating during training

% epochs to run 
epochsToRun   = 5;   % number epochs to run - default = 45
epochInterval = 5;    % interval of epochs to be finalized; default = 5; (e.g., =5 -> every 5 iteration to get a finalized Net), including the first and last epochs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('define parameters for the U-Net model structure and training ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParamUNet                           = DefaultParamUNet;
ParamUNet.eyeDisease                = 'RP';
ParamUNet.epochsToRun               = epochsToRun;
ParamUNet.convFilterSize            = [5, 5];
ParamUNet.encoderDepth              = 4;
ParamUNet.initialEncoderNumChannels = 8;
ParamUNet.Verbose                   = 1;
ParamUNet.VerboseFrequency          = 100;
ParamUNet.validation                = 1;        % validation takes long time for large validataion dataset
ParamUNet.numValidImages            = inf;      % number of images in validation set used; =inf -> use all validation images
ParamUNet.initialLearningRate       = 0.01;     % using 0.01 for UNet; if using 0.001 as in the SW model training, the accuracy is slightly worse at beginning
ParamUNet.learnRateDropPeriod       = 10;       % how many epochs to change learting rate
ParamUNet.maxEpochs                 = 1;        % number of epochs to run at a time by TrainOneEpoch or TrainMultiEpochs; =1->epoch to be finalized; otherwise # epochs between finalized epochs
ParamUNet.miniBatchSize             = 128;      % number of patches in a batch
ParamUNet.Plots                     = 'none';   % 'training-progress' or 'none'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('get parent directory of datasets and list of datasets ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parentDataDir = uigetdir;
if parentDataDir==0
  disp('no directory selected ...');
  return;
end
dataSetList = dir(parentDataDir);    	% get the list of data folders

% get indices to dataset folders; folder name contains string '_Train'
count = 0;
IndxDataSetDir  = [];
for i=1:length(dataSetList)
  if contains(dataSetList(i).name, '_Train') && dataSetList(i).isdir
    count = count+1;
    IndxDataSetDir(count) = i;
  end
end
numDatasets = length(IndxDataSetDir);
disp(['number of datasets = ', num2str(numDatasets)]);
fprintf('\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop through data batch folders to train CNN model using each dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:numDatasets
  % get directories to the folders of training data batches and results
  dataDir = [parentDataDir, filesep, dataSetList(IndxDataSetDir(j)).name];

  % get the type of eye disease and training dataset from dataSetDir
  idx = find(dataDir==filesep);
  dataSetFolderName = dataDir(idx(end)+1:length(dataDir));
  Indx = find(dataSetFolderName=='_');
  eyeDisease = dataSetFolderName(1:Indx(1)-1);
  ParamUNet.eyeDisease = eyeDisease;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('load in metadata and OCTScan file names and update parameters for UNet CNN ...');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  load([dataDir, filesep, 'metadata.mat'], 'metadata');
  load([dataDir, filesep, 'All_OCTScan_Files.mat'], 'OCTScanFileNames');    % note that All_OCTScan_Files.mat contains 'OCTScanFileNames' and other variables
  PatchParam           = metadata.Param;
  ParamUNet.imageSize  = [metadata.Param.m, metadata.Param.n];
  ParamUNet.numClasses = length(metadata.ClassName);
  ParamUNet.classNames = metadata.ClassName;
  ParamUNet.labelIDs	 = metadata.ClassLabel;

  % get number of iterations in an epoch and set validateFreq
  trainLabelDir          = fullfile(dataDir, 'TrainingLabels');
  numTrainImages         = length(dir([trainLabelDir, filesep, '*.png']));
  validateFreq           = floor(numTrainImages/ParamUNet.miniBatchSize);
  ParamUNet.validateFreq = validateFreq;

  % from PatchParam.tvRatio to percent training image patches
  tvRatio = PatchParam.tvRatio;
  if tvRatio<inf
    PerCTrain = ['_T0', num2str(round(100 * tvRatio ./ (tvRatio + 1)))];
  else
    PerCTrain = '_T100';
  end

  % strings as part of subfolder/file names
  ImageSizeStr = [PerCTrain, '_UNET', num2str(ParamUNet.imageSize(1)), 'x', num2str(ParamUNet.imageSize(2)), '_'];
  
  % convert Param.overlap in fraction or pixels to overlap in pixels
  % MOD by YZW on 08/02/2021 for overlap>1
  overlap = PatchParam.overlap;
  n = ParamUNet.imageSize(2);
  [overlapPixel, OverlapStr] = ConvertOverlapToPixel(overlap, n);
  % end of MOD

  FltSizeStr      = ['f', num2str(ParamUNet.convFilterSize(1))];
  EncoderDepthStr = ['e', num2str(ParamUNet.encoderDepth)];
  if ParamUNet.initialEncoderNumChannels<10
    InitChanNumStr  = ['i0', num2str(ParamUNet.initialEncoderNumChannels)];
  else
    InitChanNumStr  = ['i', num2str(ParamUNet.initialEncoderNumChannels)];
  end

  % get sub-folder name for U-Net training results and make the subfolder
  % the naming convention for U-Net result folder:
  %   [disease_imageSize_filterSize_encoderDepth_initialChannels]
  netFolderName = [ParamUNet.eyeDisease, ImageSizeStr, OverlapStr, '_', FltSizeStr, EncoderDepthStr, InitChanNumStr];
  netDir = [dataDir, filesep, netFolderName];
  if ~isfolder(netDir)
    mkdir(netDir);
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp(['check if there is any checkpoint and trained net for ''', dataSetFolderName, ''' ...']);
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~, EpochsDone] = GetTrainedNet(netDir);

  % get starting epoch #
  if ~isempty(EpochsDone)
    startEpoch = EpochsDone.numEpochsDone+1;
  else
    startEpoch = 1;
  end

  if startEpoch>ParamUNet.epochsToRun
    disp('all epochs are done;  check next dataset');
    
  else        % still having epochs to run
    % assemble Param that includes ParamUNet and PatchParam
    Param.PatchParam       = PatchParam;
    Param.ModelParam       = ParamUNet;
    Param.OCTScanFileNames = OCTScanFileNames;    % added on 12/19/2021; store names of OCTScan files used to generate the training dataset
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('create data source for network training ...');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % directories to training OCT image patches and corresponding pixel labels
    trainImageDir = fullfile(dataDir, 'TrainingImages');
    trainLabelDir = fullfile(dataDir, 'TrainingLabels');

    % create an imageDatastore holding the training images
    imds = imageDatastore(trainImageDir);

    % create a pixelLabelDatastore holding the ground truth pixel labels for
    % the training images
    pxds = pixelLabelDatastore(trainLabelDir,ParamUNet.classNames,ParamUNet.labelIDs);

    % create data source for training a semantic segmentation network
    ds = pixelLabelImageDatastore(imds,pxds);

    % determine validation dataset
    ds_v = [];            % init validation dataset; default -> no validation
    if ParamUNet.validation~=0
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      disp('create data source for network validation ...');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % get test data directory
      if sepTestDir==0
        % test data directory is the same as training data directory
        testDataDir = dataDir;
      else  
        % testing batches in a separate directory, and its default location
        % is under the parent data directory and having 'Test' in the folder
        % name. Search dataSetList to find the 'Test' subfolder
        count = 0;
        IndxTestDataSetDir  = [];
        for k=1:length(dataSetList)
          if contains(dataSetList(k).name, '_Test')
            count = count+1;
            IndxTestDataSetDir(count) = k;
          end
        end
        numTestDatasets = length(IndxTestDataSetDir);
        disp(['number of test datasets = ', num2str(numTestDatasets)]);

        if ~isempty(IndxValidateDataSetDir)
          testDataDir = [parentDataDir, filesep, dataSetList(IndxTestDataSetDir(1)).name];
        else
          testDataDir = parentDataDir;
        end
      end

      % directories to validating image patches and corresponding pixel labels
      validateImageDir = fullfile(testDataDir, 'ValidateImages');
      validateLabelDir = fullfile(testDataDir, 'ValidateLabels');

      if exist(validateImageDir, 'dir') && exist(validateLabelDir, 'dir')
        if length(dir(validateLabelDir))>2
          % create an imageDatastore holding the validating images
          imds_v = imageDatastore(validateImageDir);

          % get images of first numValidImages
          imds_v.Files = imds_v.Files(1:min(ParamUNet.numValidImages, length(imds_v.Files)));

          % create a pixelLabelDatastore holding the ground truth pixel labels for
          % the validating images
          pxds_v = pixelLabelDatastore(validateLabelDir, ParamUNet.classNames, ParamUNet.labelIDs);
          Files_P = pxds_v.Files(1:min(ParamUNet.numValidImages, length(pxds_v.Files)));
          pxds_v = pixelLabelDatastore(Files_P, ParamUNet.classNames, ParamUNet.labelIDs);

          % create data source for validating a semantic segmentation network
          ds_v = pixelLabelImageDatastore(imds_v, pxds_v);
        end
      end
    end
	  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get a list of remaining epochs to run - added on 09/28/2022 by YZW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EpochRunList = GetEpochRunList(startEpoch, epochsToRun, epochInterval);
    numEpochRunGroups = EpochRunList.numEpochRunGroups;
    RemainEpochs      = EpochRunList.RemainEpochs;
    FinalizeFlags     = EpochRunList.FinalizeFlags;
    NumEpochsInGroup  = EpochRunList.NumEpochsInGroup;
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop through number of remaining epochs to train UNet
    % MOD on 09/28/2022 by YZW for group epochs training
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:numEpochRunGroups
      indxToRemainEpochsGroupStart = sum(NumEpochsInGroup(1:i)) - NumEpochsInGroup(i) + 1;
      indxToRemainEpochsGroupEnd   = sum(NumEpochsInGroup(1:i));        % accumulated to the last epoch in the group
      currentEpochNumStart  = RemainEpochs(indxToRemainEpochsGroupStart);
      currentEpochNumEnd    = RemainEpochs(indxToRemainEpochsGroupEnd);  % currentEpochNum is the last Epoch in the group

      % update epoch number for the current model training
      Param.ModelParam.currentEpochNum = currentEpochNumEnd;    % this is the number for the finalized epoch in the group

      % update .maxEpochs in ParamUNet and Param.ModelParam with NumEpochsInGroup(i)
      ParamUNet.maxEpochs        = NumEpochsInGroup(i);
      Param.ModelParam.maxEpochs = NumEpochsInGroup(i);

      % run validation for the groups with odd anchor epochs (1, 5, 15,...)
      % epochs in the same group either with or without validation
      % Also make sure to run validation for the last epoch group
      Param.ModelParam.validation = mod(currentEpochNumEnd, 2);
      if i==numEpochRunGroups
        Param.ModelParam.validation = 1;
      end
  
      % get starting and ending epoch number of the current group
      if currentEpochNumStart<10
        EpochStrStart = ['epoch0', num2str(currentEpochNumStart)];
      else
        EpochStrStart = ['epoch', num2str(currentEpochNumStart)];
      end
      if currentEpochNumEnd<10
        EpochStrEnd = ['-0', num2str(currentEpochNumEnd)];
      else
        EpochStrEnd = ['-', num2str(currentEpochNumEnd)];
      end
      if currentEpochNumStart==currentEpochNumEnd
        EpochStr = EpochStrStart;
      else
        EpochStr = [EpochStrStart, EpochStrEnd];
      end
 
      % open diary for the current epoch of the training
      clc;
      diaryFileName = [netDir, filesep, 'Progress_Report_', ParamUNet.eyeDisease, ImageSizeStr, OverlapStr, '_', FltSizeStr, EncoderDepthStr, InitChanNumStr, EpochStr, '.txt'];
      diary(diaryFileName);

      % display computer name for diary recording
      disp(['computer name: ', ComputerName]);
      disp(['running script name: ', RunningScriptName]);
      fprintf('\n');

      % run one epoch training at a time
      disp(['Data set directory = ', dataDir]);
      disp(['run epoch #', num2str(currentEpochNumStart), ' ...']);
      tic
      if epochInterval==1
        net = RPModelUNetTrainOneEpoch(ds, ds_v, netDir, Param);
      else
        net = RPModelUNetTrainMultiEpochs(ds, ds_v, netDir, Param);
      end
      toc

      diary off;
    end

  end
  
	fprintf('\n');

end

% end of train_RP_model_batch_UNet.m
