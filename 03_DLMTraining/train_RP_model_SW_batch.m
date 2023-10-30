% train_RP_model_SW_batch.m
%   script to train sliding-window (SW) CNN model with image patch datasets
%   generated from OCTScan*.mat files using script_generate_datasets_SW.m
%
% by YZW, 12/09/2019; additional comments on 10/30/2023
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
epochsToRun   = 5;   % number epochs to run; default: 45
epochInterval = 5;   % interval of epochs to be finalized (e.g., =5 -> every 5 iteration to get a finalized Net), including the first and last epochs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('define parameters for the SW CNN model structure and training ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParamSW                           = DefaultParamSW;
ParamSW.eyeDisease                = 'RP';
ParamSW.epochsToRun               = epochsToRun;
ParamSW.convFilterSize            = [5, 5];
ParamSW.initialEncoderNumChannels = 32;
ParamSW.Verbose                   = 1;
ParamSW.VerboseFrequency          = 100;
ParamSW.validation                = 1;      % validation takes long time for large validataion dataset
ParamSW.numValidImages            = inf;    % number of images in validation set used; = inf -> use all validation images
ParamSW.learnRateDropPeriod       = 10;     % how many epochs to change learting rate
ParamSW.initialLearningRate       = 0.05;   % 0.05 -> match MatConvNet
ParamSW.maxEpochs                 = 1;      % number of epochs to run at a time by TrainOneEpoch or TrainMultiEpochs; =1->epoch to be finalized; otherwise # epochs between finalized epochs
ParamSW.miniBatchSize             = 128;    % 100 -> match MatConvNet; =128 -> match UNet
ParamSW.batchNormalization        = 1;      % =1 -> add batch normalization layer; =0 -> don't add
ParamSW.Plots                     = 'none'; % 'training-progress' or 'none'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('get parent directory of SW datasets and the list of datasets ...');
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
  
  % get dataset name
  Indx = find(dataDir==filesep);
  DataSetName = dataDir(Indx(end)+1:end);
  fprintf('\n');
  disp(['Training CNN using dataset ''', DataSetName, '''']);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('load in metadata and OCTScan file names and update parameters for SW CNN ...');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  load([dataDir, filesep, 'metadata.mat'], 'metadata');
  load([dataDir, filesep, 'All_OCTScan_Files.mat'], 'OCTScanFileNames');    % note that All_OCTScan_Files.mat contains 'OCTScanFileNames' and other variables
  PatchParam  = metadata.Param;
  classNames  = metadata.ClassName;
  classLabels = metadata.ClassLabel;      % ClassLabel starts from 0
  numClasses  = length(classNames);
  
  % update ParamSW
  ParamSW.numClasses  = numClasses;
  ParamSW.classNames  = classNames;
  ParamSW.classLabels = classLabels;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('load training images and labels ...');    % to get image patch size
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  load([dataDir, filesep, 'TrainImages']);
  load([dataDir, filesep, 'TrainLabels']);
  % TrainLabels = TrainLabels + 1;              % index starts from 1

  % get image size and update ParamSW
  [height, width, numChannels, ~] = size(TrainImages);
  imageSize = [height width numChannels];
  ParamSW.imageSize  = imageSize;

  % get the type of eye disease from dataSetDir and update ParamSW
  idx = find(dataDir==filesep);
  dataSetFolderName = dataDir(idx(end)+1:length(dataDir));
  Indx = find(dataSetFolderName=='_');
  eyeDisease = dataSetFolderName(1:Indx(1)-1);
  ParamSW.eyeDisease = eyeDisease;

  % get number of iterations in an epoch and set validateFreq
  numIterationsPerEpoch   = length(TrainLabels)/ParamSW.miniBatchSize;
  ParamSW.validateFreq = floor(numIterationsPerEpoch);

  % from PatchParam.tvRatio to percent training image patches
  tvRatio = PatchParam.tvRatio;
  if tvRatio<inf
    PerCTrain = ['_T0', num2str(round(100 * tvRatio ./ (tvRatio + 1)))];
  else
    PerCTrain = '_T100';
  end

  % strings as part of subfolder/file names - include '_SW_'
  ImageSizeStr    = [PerCTrain, '_SW', num2str(ParamSW.imageSize(1)), 'x', num2str(ParamSW.imageSize(2)), '_'];
  
  % convert Param.overlap in fraction or pixels to overlap in pixels
  % MOD by YZW on 08/02/2021 for overlap>1
  overlap = PatchParam.overlap;
  n = ParamSW.imageSize(2);
  [overlapPixel, OverlapStr] = ConvertOverlapToPixel(overlap, n);

  FltSizeStr      = ['f', num2str(ParamSW.convFilterSize(1))];
  if ParamSW.initialEncoderNumChannels<10
    InitChanNumStr  = ['i0', num2str(ParamSW.initialEncoderNumChannels)];
  else
    InitChanNumStr  = ['i', num2str(ParamSW.initialEncoderNumChannels)];
  end

  % get sub-folder name for SW CNN training results and make the subfolder
  % the naming convention SW CNN training result folder:
  %   [disease_imageSize_filterSize_initialChannels]
  netFolderName = [ParamSW.eyeDisease, ImageSizeStr, OverlapStr, '_', FltSizeStr, InitChanNumStr];
  netDir = [dataDir, filesep, netFolderName];
  if ~isfolder(netDir)
    mkdir(netDir);
  end
  disp(['Data set directory = ', dataDir]);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('check if there is any checkpoint and trained net ...');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  [~, EpochsDone] = GetTrainedNet(netDir);

  % get starting epoch #
  if ~isempty(EpochsDone)
    startEpoch = EpochsDone.numEpochsDone+1;
  else
    startEpoch = 1;
  end

  if startEpoch>ParamSW.epochsToRun
    disp(['all epochs are done for ''', dataSetFolderName, ''';  next dataset']);
  else              % still having epochs to run
    % assemble Param that includes ParamSW and PatchParam
    Param.PatchParam       = PatchParam;
    Param.ModelParam       = ParamSW;
    Param.OCTScanFileNames = OCTScanFileNames;    % added on 12/19/2021; store names of OCTScan files used to generate the training dataset
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('convert TrainLabels to categorical data type for training ...')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TrainLabelsC    = categorical(TrainLabels, classLabels, classNames);
    classCategories = categories(TrainLabelsC);	% check the order of categories

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
      for i=1:length(dataSetList)
        if contains(dataSetList(i).name, '_Test')
          count = count+1;
          IndxTestDataSetDir(count) = i;
        end
      end
      numTestDatasets = length(IndxTestDataSetDir);
      disp(['number of test datasets = ', num2str(numTestDatasets)]);

      testDataDir = [parentDataDir, filesep, dataSetList(IndxTestDataSetDir(1)).name];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('load validation images and labels ...');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist([testDataDir, filesep, 'ValidImages.mat'], 'file') && exist([testDataDir, filesep, 'ValidLabels.mat'], 'file')
      % load testing/validation images and labels
      load([testDataDir, filesep, 'ValidImages']);
      load([testDataDir, filesep, 'ValidLabels']);
      % ValidLabels = ValidLabels + 1;                  % index starts from 1

      if ~isempty(ValidImages) && ~isempty(ValidLabels)
        disp('convert ValidLabels to categorical data type for validation ...')
        ValidLabelsC = categorical(ValidLabels, classLabels, classNames);
        numValidImages = min(ParamSW.numValidImages, length(ValidLabelsC));

        ValidImages  = ValidImages(:,:,1,1:numValidImages);
        ValidLabelsC = ValidLabelsC(1:numValidImages);
      else
        ValidImages  = [];
        ValidLabelsC = [];
      end

    else
      ValidImages  = [];
      ValidLabelsC = [];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('display training image examples in thumbnails and their labels ...');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DisplayPatchThumbnails(TrainImages, TrainLabels, PatchParam, 100, 8);
    title('Examples of Training Image Patches');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('display testing image examples in thumbnails and their labels ...');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(ValidImages)
      DisplayPatchThumbnails(ValidImages, ValidLabels, PatchParam, 100, 8);
      title('Examples of Validation Image Patches');
      pause(0.1);
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
    % loop through the list of remaining epoch groups to train SW model
    % MOD on 09/28/2022 by YZW for group epochs training
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:numEpochRunGroups
      indxToRemainEpochsGroupStart = sum(NumEpochsInGroup(1:i)) - NumEpochsInGroup(i) + 1;
      indxToRemainEpochsGroupEnd   = sum(NumEpochsInGroup(1:i));        % accumulated to the last epoch in the group
      currentEpochNumStart  = RemainEpochs(indxToRemainEpochsGroupStart);
      currentEpochNumEnd    = RemainEpochs(indxToRemainEpochsGroupEnd);  % currentEpochNum is the last Epoch in the group

      % update epoch number for the current model training
      Param.ModelParam.currentEpochNum = currentEpochNumEnd;    % this is the number for the finalized epoch in the group

      % update .maxEpochs in ParamSW or Param.ModelParam with NumEpochsInGroup(i)
      ParamSW.maxEpochs          = NumEpochsInGroup(i);
      Param.ModelParam.maxEpochs = NumEpochsInGroup(i);

      % only do validation for odd epochs
      % note for multi-epoch groups, validataion only for the groups ended
      % with odd epoch number -> epochs in the same group either with or
      % without validation
      % also make sure to run validation for the last epoch group, 
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
      diaryFileName = [netDir, filesep, 'Progress_Report_', ParamSW.eyeDisease, ImageSizeStr, OverlapStr, '_', FltSizeStr, InitChanNumStr, EpochStr, '.txt'];
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
        net = RPModelSWTrainOneEpoch(TrainImages, TrainLabelsC, ValidImages, ValidLabelsC, netDir, Param);
      else
        net = RPModelSWTrainMultiEpochs(TrainImages, TrainLabelsC, ValidImages, ValidLabelsC, netDir, Param);
      end
      toc

      diary off;
    end

  end
  
  fprintf('\n');
  
end

% end of train_RP_model_SW_batch.m
