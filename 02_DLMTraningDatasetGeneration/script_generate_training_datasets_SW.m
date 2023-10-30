% script_generate_training_datasets_SW.m
%   script to extract image patches as well as their classes from a set of 
%   OCTScan*.mat files containg manual segmentation (ground truth) of 
%   retinal layer boundary lines for deep learning model training and 
%   validation. The image patches are saved into datasets in specified 
%   folders. This script will generate datasets in .mat files for 
%   sliding-window (SW) CNN model. Multiple datasets can be generated based 
%   on percent patients included. Note that the validation set generated 
%   using this script is for internal validation process. While there is no
%   overlapping between training and validation image patches, the patches
%   used for validation are likely from the patients whose data are also
%   used to generate training patches.
%
% Yi-Zhong Wang, 01/30/2019; modified on 05/28/2020; new comments 10/29/23
% Retina Foundation of the Southwest
% Copyright. All Rights Reserved.
%

clc; clearvars; close all;

% make sure the current directory is accessible so that uigetdir works
cd(fileparts(mfilename('fullpath')));

% added on 05/28/2020
dataForMatConvNet = 0;      % = 1 -> generating data for MatConvNet
dataForSWMATLAB   = 1;      % = 1 -> generating data for SW MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('define parameters - common and specific to the model ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tvRatio   = 4;        % training/validation ratio = 4 -> 80/20 split; inf -> all for training
batchSize = 128;      % minibatch size for training

% list of graders - for recording only. All manual graders' segmentation in
% an OCTScan will be used to generate training/validation data
Graders = {'YZW'};    % can add other graders

% generate patches from percent number of patients in the list
% PerNum = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
PerNum = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('specify a parent directory to save the datasets generated ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parentDataDir = uigetdir;
if parentDataDir==0
  return;
end

% retinal layer boundary lines manually graded
LineName = {'ILM' 'INL' 'PR1' 'RPE' 'BM'};

% make sure LineName list is in the correct order (from top to bottom)
LineName = HEOCTCheckLineNameListOrder(LineName);

% get the list of class names, which includes 'NONE' as the first element
% to represent class 0 - negative patches; need to add '1' to patch classes
% for addressing the list of ClassName
ClassName = cat(2, {'NONE'}, LineName);
ClassLabel = 0:length(LineName);          % class labels 0 to # lines

% additional negative patch segment lines other than LineName
NegativeLN = {};

% pre-process or not
preProcess = 0;

% default parameters for generating image patches for SW model training
Param = DefaultParamPatchExtractionSW;

% modify parameters
Param.patchType	     = 1;           % training patches
Param.LineName 	     = LineName;    % positive patch segment lines
Param.ClassName      = ClassName;
Param.ClassLabel     = ClassLabel;
Param.NegativeLN     = NegativeLN;  % additional negative patch segment line
Param.Graders        = Graders;
Param.tvRatio        = tvRatio;   	% 4->training/validating ratio 80/20; inf->all patches for training
Param.batchSize      = batchSize;   % size of mini-batch for training
Param.Pre.preProcess = preProcess;

% convert Param.overlap in fraction or pixels to overlap in pixels
overlap                    = Param.overlap;
n                          = Param.n;
[overlapPixel, OverlapStr] = ConvertOverlapToPixel(overlap, n);

% from tvRatio to percent training
if tvRatio<inf
  PerCTrain = ['_0', num2str(round(100 * tvRatio ./ (tvRatio + 1)))];
else
  PerCTrain = '_100';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('specify a folder of OCTScan files for training/validation ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify a directory that contains OCTScan*.mat files used to generate
% training images and labels, and get a list of file names and other info
% note that here new naming convention for OCTScan*.mat files is used, i.e.
% scanType in the 2nd place, ID in 4th place
D = GetOCTScanFiles;
if isempty(D)
  return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('make sub-data directories for SW MatConvNet and SW MATLAB ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DirName    = D.DirName;                % directory to OCTScan*.mat files
Indx       = find(DirName==filesep);
SubDataDir = DirName(Indx(end)+1:length(DirName));

% check if SubDataDir name contains '_Train' or not, if not, add '_Train'
if ~contains(SubDataDir, '_Train')
  SubDataDir = [SubDataDir, '_Train'];
end

NameApend  = [PerCTrain, '_SW_', num2str(Param.m), '_', num2str(Param.n), '_', OverlapStr];
SubDataDirSWMatConvNet     = [SubDataDir, NameApend, '_MatConvNet'];
SubDataDirSWMATLAB         = [SubDataDir, NameApend];
SubDataDirSWMatConvNetFull = [parentDataDir, filesep, SubDataDirSWMatConvNet];
SubDataDirSWMATLABFull     = [parentDataDir, filesep, SubDataDirSWMATLAB];

if dataForMatConvNet~=0
  mkdir(SubDataDirSWMatConvNetFull);
end

if dataForSWMATLAB~=0
  mkdir(SubDataDirSWMATLABFull);
end


% open diary to record datasets generated
clc;
diaryFileName = [SubDataDirSWMATLABFull, filesep, SubDataDirSWMATLAB, '_Datasets_Generated.txt'];
diary(diaryFileName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('loop through % patients to generate training/validating data ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% default random sequence for patients selection
rng('default');   % fixed random sequence for the list of patients
RandPtOrder = RandomSequence(length(D.UniqueID));
for p=1:length(PerNum)
  tic;
	% get number of patients for generating image patches
  fprintf('\n');
  disp(['  get ', num2str(PerNum(p)), '% of patients from the list of OCTScan files for generating image patches ...']);
  numPts = round(length(D.UniqueID)*PerNum(p)/100);
  PtIDs  = D.UniqueID(RandPtOrder(1:numPts));
  
  % search for PtIDs in D.ID to get OCTFileNames
  Indx             = contains(D.ID, PtIDs);
  OCTScanFileNames = D.FileNames(Indx);
  ID               = D.ID(Indx);
  Eye              = D.Eye(Indx);

  disp('loop through all OCTScan*.mat files to create patches with labels ...');
  Patches  = [];                  % array for patches
  PatchesC = [];                  % array for corresponding labels
  for j=1:length(OCTScanFileNames)    % loop through all training OCTScan*.mat
    % get OCTScan*.mat file name & load OCTScan variable
    FileName = OCTScanFileNames{j};
    load([D.DirName, filesep, FileName]);	% the variable loaded is OCTScan

    % loop through all B scans in OCTScan to extract image patches
    numBscans = OCTScan.numImages - 1;    % first image is infrared
    for k=1:numBscans
      BScanImage = OCTScan.Images(k+1);

      % get ImageArray and sort different types of segmentations
      ImageArray = BScanImage.ImageArray;
      Segments   = HEOCTSortOCTScanSegmentation(BScanImage, FileName);

      % get B-scan image scan parameters
      Param.BScanParam.ScaleXY      = BScanImage.ScaleXY;
      Param.BScanParam.StartCoordXY = BScanImage.StartCoordXY;
      Param.BScanParam.EndCoordXY   = BScanImage.EndCoordXY;

      if ~isempty(Segments.Manual)    	% manual correction - ground truth
        for i=1:length(Segments.Manual)	% could have multiple graders
          if ~isempty(contains(Graders, Segments.Manual(i).Grader))
            % get Segment and create patches
            Segment = Segments.Manual(i).Segment;
            [Patch, PatchC] = GetBscanImagePatchesSW(ImageArray, Segment, Param);

            % note that Patch is a 4-D array, reshape it to 2-D array by
            % keeping the last dimension, transpose row and column before
            % reshape, then transpose reshaped array
            l = size(Patch, 4);
            Patch_r  = reshape(permute(Patch, [2 1 3 4]), [], l)';
            PatchC_r = PatchC';

            % concatenate patches and classes
            Patches  = [Patches; Patch_r];
            PatchesC = [PatchesC; PatchC_r];

            if Param.dataAug==1
              % data augmentation - horizontal flip
              Patch_f    = fliplr(Patch);
              PatchC_f   = PatchC;
              Patch_f_r  = reshape(permute(Patch_f, [2 1 3 4]), [], l)';
              PatchC_f_r = PatchC_f';

              Patches  = [Patches; Patch_f_r];
              PatchesC = [PatchesC; PatchC_f_r];
            end
          end
        end
      end
      
    end 
  end

  % randomize traing patches
  numPatches = length(PatchesC);            % total number of patches
  disp([num2str(numPatches), ' total number of patches generated for ', num2str(PerNum(p)), '% files ...']);
  if numPatches>0
    % note that Patches is a 2-D array
    RandSeq    = RandomSequence(numPatches);	% generate random sequence
    Patches    = Patches(RandSeq, :);         % randomize patches
    PatchesC   = PatchesC(RandSeq);           % same randomization for class
    
    % reshape Patches from 2-D to 4-D for SW MATLAB, keep last dimension,
    % transpose row & col before reshape, then permute reshaped array
    Patches_r = permute(reshape(Patches', 33, 33, 3, []), [2 1 3 4]);
    PatchesC_r = PatchesC';

    % divide Patches into training and validation sets according to .tvRatio,
    batchSize   = Param.batchSize;
    numBatches  = floor(numPatches./batchSize);
    tvRatio     = Param.tvRatio;
    if tvRatio==inf
      nTrainBatch = numBatches;
      nTestBatch  = 0;
    else
      nTrainBatch = floor(numBatches * tvRatio ./ (tvRatio+1));
      nTestBatch  = numBatches - nTrainBatch;
    end
    numTrainPatches = nTrainBatch * batchSize;
    numValidPatches  = nTestBatch * batchSize;
    disp([num2str(numTrainPatches), ' patches for training']);
    disp([num2str(numValidPatches), ' patches for validation']);

    % display example patches
    if Param.checkPatches==1
      disp('display examples of extracted patches in thumbnails ...');
      DisplayPatchThumbnails(Patches_r, PatchesC_r, Param, 100, 8);      
    end

    % get sub-directory names; add % patches to the sub-directory name; add
    % leading zero to % when needed so dataset names have correct order
    if PerNum(p)<10
      SubSubDataDirSWMatConvNet	= [SubDataDirSWMatConvNet, '_00', num2str(PerNum(p)), 'PCT_batches'];
      SubSubExpDirSWMatConvNet = [SubDataDirSWMatConvNet, '_00', num2str(PerNum(p)), 'PCT_lenet'];
      SubSubDataDirSWMATLAB = [SubDataDirSWMATLAB, '_00', num2str(PerNum(p)), 'PCT'];
    elseif PerNum(p)<100
      SubSubDataDirSWMatConvNet	= [SubDataDirSWMatConvNet, '_0', num2str(PerNum(p)), 'PCT_batches'];
      SubSubExpDirSWMatConvNet = [SubDataDirSWMatConvNet, '_0', num2str(PerNum(p)), 'PCT_lenet'];
      SubSubDataDirSWMATLAB = [SubDataDirSWMATLAB, '_0', num2str(PerNum(p)), 'PCT'];
    else
      SubSubDataDirSWMatConvNet	= [SubDataDirSWMatConvNet, '_', num2str(PerNum(p)), 'PCT_batches'];
      SubSubExpDirSWMatConvNet = [SubDataDirSWMatConvNet, '_', num2str(PerNum(p)), 'PCT_lenet'];
      SubSubDataDirSWMATLAB = [SubDataDirSWMATLAB, '_', num2str(PerNum(p)), 'PCT'];
    end
    SubSubDataDirSWMatConvNetFull = [parentDataDir, filesep, SubDataDirSWMatConvNet, filesep, SubSubDataDirSWMatConvNet];
    SubSubExpDirSWMatConvNetFull  = [parentDataDir, filesep, SubDataDirSWMatConvNet, filesep, SubSubExpDirSWMatConvNet];
    SubSubDataDirSWMATLABFull     = [parentDataDir, filesep, SubDataDirSWMATLAB, filesep, SubSubDataDirSWMATLAB];
    
    % save dataset for MatConvNet
    if dataForMatConvNet ~=0
      disp('save dataset batches in the format that MatConvNet uses ...');
      mkdir(SubSubDataDirSWMatConvNetFull);   % make dataset sub-directory
      mkdir(SubSubExpDirSWMatConvNetFull);    % make exp sub-directory

      % save all OCTScan*.mat file names
      OCTScanFilesOutputName = 'All_OCTScan_Files.mat';
      save([SubSubDataDirSWMatConvNetFull, filesep, OCTScanFilesOutputName], 'OCTScanFileNames', 'ID', 'Eye');

      % save class names - label_names and Param to 'batches.meta.mat'
      label_names = Param.ClassName';
      BatchMetaFileName = 'batches.meta.mat';
      save([SubSubDataDirSWMatConvNetFull, filesep, BatchMetaFileName], 'label_names', 'Param');

      % save training & validating batches following MatConvNet convention
      % note that Patches is a 2-D array
      for i=1:numBatches
        Indx = (1+(i-1)*batchSize) : (i*batchSize);
        data = uint8(Patches(Indx, :));
        labels = uint8(PatchesC(Indx));
        if i<=nTrainBatch
          batch_label=['training batch ',num2str(i),' of ',num2str(nTrainBatch)];
          BatchFN = ['data_batch_', num2str(i), '.mat'];
        else
          j = i-nTrainBatch;
          batch_label=['testing batch ',num2str(j),' of ',num2str(nTestBatch)];
          BatchFN = ['test_batch_', num2str(j), '.mat'];
        end

        % save batch file
        save([SubSubDataDirSWMatConvNetFull, filesep, BatchFN], 'data', 'labels', 'batch_label');
      end
    end
    
    % save dataset for SW MATLAB
    if dataForSWMATLAB~=0
      disp('save dataset batches in the format that MATLAB uses ...');
      mkdir(SubSubDataDirSWMATLABFull);     % make dataset sub-directory
      
      % get training and validation patches and their labels
      % using the same numTrainPatches & numTestBatches as for MatConvNet
      TrainImages = Patches_r(:, :, 1, 1:numTrainPatches);
      TrainLabels = PatchesC_r(1:numTrainPatches);
      ValidImages = Patches_r(:, :, 1, numTrainPatches+1:numTrainPatches+numValidPatches);
      ValidLabels = PatchesC_r(numTrainPatches+1:numTrainPatches+numValidPatches);

      % save all OCTScan*.mat file names
      OCTScanFilesOutputName = 'All_OCTScan_Files.mat';
      save([SubSubDataDirSWMATLABFull, filesep, OCTScanFilesOutputName], 'OCTScanFileNames', 'ID', 'Eye');

      % save class names and labels, as well as Param to 'metadata.mat'
      metadata.ClassName  = Param.ClassName;
      metadata.ClassLabel = Param.ClassLabel;
      metadata.Param      = Param;
      save([SubSubDataDirSWMATLABFull, filesep, 'metadata.mat'], 'metadata');
      
      % save trianing image patches & labels arrays in .mat files
      save([SubSubDataDirSWMATLABFull, filesep, 'TrainImages.mat'], 'TrainImages', '-v7.3');
      save([SubSubDataDirSWMATLABFull, filesep, 'TrainLabels.mat'], 'TrainLabels', '-v7.3');
      
     	% save valiation/testing image patches & labels arrays in .mat files
      save([SubSubDataDirSWMATLABFull, filesep, 'ValidImages.mat'], 'ValidImages', '-v7.3');
      save([SubSubDataDirSWMATLABFull, filesep, 'ValidLabels.mat'], 'ValidLabels', '-v7.3');
    end

  end

  toc;

  % clear large variable
  clear Patch Patch_r PatchC PatchC_r Patches Patches_r PatchesC PatchesC_r RandSeq TrainImages TrainLabels ValidImages ValidLabels;

end

% close diary
fprintf('\n');
fprintf('\n');
diary off;

% end of script_generate_training_datasets_SW.m
