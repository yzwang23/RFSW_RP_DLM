% script_generate_datasets_UNet.m
%   script to extract image patches and their classes from a single set of 
%   OCTScan*.mat files for DML model training and validation, and save them
%   into datasets in specified folders. This script will generate datasets 
%   in .mat files for UNet. Multiple datasets can be generated based on %
%   patients included. Note that the validation set generated using this
%   script is for internal validation process. While there is no
%   overlapping between training and validation image patches, the patches
%   used for validation are likely from the patients whose data are also
%   used to generate training patches.
%
% Yi-Zhong Wang, 06/02/2020
% Retina Foundation of the Southwest
% Copyright. All Rights Reserved.
%

% Note that this script is based on script_generate_datasets_oneset.m
% Note that a fixed random sequence is used for the list of patients by
%	employing the following method to divide patients into percent
% 	rng('default');
% 	A = RandomSequence(numPatients, [0 1], 1);

clc; clearvars; close all;

% save scaled patch images or not
saveScaledPatchImage = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('define parameters - common and specific to the model ...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% unet specific parameters
mUNet           = 256;    	% [mUNet, nUNet] -> patch size
nUNet           = 32;
overlapUNet     = 0.875;    % patch overlap, e.g., 0.50 -> 50% overlapping
patchFormatUNet = 1;        % default patchFormat for UNet is 1

% training/validation ratio
tvRatio = 4;              % training/validation ratio = 4 -> 80/20 split; inf -> all for training

% display sample patches or not
checkSamplePatches = 1;

% datasets generated with percent # of patients from the full training data
% PerNum = [10, 20, 30, 40, 50, 60, 70, 80, 90. 100];
PerNum = 100;

% training mini-batch size
miniBatchSize = 128;

% specify OCT segmentation boundary lines defined by gold standard, as well
% as define class names and class labels, including 'NONE' as the first 
% element with class 0 (negative patch or background); Note that for area
% segmentation (unet), number of classes is equal to # of boundary lines
LineName       = {'ILM' 'INL' 'PR1' 'RPE' 'BM'};         % segment lines
ClassNameUNet	 = {'None' 'ILM_INL','INL_PR1','PR1_RPE','RPE_BM'};
ClassLabelUNet = [0 1 2 3 4];

% make sure the order of line names from top to bottom
LineName = HEOCTCheckLineNameListOrder(LineName);	% check order

% list of graders - for gold-standard classifications of patches
% Graders = {'YZW', 'DG'};      % using labels from multiple graders' 
Graders = {'YZW'};

% from tvRatio to percent training
if tvRatio<inf
  PerCTrain = ['_0', num2str(round(100 * tvRatio ./ (tvRatio + 1)))];
else
  PerCTrain = '_100';
end
                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('specify a parent directory to save the datasets generated ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parentDataDir=[getenv('HOME'),'/Documents/DML_Models/RPModelUNet/Data'];
parentDataDir = uigetdir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('specify a directory of OCTScan files for training/validation ...');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify a directory that contains OCTScan*.mat files used to generate
% training images and labels, and get a list of file names and other info
% note that here new naming convention for OCTScan*.mat files is used
D = GetOCTScanFiles;                % scanType in the 2nd place, ID 4th

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate datasets for UNet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('generate datasets for UNet ...');

% modify parameters specific to unet
Param            = DefaultParamPatchExtraction;
Param.segType    = 0;       % =0 -> area; =1 -> line
Param.m          = mUNet;
Param.n        	 = nUNet;
Param.LineName   = LineName;
Param.ClassName  = ClassNameUNet;
Param.ClassLabel = ClassLabelUNet;
Param.Graders    = Graders;
Param.overlap    = overlapUNet;
Param.tvRatio    = tvRatio;
Param.batchSize  = miniBatchSize;

% pre-process or not
Param.Pre.preProcess = 0;

% make sub-data directory to save datasets for unet training
DirName = D.DirName;            % training OCTScan files folder
Indx = find(DirName==filesep);
SubDataDir = DirName(Indx(end)+1:length(DirName));
NameApend = [PerCTrain, '_UNET_', num2str(Param.m), '_', num2str(Param.n), '_O', num2str(round(Param.overlap * Param.n))];
SubDataDirUNet = [SubDataDir, NameApend];
SubDataDirUNetFull = [parentDataDir, filesep, SubDataDirUNet];
mkdir(SubDataDirUNetFull);


% open diary to record datasets generated
clc;
diaryFileName = [SubDataDirUNetFull, filesep, SubDataDirUNet, '_Datasets_Generated.txt'];
diary(diaryFileName);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('loop through percent # of patients to generate datasets of image patches ...');
rng('default');     % fixed random sequence for the list of patients
RandPtOrder = RandomSequence(length(D.UniqueID)); % random sequence for patients selection
Param.patchType   = 1;                % training patches
Param.patchFormat = patchFormatUNet;  % the way training patch generated
Param.bgPatches 	= 0;                % including background patches or not
Param.overlap     = overlapUNet;      % restore overlap value
for p=1:length(PerNum)
  tic;
  % get number of patients for generating image patches
  fprintf('\n');
  disp(['  get ', num2str(PerNum(p)), '% of patients from the list of OCTScan files for generating image patches ...']);
  numPts = round(length(D.UniqueID)*PerNum(p)/100);
  PtIDs = D.UniqueID(RandPtOrder(1:numPts));

  % search for PtIDs in D.ID to get OCTFileNames
  Indx             = contains(D.ID, PtIDs);
  OCTScanFileNames = D.FileNames(Indx);
  ID               = D.ID(Indx);          % for saving
  Eye              = D.Eye(Indx);         % for saving

  disp('  extract image patches & their labels for unet training ...');
  [Patches, PatchesC] = ExtractPatchesFromOCTScans(OCTScanFileNames, DirName, Param);
  numPatches = size(Patches, 4);
  disp(['  ', num2str(numPatches), ' total number of patches generated from ', num2str(PerNum(p)), '% of patients ...']);
  
  % divide Patches into training & validating sets according to tvRatio
  batchSize	 = Param.batchSize;
  numBatches = floor(numPatches./batchSize);
  tvRatio    = Param.tvRatio;
  if tvRatio==inf
    nTrainBatch = numBatches;
    nTestBatch  = 0;
  else
    nTrainBatch = floor(numBatches * tvRatio ./ (tvRatio+1));
    nTestBatch  = numBatches - nTrainBatch;
  end
  numTrainPatches = nTrainBatch * batchSize;
  numTestPatches  = nTestBatch * batchSize;

  % get training and validation patches and their labels - only one channel
  % using the same numTrainPatches & numTestBatches as for MatConvNet
  TrainImages = Patches(:, :, 1, 1:numTrainPatches);
  TrainLabels = PatchesC(:, :, 1, 1:numTrainPatches);
  ValidImages = Patches(:, :, 1, numTrainPatches+1:numTrainPatches+numTestPatches);
  ValidLabels = PatchesC(:, :, 1, numTrainPatches+1:numTrainPatches+numTestPatches);

  % exclude all background patches and their labels for validation images
  % if it is not already done in 'ExtractPatchesFromOCTScans'
  if Param.bgPatches==1
    IndxP = zeros(numTestPatches, 1);
    for k=1:numTestPatches
      PatchLabel = ValidLabels(:, :, 1, k);
      IndxP_temp = find(PatchLabel>0, 1);
      if ~isempty(IndxP_temp)
        IndxP(k) = 1;
      end
    end

    % extract non-background patches and their labels from validation dataset
    Indx_p = find(IndxP>0);					% patches having positive pixel classes
    Patches_p = ValidImages(:, :, :, Indx_p);
    PatchesC_p = ValidLabels(:, :, :, Indx_p);

    % reassign non-background patches and labels back to ValidImages & ValidLabels
    ValidImages = Patches_p;
    ValidLabels = PatchesC_p;
  end
  
  % actual number of patches for training and validation
	disp(['  ', num2str(numTrainPatches), ' actual number of training patches generated']);
	disp(['  ', num2str(size(ValidLabels,4)), ' actual number of validation patches generated']);
  
	if checkSamplePatches==1 && size(TrainLabels,4)>0
    disp('  display training image patch examples in thumbnails ...');
    DisplayPatchThumbnails(TrainImages, TrainLabels, Param, 9);
  end

  if checkSamplePatches==1 && size(ValidLabels,4)>0
    disp('  display validatoin image patch examples in thumbnails ...');
    DisplayPatchThumbnails(ValidImages, ValidLabels, Param, 9);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  disp('  save meta data, training/validation image and label patches ...');
  % save patches as .jpg files in SubSubDataDir under SubDataDir
  %   image patches are saved in the subfolder 'TrainingImages
  %     image file naming convention: image_0001.jpg, ...
  %   class patches are saved in the subfolder 'TrainingLabels
  %     label file naming convention: label_0001.jpg, ...

  % for saving files
  numDigits = length(num2str(numPatches));  % # digits in file name

  % add percent patches to the directory name then create sub-sub-directory
  % leading zeros for percent -> so the dataset names have the order needed
  if PerNum(p)<10
    SubSubDataDirUNet	= [SubDataDirUNet, '_00', num2str(PerNum(p)), 'PCT'];
  elseif PerNum(p)<100
    SubSubDataDirUNet	= [SubDataDirUNet, '_0', num2str(PerNum(p)), 'PCT'];
  else
    SubSubDataDirUNet	= [SubDataDirUNet, '_', num2str(PerNum(p)), 'PCT'];
  end
  SubSubDataDirUNetFull  = [parentDataDir, filesep, SubDataDirUNet, filesep, SubSubDataDirUNet];   % path for saving
  mkdir(SubSubDataDirUNetFull);

  % save all names of OCTScan files used for extracting training patches
  OCTScanFilesOutputName = 'All_OCTScan_Files.mat';
  save(fullfile(SubSubDataDirUNetFull, OCTScanFilesOutputName), 'OCTScanFileNames', 'ID', 'Eye');

  % save class names and labels, as well as Param to 'metadata.mat'
  label_names = Param.ClassName';
  metadata.ClassName  = Param.ClassName;
  metadata.ClassLabel = Param.ClassLabel;
  metadata.Param      = Param;
  save(fullfile(SubSubDataDirUNetFull, 'metadata.mat'), 'metadata');
  
  
%   % save trianing image patches & labels arrays in .mat files
%   save(fullfile(SubSubDataDirUNetFull, 'TrainImages.mat'), 'TrainImages', '-v7.3');
%   save(fullfile(SubSubDataDirUNetFull, 'TrainLabels.mat'), 'TrainLabels', '-v7.3');
%   
%   % save validation image patches & labels arrays in .mat files  
%   save([SubSubDataDirSWMATLABFull, filesep, 'ValidateImages.mat'], 'ValidImages', '-v7.3');
%   save([SubSubDataDirSWMATLABFull, filesep, 'ValidateLabels.mat'], 'ValidLabels', '-v7.3');


  % save training and validating image patches & label arrays in files
  SavePatchesToFile(SubSubDataDirUNetFull, 'TrainingImages', TrainImages, 'Image');
  SavePatchesToFile(SubSubDataDirUNetFull, 'TrainingLabels', TrainLabels, 'Label');

  if saveScaledPatchImage~=0
    TrainLabels_scaled = uint8(255 * double(TrainLabels) / max(Param.ClassLabel));
    SavePatchesToFile(SubSubDataDirUNetFull, 'TrainingLabelsScale', TrainLabels_scaled, 'LabelScale');
  end

  % save validation image patches and labels in files for unet
  SavePatchesToFile(SubSubDataDirUNetFull, 'ValidateImages', ValidImages, 'Image');
  SavePatchesToFile(SubSubDataDirUNetFull, 'ValidateLabels', ValidLabels, 'Label');

  if saveScaledPatchImage~=0
    ValidLabels_scaled = uint8(255 * double(ValidLabels) / max(Param.ClassLabel));
    SavePatchesToFile(SubSubDataDirUNetFull, 'ValidateLabelsScale', ValidLabels_scaled(:,:,1,:), 'LabelScale');
  end
  
  toc;

  % clear large variable
  clear Patches PatchesC Patches_p PatchesC_p TrainImages TrainLabels ValidImages ValidLabels;

end

% close diary
fprintf('\n');
fprintf('\n');
diary off;

% end of script_generate_datasets_UNet.m
