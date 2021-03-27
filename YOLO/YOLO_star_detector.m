% Load the labeled data
data = load('spots2.mat').gTruth;

% Generate the appropriate path names
[filepath,name,ext] = fileparts(string(data.DataSource.Source));
imageFileNames = strings(length(name),1);
for k = 1:length(name)
   imageFileNames(k) = fullfile(pwd,'spots_images_2',append(name(k),ext(k)));
end

% Create structured table
dataset = table(imageFileNames);
dataset.gTruth = data.LabelData.path;

% Split Data
rng(0);
shuffledIndices = randperm(length(name));
idx = floor(0.8 * length(shuffledIndices)); % 80% data for training
trainingDataTbl = dataset(shuffledIndices(1:idx), :);
testDataTbl = dataset(shuffledIndices(idx+1:end), :);

% Image Datastore
networkInputSize = [256 256 3];
imdsTrain = imageDatastore(trainingDataTbl.imageFileNames);
imdsTest = imageDatastore(testDataTbl.imageFileNames);

% Ground truth Datastore
bldsTrain = boxLabelDatastore(trainingDataTbl(:, 2:end));
bldsTest = boxLabelDatastore(testDataTbl(:, 2:end));

trainingData = combine(imdsTrain, bldsTrain);
testData = combine(imdsTest, bldsTest);
% validateInputData(trainingData);
% validateInputData(testData);

% Training =================
augmentedTrainingData = transform(trainingData, @augmentData);

% Visualize the augmented images.
% augmentedData = cell(4,1);
% for k = 1:4
%     data = read(augmentedTrainingData);
%     augmentedData{k} = insertShape(data{1,1}, 'Rectangle', data{1,2});
%     reset(augmentedTrainingData);
% end
% figure
% montage(augmentedData, 'BorderSize', 10)

% Training Setup ==========

% Preprocessing testing ===
preprocessedTrainingData = transform(augmentedTrainingData, @(data)preprocessData(data, networkInputSize));
preprocessedValidationData = transform(testData,@(data)preprocessData(data,networkInputSize));

rng(0)
trainingDataForEstimation = transform(preprocessedTrainingData, @(data)preprocessData(data, networkInputSize));
numAnchors = 6;
[anchorBoxes, meanIoU] = estimateAnchorBoxes(trainingDataForEstimation, numAnchors);

% Create network
featureExtractionNetwork = darknet53; % Use darknet-53
featureLayer = 'conv52';
% featureExtractionNetwork = resnet50;
% featureLayer = 'activation_40_relu';

classNames = trainingDataTbl.Properties.VariableNames(2:end);
numClasses = 1;

lgraph = yolov2Layers(networkInputSize,numClasses,anchorBoxes,featureExtractionNetwork,featureLayer);

options = trainingOptions('sgdm', ...
        'MiniBatchSize',8, ....
        'InitialLearnRate',1e-3, ...
        'MaxEpochs',20, ... 
        'CheckpointPath',tempdir, ...
        'ValidationData',preprocessedValidationData);
      
[detector,info] = trainYOLOv2ObjectDetector(preprocessedTrainingData,lgraph,options);
