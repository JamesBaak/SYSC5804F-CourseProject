% Load the labeled data
data = load('spots.mat').gTruth;

% Generate the appropriate path names
[filepath,name,ext] = fileparts(string(data.DataSource.Source));
imageFileNames = strings(length(name),1);
for k = 1:length(name)
   imageFileNames(k) = fullfile(pwd,'spots_images',append(name(k),ext(k)));
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
networkInputSize = [227 227 3];

% Preprocessing testing ===
% preprocessedTrainingData = transform(augmentedTrainingData, @(data)preprocessData(data, networkInputSize));
% data = read(preprocessedTrainingData);
% I = data{1,1};
% bbox = data{1,2};
% annotatedImage = insertShape(I, 'Rectangle', bbox);
% annotatedImage = imresize(annotatedImage,2);
% figure
% imshow(annotatedImage)
% reset(preprocessedTrainingData);

rng(0)
trainingDataForEstimation = transform(trainingData, @(data)preprocessData(data, networkInputSize));
numAnchors = 6;
[anchors, meanIoU] = estimateAnchorBoxes(trainingDataForEstimation, numAnchors);

area = anchors(:, 1).*anchors(:, 2);
[~, idx] = sort(area, 'descend');
anchors = anchors(idx, :);
anchorBoxes = {anchors(1:3,:)
    anchors(4:6,:)
    };
% Create network
net = darknet53('Weights','None'); % Use darknet-53
classNames = trainingDataTbl.Properties.VariableNames(2:end);
lgraph = load('darknet_yolo.mat').lgraph_2;

options = trainingOptions('sgdm',...
          'InitialLearnRate',0.001,...
          'Verbose',true,...
          'MiniBatchSize',16,...
          'MaxEpochs',30,...
          'Shuffle','never',...
          'VerboseFrequency',30,...
          'CheckpointPath',tempdir);
      
[detector,info] = trainYOLOv2ObjectDetector(trainingData,lgraph,options);
