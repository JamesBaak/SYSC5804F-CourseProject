train = false;
threshold = 0.3;
networkInputSize = [256 256 3];
numClasses = 1;
allImageFileNames = [];
allDataPaths = {};

if train
    % Load the labeled data
    GroundTruthData = {'spots3.mat','spots4.mat','spots5.mat'};
    ImageFolders    = {'spots_images_3','spots_images_4','spots_images_5'};
    for d = 1:length(GroundTruthData)
        data = load(GroundTruthData{d}).gTruth;
        
        % Generate the appropriate path names
        [filepath,name,ext] = fileparts(string(data.DataSource.Source));
        imageFileNames = strings(length(name),1);
        for k = 1:length(name)
           imageFileNames(k) = fullfile(pwd,ImageFolders(d),append(name(k),ext(k)));
        end

        allImageFileNames = [allImageFileNames;imageFileNames];
        allDataPaths      = [allDataPaths;data.LabelData.path];
    end

    % Create structured table
    dataset = table(allImageFileNames);
    dataset.gTruth = allDataPaths;
    
    % Split Data
    rng(0);
    shuffledIndices = randperm(length(name));
    idx = floor(0.8 * length(shuffledIndices)); % 80% data for training
    trainingDataTbl = dataset(shuffledIndices(1:idx), :);
    testDataTbl = dataset(shuffledIndices(idx+1:end), :);

    % Image Datastore
    imdsTrain = imageDatastore(trainingDataTbl.allImageFileNames);
    imdsTest = imageDatastore(testDataTbl.allImageFileNames);

    % Ground truth Datastore
    bldsTrain = boxLabelDatastore(trainingDataTbl(:, 2:end));
    bldsTest = boxLabelDatastore(testDataTbl(:, 2:end));

    trainingData = combine(imdsTrain, bldsTrain);
    testData = combine(imdsTest, bldsTest);
    % validateInputData(trainingData);
    % validateInputData(testData);

    % Training =================

    % Training Setup ==========
    
     augmentedTrainingData = transform(trainingData, @augmentData);

    % Preprocessing testing ===
    preprocessedTrainingData = transform(augmentedTrainingData, @(data)preprocessData(data, networkInputSize));
    preprocessedValidationData = transform(testData,@(data)preprocessData(data,networkInputSize));

    rng(0)
    trainingDataForEstimation = transform(preprocessedTrainingData, @(data)preprocessData(data, networkInputSize));
    numAnchors = 6;
    [anchorBoxes, meanIoU] = estimateAnchorBoxes(trainingDataForEstimation, numAnchors);
    
    % Create network
    featureExtractionNetwork = darknet53; % Use darknet-53
    featureLayer = 'conv33';
%     featureExtractionNetwork = resnet50('Weights', 'none');
%     featureLayer = 'activation_40_relu';

    lgraph = yolov2Layers(networkInputSize,numClasses,anchorBoxes,featureExtractionNetwork,featureLayer);

    options = trainingOptions('sgdm', ...
            'MiniBatchSize',6, ....
            'InitialLearnRate',1e-3, ...
            'MaxEpochs',20, ... 
            'CheckpointPath',tempdir, ...
            'ValidationData',preprocessedValidationData, ...
            'Plots','training-progress');
        %             'VerboseFrequency',5, ...

    [detector,info] = trainYOLOv2ObjectDetector(preprocessedTrainingData,lgraph,options);
end

testImages = {'spots_images\1.png','spots_images\2.png','spots_images\3.png','spots_images\4.png'};

for ti = 1:length(testImages)
   [img, map] = imread(testImages{ti});
    I = cat(3,img,img,img); % Convert to RGB
    I = imresize(I,networkInputSize(1:2));
    [bboxes,scores] = detect(detector,I,'Threshold', threshold);
    if ~isempty(bboxes) && ~isempty(scores)
        I = insertObjectAnnotation(I,'rectangle',bboxes,scores);
    end
    figure
    imshow(I) 
end

% Testing
preprocessedTestData = transform(testData,@(data)preprocessData(data,networkInputSize));
detectionResults = detect(detector, preprocessedTestData, 'Threshold', threshold);
[ap,recall,precision] = evaluateDetectionPrecision(detectionResults, preprocessedTestData, threshold);

% Plot results
figure
plot(recall,precision)
xlabel('Recall')
ylabel('Precision')
grid on
title(sprintf('Average Precision = %.2f',ap))


