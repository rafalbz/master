%This script perform the PIV analysis

clear

totalTime = tic;


javaaddpath('HydrolabPIV-master/src/measures');
javaaddpath('HydrolabPIV-master/src/interp');
addpath(genpath('HydrolabPIV-master/'));


if isempty(gcp('nocreate'))
    parpool('local', 32);
end
pctRunOnAll javaaddpath('HydrolabPIV-master/src/measures');
pctRunOnAll javaaddpath('HydrolabPIV-master/src/interp');


startIndex = 3; % starting index of images start from 3
endIndex = 1503; % number of images to process
stepSize = 1; % step size (1 = consecutive, 2 = every second image, etc.)

filePattern = '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/final/29.03/smooth_1f/run1/run1_%06d.png';
matFilename = '/uio/hypatia/geofag-personlig/students/metos/rafalbz/master/final/29.03/smooth_1f/PIV_run1.mat';


numImages = (endIndex - startIndex) / stepSize + 1;


a1 = 81;
b1 = 243;
opt1 = setpivopt('range', [-a1 a1 -a1 a1], 'subwindow', b1, b1, .75);

a2 = 27;
b2 = 81;
opt2 = setpivopt('range', [-a2 a2 -a2 a2], 'subwindow', b2, b2, .75);


pivResults1 = cell(1, numImages - 1);
pivResults2 = cell(1, numImages - 1);
U1 = cell(1, numImages - 1);
V1 = cell(1, numImages - 1);
U2 = cell(1, numImages - 1);
V2 = cell(1, numImages - 1);


batchSize = 500;
numBatches = ceil((numImages - 1) / batchSize);


for batch = 1:numBatches
    fprintf('Processing batch %d of %d\n', batch, numBatches);
    

    startPairIdx = (batch - 1) * batchSize + 1;
    endPairIdx = min(batch * batchSize, numImages - 1);
    

    startImgIdx = startPairIdx;
    endImgIdx = endPairIdx + 1;
    

    batchImages = cell(1, endImgIdx - startImgIdx + 1);
    for i = startImgIdx:endImgIdx
        imgFileIdx = startIndex + (i - 1) * stepSize;
        batchImages{i - startImgIdx + 1} = imread(sprintf(filePattern, imgFileIdx));
    end
    

    parfor i = startPairIdx:endPairIdx
        try

            localIdx = i - startImgIdx + 1;
            
            im1 = batchImages{localIdx};
            im2 = batchImages{localIdx + 1};
            
            % Perform first pass PIV
            tempRes1 = normalpass([], im1, [], im2, [], opt1);
            [U1{i}, V1{i}] = replaceoutliers(tempRes1);
            pivResults1{i} = tempRes1;
            
            % Perform second pass PIV
            tempRes2 = normalpass(tempRes1, im1, [], im2, [], opt2);
            [U2{i}, V2{i}] = replaceoutliers(tempRes2);
            pivResults2{i} = tempRes2;
            
        catch ME
            warning('Error processing image pair %d: %s', i, ME.message);
        end
    end
    

    clear batchImages;
end


x1 = pivResults1{1}.x;
y1 = pivResults1{1}.y;
x2 = pivResults2{1}.x;
y2 = pivResults2{1}.y;


save(matFilename, 'U2', 'V2', 'x2', 'y2', 'U1', 'V1', 'x1', 'y1', 'pivResults1', 'pivResults2'); %U2, V2, x2, y2 are the final results

disp('PIV analysis and saving complete.');


delete(gcp('nocreate'));


elapsedTotalTime = toc(totalTime);
fprintf('Total time for PIV analysis and saving: %.2f seconds\n', elapsedTotalTime);