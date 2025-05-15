% Transformation from image coord. to real-world coord.
%{
% Add necessary paths
javaaddpath('C:/Users/Rafal/Documents/Master/main/code/HydrolabPIV-master//src/measures');
javaaddpath('C:/Users/Rafal/Documents/Master/main/code/HydrolabPIV-master/src/interp');
addpath(genpath('C:/Users/Rafal/Documents/Master/main/code/HydrolabPIV-master/'));
addpath(genpath('C:/Users/Rafal/Documents/Master/code'));


%%{
%coord = imread('C:/Users/Rafal/Documents/Master/piv_source/final/28.03/scaling/scaling/scaling_000175.png');
coord = imread('C:/Users/Rafal/Documents/Master/piv_source/final/29.03/scaling/scaling/scaling_000166.png');
imshow(coord);
axis equal;
hold on;
imagesc(coord);

h = impoly;
pixel = h.getPosition;

c = graythresh(coord);
bw = im2bw(coord);
cc = bwconncomp(bw);
stats = regionprops(cc, 'Centroid');
xc = vertcat(stats.Centroid);
idx = knnsearch(xc, pixel);
pixel = xc(idx, :);

%save('C:/Users/Rafal/Documents/Master/piv_source/final/coord_28.mat', 'pixel');
save('C:/Users/Rafal/Documents/Master/piv_source/final/coord_29.mat', 'pixel');
%}
 


%{
coord = imread('C:/Users/Rafal/Documents/Master/piv_source/final/28.03/scaling/canyon/canyon.png');
imshow(coord);
axis equal;
hold on;


% Load existing points
load('C:/Users/Rafal/Documents/Master/piv_source/final/coord_28.mat', 'pixel');

% Plot using impoly
h = impoly(gca, pixel);

% Wait for user to finish (use h.getPosition directly)
uiwait(msgbox('Adjust points, then click OK to save.'));  % Message box to pause execution

pixel = h.getPosition();  % Get updated position
save('C:/Users/Rafal/Documents/Master/piv_source/final/coord_28.mat', 'pixel');  % Save the new position
%disp('Updated points saved to file.');
%}

%{
coord = imread('C:/Users/Rafal/Documents/Master/piv_source/final/29.03/scaling/canyon\canyon.png');
imshow(coord);
axis equal;
hold on;


% Load existing points
load('C:/Users/Rafal/Documents/Master/piv_source/final/coord_29.mat', 'pixel');

% Plot using impoly
h = impoly(gca, pixel);

% Wait for user to finish (use h.getPosition directly)
uiwait(msgbox('Adjust points, then click OK to save.'));  % Message box to pause execution

pixel = h.getPosition();  % Get updated position
save('C:/Users/Rafal/Documents/Master/piv_source/final/coord_29.mat', 'pixel');  % Save the new position
%disp('Updated points saved to file.');
%}