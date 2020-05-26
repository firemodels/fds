% McDermott
% 26 May 2020
% image_compare.m

close all
clear all

fds_repo = '/Users/rmcdermo/GitHub/FireModels_rmcdermo/fds/';
fig_repo = '/Users/rmcdermo/GitHub/FireModels_rmcdermo/fig/';

% test two images that should be the same

I1 = imread([fds_repo,'Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/symmetry_test.png']);
I2 = imread([fig_repo,'fds/Reference_Figures/symmetry_test.png']);

figure
imshow(I1)

figure
imshow(I2)

Diff1 = double(I2(:,:,1)-I1(:,:,1));
Diff2 = double(I2(:,:,2)-I1(:,:,2));
Diff3 = double(I2(:,:,3)-I1(:,:,3));

figure
surf(Diff1)

figure
imshow(Diff1)

norm(Diff1)/numel(Diff1)
norm(Diff2)/numel(Diff2)
norm(Diff3)/numel(Diff3)

% test two images we know to be different
clear all

fds_repo = '/Users/rmcdermo/GitHub/FireModels_rmcdermo/fds/';
fig_repo = '/Users/rmcdermo/GitHub/FireModels_rmcdermo/fig/';

I1 = imread([fds_repo,'Manuals/FDS_Verification_Guide/SCRIPT_FIGURES/symmetry_test.png']);
I2 = imread([fig_repo,'fds/Reference_Figures/ribbed_channel_80_s1000.png']);

% images must be the same size for comparison
% so, take max element size and fill values with 0

I1_SIZE = size(I1)
I2_SIZE = size(I2)

I1_NEW = zeros(max(I1_SIZE(1),I2_SIZE(1)),max(I1_SIZE(2),I2_SIZE(2)),3);
I2_NEW = 255*ones(max(I1_SIZE(1),I2_SIZE(1)),max(I1_SIZE(2),I2_SIZE(2)),3);

% NOTE: Matlab is 1 based for arrays
I1_NEW(1:I1_SIZE(1),1:I1_SIZE(2),:) = I1;
I2_NEW(1:I2_SIZE(1),1:I2_SIZE(2),:) = I2;

figure
imshow(I1_NEW)

figure
imshow(I2_NEW)

figure
imshow(I2_NEW-I1_NEW)

Diff1 = double(I2_NEW(:,:,1)-I1_NEW(:,:,1));
Diff2 = double(I2_NEW(:,:,2)-I1_NEW(:,:,2));
Diff3 = double(I2_NEW(:,:,3)-I1_NEW(:,:,3));

figure
surf(Diff1)

figure
imshow(Diff1)

norm(Diff1)/numel(Diff1)
norm(Diff2)/numel(Diff2)
norm(Diff3)/numel(Diff3)
