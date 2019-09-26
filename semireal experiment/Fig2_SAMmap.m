clc; clear; close all;
addpath([pwd,'/data and specifications']);
load('record.mat');
SAM_map = record.SAM_map;
count = 1;
for i = [1:4 6 12]
    img = imgray2pcolor(reshape(SAM_map(i,:,:),1070,1070)/16,'jet',255);
    imshow(img);
    pause(1)
end
