%% setup data files 

%% save data for easier retrieval + changing
clear all; close all; clc;
dat = struct;
dat.bins = 30:20:290; % center of bins, first bin is [20,40], last [280,300]
dat.f300 = [5, 4, 5, 5, 5, 14, 25, 9, 2, 2, 1, 1, 1, 0];
dat.f3000 = [8, 7, 8, 5, 3, 7, 11, 27, 9, 3, 2, 0, 0, 0];
dat.f30000 = [4, 5, 4, 4, 6, 5, 8, 8, 20, 11, 4, 1, 2, 2];
dat.ctrl = [9, 6, 6, 4, 3, 3, 3, 3, 2, 1, 1, 2, 1, 1];
dat.f = [300,3000,30000,0];

save('dat.mat','dat')