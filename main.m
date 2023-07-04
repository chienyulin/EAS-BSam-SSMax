% =========================================================================
% This main script is to implement ISSMax, EAS-ISSMax, BSam-SSMax, 
% EAS-BSam, and EAS-BSam-SSMax mentioned in the paper 
% "Band Sampling of Hyperspectral Anomaly Detection in Effective Anomaly Space."
%
% INPUTS:
%   - HIM:  hyperspectral imagecube (rows by columns by bands)
%   - GT:   ground truth (rows by columns)
% OUTPUTS:
%   - final_map:    final detection map (rows by columns)
%   - qdm_report:   quantitative detection measures report
%
% REFERENCE:
%   To-Be-Added
%
% Written and sorted by © 2023 Chien-Yu Lin.
% =========================================================================

clc
clear
close all;
dbstop if error

% find where this main.m file is.
folder = fileparts(which(mfilename)); 
% add the folder plus all subfolders to the path.
addpath(genpath(folder));
% random seed setting for reproducibility
my_rng_seed = 304;
rng(my_rng_seed,'twister')

%% hyperspectral image and parameters setting
load('panelHIM.mat','-mat','HIM')
load('hydice_GT_binary.mat','-mat','GT')

para_set.max_iter = 20;
para_set.e = 0.95;
para_set.loop_rng_seed = 304;

%% ISSMax
[final_map, qdm_report] = ISSMax(HIM,GT,para_set);
qdm_report

%% EAS-ISSMax
[final_map, qdm_report] = EAS_ISSMax(HIM,GT,para_set);
qdm_report

%% BSam-SSMax
[final_map, qdm_report] = BSam_SSMax(HIM,GT,para_set);
qdm_report

%% EAS-BSam
[final_map, qdm_report] = EAS_BSam(HIM,GT,para_set);
qdm_report

%% EAS-BSam-SSMax
[final_map, qdm_report] = EAS_BSam_SSMax(HIM,GT,para_set);
qdm_report

