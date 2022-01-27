% Kim et al., Cell, 2019 - slow-oscillations and delta-waves detection code
clear; close all;

load('example_data.mat'); %load example data
data

% detection of slow-oscillations and delta-waves
% you need to input -   LFP and sampling frequency
%                       logicals or binary index of NREM sleep detection (you would have your detections for your data)
%                       logicals or binary index of artifact detection. you may set it zeros to ignore it
%                       plotting results 1|0
%                       parameters for [peak-thr trough-thr dur-min dur-max]
so_delta=detect_so_delta(data.LFP,data.Fs_LFP,...
            'sleep_idx',data.sleep_idx,...
            'artifact_idx',data.artifact_idx,...
            'PLOT',1,...
            'mnl_parm',[85 40 .15 .5]); 