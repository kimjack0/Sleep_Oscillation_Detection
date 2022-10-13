% Kim et al., 2022 - ripple detection code
clear; close all;

load('example_data.mat'); %load example data
data

% detection of ripples
% you need to input -   LFP and sampling frequency
%                       detection prameters: low threshold / high threshold / min duration / max duration
%                       logicals or binary index of NREM sleep detection (you would have your detections for your data)
%                       logicals or binary index of artifact detection. you may set it zeros to ignore it
%                       plotting results
%                       using only sleep period
%                       ripple frequency band

% outpus -   pks - time peak of filtered lfp signal for each thresh crossing
%            start - upward thresh crossing
%            finish - downward thresh crossing
%            amp - peak amplitude of envelope
%            dur - duration (secs)


%%
fpass = [150,250];
session_size=[size(data.LFP,1)];
th_dur=[1.0 4.0 0.03 100]; % low threshold / high threshold / min duration / max duration
ripples = detect_ripples(mat2cell(data.LFP, session_size, [1]),...
    th_dur,...
    'Fs',data.Fs_LFP,...
    'sleep_idx',mat2cell(data.sleep_idx, session_size, [1]),...
    'artifact_idx',mat2cell(data.artifact_idx, session_size, [1]),...
    'PLOT',1,...
    'sleep_classify',1,...
    'fpass',fpass);