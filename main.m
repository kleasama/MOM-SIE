clear;clc;close all
%% Process user input
UserInput;
%% Calculate EM parameters
EMparams;
%% Load model and create connectivity
load(['models\',model]);
if show_mesh == 1;viewer(p,t);end
rwg = RWG(t);
tcn = TCN(rwg, length(t(:,1)));