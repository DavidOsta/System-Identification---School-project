%% Script pro testovani vytvorene trydi apod.

clear;clc;close all;
checkcode test_script.m


%% Loading data
addpath(genpath('datasets'));
% Impulses

cart_impulses_A_loaded = load('A_last_week/pulses_struct.mat');
measured_impulses_A = cart_impulses_A_loaded.pulses_cell;

cart_impulses_B_loaded = load('B_last_week/pulses_struct.mat');
measured_impulses_B = cart_impulses_B_loaded.pulses_cell;

% Sinus

cart_sinus_A_loaded = load('A_last_week/sinus_struct.mat');
measured_sinus_A = cart_sinus_A_loaded.sinuses_cell;

cart_sinus_B_loaded = load('B_last_week/sinus_struct.mat');
measured_sinus_B = cart_sinus_B_loaded.sinuses_cell;


%% Processing

load_system('sim_friction_model');
estimators_struct = init_parametrs();

% imp_len = length(measured_impulses_A);
% sin_len = length(measured_sinus_A);
% 
% sim_results_pulses = cell(imp_len ,2);
% sim_results_sinuses = cell(sin_len ,2);
% 
% sim_errors_pulses = zeros(imp_len,2);
% sim_errors_sinuses = zeros(sin_len,2);
% 
% 
% sample = cart_impulses_A_loaded.pulses_cell{5};
% input_signal = sample.force_input;
% measured_position = sample.measured_position;
% sim('sim_friction_model');
%-------------------------------------------------------------------------
% GA with 3 generations, 10 trials 2 best keep
tic;
fm1 = Friction_model(measured_impulses_A,20);
fm1 = fm1.find_friction_model();
fm1.print_sim_results();

toc;
