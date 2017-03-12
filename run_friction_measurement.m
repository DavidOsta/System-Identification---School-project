clear;clc;close all;
checkcode run_friction_measurement.m
 
addpath(genpath('models'));
addpath(genpath('data'));
addpath(genpath('src'));



fm = FrictionMeasurement();
fm = fm.measure_response();
fm.plot_results();
fm.save_data();

