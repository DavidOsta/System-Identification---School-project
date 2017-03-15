clear;clc; close all
add_paths();

lol = DampingMeasurement;
lol = lol.measure_response;
saved_data = lol.get_results_path;