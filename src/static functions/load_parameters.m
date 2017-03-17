function [ loaded_parameters ] = load_parameters( file_path )
%LOAD_PARAMETERS Summary of this function goes here
%   Detailed explanation goes here

% load stored parameters
load_data = load(file_path);
fnames = fieldnames(load_data);

loaded_parameters = load_data.(fnames{1});

end

