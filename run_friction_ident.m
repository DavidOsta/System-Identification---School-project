clear;clc;close all;
addpath(genpath('src'));

add_paths();


% This is temporary
 cart_impulses_A_loaded = load('A_last_week/pulses_struct.mat');
measured_impulses_A = cart_impulses_A_loaded.pulses_cell;


tic;

% there are some issues with parallel parloop and simulink's slprj files on Ubuntu
% 16.4, need to fix it later
try 
fm1 = FrictionModel(measured_impulses_A,10);
fm1 = fm1.find_friction_model();
fm1.print_sim_results();
catch ME
    switch ME.identifier
        case {'MATLAB:load:cantReadFile', 'MATLAB:nonExistentField',...
              'MATLAB:load:notBinaryFile', 'MATLAB:load:unableToReadMatFile'}
            rmdir('slprj','s');
            run_friction_ident;
        case 'Simulink:Commands:OpenSystemUnknownSystem'
            load_system('sim_friction_model');
            run_friction_ident;
        otherwise
            rethrow(ME)
    end
end

toc;
