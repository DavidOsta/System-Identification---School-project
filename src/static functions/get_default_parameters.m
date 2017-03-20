function [ parameters ] = get_default_parameters()
%GET_DEFAULT_PARAMETERS Summary of this function goes here
%   Detailed explanation goes here

parameters.mass_cart_A = 2.7000;
parameters.visc_friction_cart_A = 7.0200;
parameters.mass_cart_B =1.5500;
parameters.visc_friction_cart_B = 0.5000;
parameters.mass_pendulum = 0.055;
parameters.damping_pendulum =  0.0016;
parameters.length_pendulum = 0.4667;
parameters.damping_spring =5.8200;
parameters.stiffnes_spring = 584.9100;
parameters.g = 9.81;
parameters.nat_freq_pendulum =  4.6757;
parameters.damping_rat_pendulum = 2.0400e-04;


end

