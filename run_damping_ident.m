% clear;clc; close all
% checkcode run_damping_ident.m
add_paths();

data =...
    load('damping_measurement-15-Mar-2017 21:21:43');

% rename to cart
cart_position = data.measured_data.position;
pendulum_angle = data.measured_data.angle;
impulse_signal = data.measured_data.impulse_signal;

pendulum_weight = 0.055 ; %[kg]
pendulum = DampingModel(pendulum_angle, impulse_signal, pendulum_weight);
pendulum_parameters = pendulum.get_parameters;
pendulum.plot_comparison

cart_weight = 1.55; %[kg]
cart = DampingModel(cart_position, impulse_signal, cart_weight);
cart_parameters = cart.get_parameters;
cart.plot_comparison


