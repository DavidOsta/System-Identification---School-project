clear;clc; close all
add_paths();

data =...
    load('damping_measurement-13-Mar-2017 19:56:16');

% rename to cart
cart_position = data.measured_data.position;
pendulum_angle = data.measured_data.angle;

tic;
damp_cart = DampingIdent(cart_position);

init_vals = [-500; 0.1; 3]; % gain, damp_ratio, nat_freq
l_bounds = [-1000; 0; 0];
u_bounds = [1000; 2; 100];

damp_cart =...
    damp_cart.estimate_parameters(init_vals,l_bounds,u_bounds);
damp_cart.get_best_pop()

% damp_pendulum = DampingIdent(pendulum_angle);
% 
% init_vals = [800; 0.115; 20]; % gain, damp_ratio, nat_freq
% l_bounds = [790; 0.05; 17];
% u_bounds = [810; 0.15; 21];
% 
% damp_pendulum =...
%     damp_pendulum.estimate_parameters(init_vals,l_bounds,u_bounds);
% damp_pendulum.get_best_pop()

toc;