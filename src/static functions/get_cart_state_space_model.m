function [ state_space_struct, observer_struct ] =...
    get_cart_state_space_model( parameters_file_path )
%GET_CART_STATE_SPACE_MODEL Summary of this function goes here
%   Detailed explanation goes here

parameters = load_parameters( parameters_file_path );

m_a = parameters.mass_cart_A;
v_f_a = parameters.visc_friction_cart_A;


A = [0, 1; 0, -v_f_a/m_a];
B = [0; 1/m_a];
C = [1 0; 0 1];
D = [0; 0];
X = [0, 0];

l1 = 10;
l2 = 10;
L = [l1; l2];
Ahat = A;
Bhat = [B, L];
Chat = [1 0; 0 1];
Dhat = [0, 0; 0, 0];
Xhat = [0.1, 0];

state_space_struct.A = A;
state_space_struct.B = B;
state_space_struct.C = C;
state_space_struct.D = D;
state_space_struct.X = X;

observer_struct.A = Ahat;
observer_struct.B = Bhat;
observer_struct.C = Chat;
observer_struct.D = Dhat;
observer_struct.X = Xhat;
observer_struct.L = L;

end

