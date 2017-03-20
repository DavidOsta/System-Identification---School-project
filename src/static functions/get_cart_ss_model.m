function [ system_ss, observer_ss ] =...
    get_cart_ss_model( parameters_file_path )
%GET_CART_SS_MODEL_FOR_PARFOR State space model of cart with friction
%   Detailed explanation goes here

if nargin == 1
    parameters = load_parameters( parameters_file_path );
    g = 9.81;
elseif nargin == 0
    parameters  = get_default_parameters();
    g = parameters.g;
else
    error('Invalid parameters');
end
m_a = parameters.mass_cart_A;
v_f_a = parameters.visc_friction_cart_A;

A = [0, 1; 0, -v_f_a/m_a];
B = [0; 1/m_a];
C = [1 0; 0 1];
D = [0; 0];
X = [0;0];

l1 = 10;
l2 = 10;
L = [l1; l2];
Ahat = A;
Bhat = [B, L];
Chat = [1 0; 0 1];
Dhat = [0, 0; 0, 0];
Xhat = X;

system_ss.A = A;
system_ss.B = B;
system_ss.C = C;
system_ss.D = D;
system_ss.X = X;

observer_ss.A = Ahat;
observer_ss.B = Bhat;
observer_ss.C = Chat;
observer_ss.D = Dhat;
observer_ss.X = Xhat;

end

