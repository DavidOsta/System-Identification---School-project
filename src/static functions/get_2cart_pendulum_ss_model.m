function [ system_ss, observer_ss] =...
    get_2cart_pendulum_ss_model( parameters_file_path )
%GET_CART_STATE_SPACE_MODEL Summary of this function goes here
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
m_b = parameters.mass_cart_B;
v_f_b = parameters.visc_friction_cart_B;
m_p = parameters.mass_pendulum;
b_p = parameters.damping_pendulum;
l_p = parameters.length_pendulum;
b_s = parameters.damping_spring;
k_s = parameters.stiffnes_spring;

A = [0, 1, 0, 0, 0, 0;
    -k_s/m_a, -(b_s + v_f_a)/m_a, k_s/m_a, b_s/m_a, 0, 0;
    0, 0, 0, 1, 0, 0;
    k_s/m_b, b_s/m_b, -k_s/m_b, -(b_s + v_f_b)/m_b, m_p*g/m_b, b_p/(l_p*m_b);
    0, 0, 0, 0, 0, 1;
    -k_s/(l_p*m_b), -b_s/(l_p*m_b), k_s/(l_p*m_b), (b_s + v_f_b)/(l_p*m_b), -(m_b + m_p)*g/(l_p*m_b), -b_p*(1+m_b/m_p)/(m_b*l_p^2)
    ];
B = [0; 1/m_a; 0; 0; 0; 0];
C = [1 0 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 0 1 0];

D = [0; 0; 0];

% for all states as outputs
% C = [1 0 0 0 0 0;
%     0 1 0 0 0 0; % 0 1 0 .... - for velocity
%     0 0 1 0 0 0;
%     0 0 0 1 0 0;
%     0 0 0 0 1 0;
%     0 0 0 0 0 1];
% D = [0; 0; 0; 0; 0; 0];

system_states = {'xa' 'xa_dot' 'xb' 'xb_dot' 'phi' 'phi_dot'};
system_inputs = {'u'};
system_outputs = {'xa'; 'xb' ; 'phi'};

system_ss = ss(A,B,C,D,'statename',system_states,...
    'inputname', system_inputs,...
    'outputname',system_outputs);

l1 = 10;
l2 = 10;
l3 = 10;
l4 = 10;
l5 = 10;
l6 = 10;

L = [l1; l2; l3; l4; l5; l6];
Ahat = A;
Bhat = [B, L];
Chat = C;
Dhat = [D,D];

observer_states = {'xa' 'xa_dot' 'xb' 'xb_dot' 'phi' 'phi_dot'};
observer_inputs = {'u', 'e'};
observer_outputs = {'xa'; 'xb' ; 'phi'};

observer_ss = ss(Ahat,Bhat,Chat,Dhat,'statename',observer_states,...
    'inputname', observer_inputs,...
    'outputname',observer_outputs);

end

