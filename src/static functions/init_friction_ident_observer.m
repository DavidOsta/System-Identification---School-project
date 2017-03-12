function [A,B,C,D,L,X,Xhat] =  init_friction_ident_observer()

%% State-Space System

f_1k_lin = 1.3;
m1 = 0.5; % orig 1.81, ale dynamika byla moc pomala

A = [0, 1;
     0, -f_1k_lin/m1];
 
B = [0; 1/m1];

C = [1 0;
     0 1];
D = [0;0];
X = [0, 0]; % init conditions

%% State-Space Observer
% cim vetsi hodnota -> rychlejsi dynamika ->
% vetsi problem s sumem z mereni
l1 = 10;
l2 = 10; 
L = [l1; l2];
   
Xhat = [0, 0]; % estimator initial conditions




% As the gain of the estimator (KF or LO)  increases,
% the noise attenuation decreases but the 'speed' of
% the transient response increases. Conversely, as the gain decreases,
% noise attenuation (i.e. smoothing) increases but the transient response becomes
% more 'sluggish'. 

% To get a high-gain Kalman Filter (KF) use a large Q and a small R.
% For a high-gain Luenberg Observer (LO) place the poles of
% the observer further to the LHS of the s plane (for the continuous-time case)
% ,or closer to the origin of the z plane (for the discrete-time case).   

% The main difference between the KF and the LO is that the LO gain
% is constant and independent of time; for the KF,
% the gain is time dependent - it decreases over time until
% steady-state is reached. As a result,
% I would expect the LO to be less computationally expensive
% (because a variance update is not required).
end


% % Navrh polu jenom pro vozik
% A_cart = A(1:2, 1:2);
% B_cart = B(1:2);
% % nejake zvolene poly
% P_cart = [-20, -0.5];
% % ziskane koeficienty
% K_cart = place(A_cart,B_cart,P_cart);
% 
% % Nova matice dynamiky
% A_cart_new = A_cart - B_cart * K_cart;
% 
% 
% [poles_test_R, poles_test_IM] = get_poles(A_cart);
% [poles_test_R2, poles_test_IM2] = get_poles(A_cart_new);
% 
% close all
% figure;grid on;hold on;
% plot(poles_test_R, poles_test_IM , 'rx' ,'linewidth',2, 'markersize', 10);
% plot(poles_test_R2, poles_test_IM2 , 'gx' ,'linewidth',2, 'markersize', 10);
% 
% 
% % stavovy regulator
% % k1 = 0.5;
% % k2 = 100; % 0.001
% 
% k1 = K_cart(1);
% k2 = K_cart(2);
% k3 = 0;
% k4 = 0;
% 
% % nova matice dynamiky celeho systemu
% A_new = A - B * [k1 k2 k3 k4];
% 
% % Poles of original system and system with new matrix of dynamics
% [poles_R, poles_IM] = get_poles(A);
% [poles_new_R, poles_new_IM, egVals] = get_poles(A_new);
% figure;grid on; hold on;
% p1 = plot(poles_R,poles_IM, 'rx' ,'linewidth',2, 'markersize', 10);hold on
% p2 = plot(poles_new_R, poles_new_IM, 'gx' ,'linewidth',2, 'markersize', 10);
% 
% 
% % Zeros of original system and system with new matrix of dynamics
% 
% Z = get_zeros(A,B,C,D);
% Z_new = get_zeros(A_new,B,C,D);
% 
% 
% z1 = plot(Z, 'bo', 'linewidth',2, 'markersize', 10);
% z2 = plot(Z_new, 'mo', 'linewidth',2, 'markersize', 10);
% title('Poles and Zeros of linearized system');
% legend('Poles','new Poles','Zeros','new Zeros');
% xlabel('$Re(s)$', 'interpreter', 'latex');
% ylabel('$Im(s)$', 'interpreter', 'latex');
% 
% % Controlability
% 
% 
% sys_ss = ss(A,B,C,D);
% controll_matrix = ctrb(sys_ss);
% % controll_matrix is 4x4 => controllability must have rank(4)
% controllability = rank(controll_matrix);
% 
% % LQR design
% Q = C'*C;
% R = 1;
% K = lqr(A,B,Q,R);
% 
% % observalibility 
% 
% ob = obsv(sys_ss);
% % controll_matrix is 4x4 => controllability must have rank(4)
% %  the observability matrix is 8x4 and has rank 4, it is full rank and our system is observable.
% observability = rank(ob);
% 
% poles_A_new = eig(A_new);
% P_obs = [-5 -6 -7 -8];
% L = place(A',C',P_obs)'
% 
% 
% Ace = [(A-B*K) (B*K);
%        zeros(size(A)) (A-L*C)];
% Bce = [B;zeros(size(B))];
% Cce = [C zeros(size(C))];
% Dce = [0;0];
% 
% states = {'x' 'x_dot' 'phi' 'phi_dot' 'e1' 'e2' 'e3' 'e4'};
% inputs = {'r'};
% outputs = {'x'; 'phi'};
% 
% sys_est_cl = ss(Ace,Bce,Cce,Dce,'statename',states,'inputname',inputs,'outputname',outputs);
% 
% figure
% t = 0:0.01:50;
% r = 0.2*ones(size(t));
% [y,t,x]=lsim(sys_est_cl,r,t);
% [AX,H1,H2] = plotyy(t,y(:,1),t,y(:,2),'plot');
% set(get(AX(1),'Ylabel'),'String','cart position (m)')
% set(get(AX(2),'Ylabel'),'String','pendulum angle (radians)')
% title('Step Response with Observer-Based State-Feedback Control')
