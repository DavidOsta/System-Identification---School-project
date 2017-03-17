classdef LinearSystem
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        system_type;
        
        m_a;    % mass of cart A
        v_f_a;  % visc friction of cart A
        m_b;    % mass of cart B
        v_f_b;  % visc friction of cart B
        m_p;    % mass of pendulum
        l_p;    % length of pendulum
        b_p;    % damping of pendulum
        b_s;    % damping of spring
        k_s;    % stiffnes of spring
        
        system_state_space_struct;
        
        observer_state_space_struct;
        
        g = 9.81;
        model_name;
        parameters_folder_name;
        system_parameters;
        system_block_name = 'state_space_system';
        observer_block_name = 'state_space_observer';
    end
    
    methods(Access = public)
        % Constructor
        function self = LinearSystem(varargin)
            
            if nargin > 0
                self.model_name = varargin{1};
                self.system_type = varargin{2};
                self.parameters_folder_name = varargin{3};
                self = self.load_parameters();
                self = self.init_parameters();
            else
                error('constructor argument is missing')
            end
            
        end
        
        function ss_par = get_state_space_parameters(self)
            ss_par = self.system_state_space_struct;
        end
        function ss_par = get_observers_parameters(self)
            ss_par = self.observer_state_space_struct;
        end
        
        
    end
    
    methods(Access = private)
        
        function self = init_parameters(self)
            % INIT_PARAMETERS initialize parameters of linear systems
            % [OUTPUTS] = INIT_PARAMETERS(ARGS)
            %
            % :param args: -
            % :type args: -
            % :returns: struct of parameters
            % :raises: invalid type of system
            
            switch self.system_type
                case 'cart'
                    [state_space_struct, observer_struct] =...
                        self.get_cart_state_space_models();
                case 'two_cart_pendulum'
                    [state_space_struct, observer_struct] =...
                         self.get_two_cart_pendulum_state_space_models();
                case 'cart_pendulum'
                    [state_space_struct, observer_struct] =...
                         self.get_cart_pendulum_state_space_models();
                otherwise
                    error('unknow type of system');
            end

            self.system_state_space_struct = state_space_struct;
            setup_statespace_block(self.model_name,...
                self.system_block_name,...
                state_space_struct);
            
            self.observer_state_space_struct = observer_struct;
            setup_statespace_block(self.model_name,...
                self.observer_block_name,...
                observer_struct);
        end
        
        function [state_space_struct, observer_struct] =...
                get_cart_state_space_models(self)
            
            A = [0, 1; 0, -self.v_f_a/self.m_a];
            B = [0; 1/self.m_a];
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
            
            state_space_struct = self.create_state_space_struct(A,B,C,D,X);
            observer_struct = self.create_observer_struct(Ahat,Bhat,Chat,Dhat,Xhat,L);
        end
        
        function [state_space_struct, observer_struct] =...
                get_two_cart_pendulum_state_space_models(self)
            
            A = [0, 1, 0, 0, 0, 0;
                -self.k_s/self.m_a, -(self.b_s + self.v_f_a)/self.m_a, self.k_s/self.m_a, self.b_s/self.m_a, 0, 0;
                0, 0, 0, 1, 0, 0;
                self.k_s/self.m_b, self.b_s/self.m_b, -self.k_s/self.m_b, -(self.b_s + self.v_f_b)/self.m_b, self.m_p*self.g/self.m_b, self.b_p/(self.l_p*self.m_b);
                0, 0, 0, 0, 0, 1;
                -self.k_s/(self.l_p*self.m_b), -self.b_s/(self.l_p*self.m_b), self.k_s/(self.l_p*self.m_b), (self.b_s + self.v_f_b)/(self.l_p*self.m_b), -(self.m_b + self.m_p)*self.g/(self.l_p*self.m_b), -self.b_p*(1+self.m_b/self.m_p)/(self.m_b*self.l_p^2)
                ];
            
            B = [0; 1/self.m_a; 0; 0; 0; 0];
            C = [1 0 0 0 0 0;
                0 0 0 0 0 0; % 0 1 0 .... - for velocity
                0 0 1 0 0 0;
                0 0 0 0 0 0;
                0 0 0 0 1 0;
                0 0 0 0 0 0];
            D = [0; 0; 0; 0; 0; 0];
            X = [0, 0, 0, 0, 0, 0];
            
            
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
            Dhat = [0, 0; 0, 0; 0, 0; 0, 0; 0, 0; 0, 0];
            Xhat = [0.1, 0, 0, 0, 0, 0];
            
            state_space_struct = self.create_state_space_struct(A,B,C,D,X);
            observer_struct = self.create_observer_struct(Ahat,Bhat,Chat,Dhat,Xhat,L);
        end
        
        function [state_space_struct, observer_struct] = get_cart_pendulum_state_space_models(self)
            A = [0, 1, 0, 0;
                0, -self.v_f_a/self.m_a, (self.m_p*self.g)/self.m_a, self.b_p/(self.m_a*self.l_p);
                0, 0, 0, 1;
                0, self.v_f_a/(self.l_p*self.m_a), -((self.m_a+self.m_p)*self.g)/(self.l_p*self.m_a), -(self.b_p*(1+self.m_a/self.m_b))/(self.m_a*self.l_p^2)];
            
            B = [0; 1/self.m_a; 0; -1/(self.l_p*self.m_a)]; % 1/m1, jsem mel na prvni pozici misto druhe
            
            C = [1 0 0 0;
                0 0 1 0];
            D = [0; 0];
            X = [0, 0, 0, 0];

            state_space_struct = self.create_state_space_struct(A,B,C,D,X);

            l1 = 10;
            l2 = 10;
            l3 = 10;
            l4 = 10;
            L = [l1; l2; l3; l4];
            Ahat = A;
            Bhat = [B, L];
            Chat = C;
            Dhat = [0, 0; 0, 0];
            Xhat = [0.1, 0, 0, 0];
            
            observer_struct = self.create_observer_struct(Ahat,Bhat,Chat,Dhat,Xhat,L);
            
        end
        
        
        
        function ss_struct = create_state_space_struct(self,A,B,C,D,X)
            ss_struct.A = A;
            ss_struct.B = B;
            ss_struct.C = C;
            ss_struct.D = D;
            ss_struct.X = X;
        end
        
        function os_struct = create_observer_struct(self,A,B,C,D,X,L)
            os_struct = self.create_state_space_struct(A,B,C,D,X);
            os_struct.L = L;
        end
        
        
        function self = load_parameters(self)
            % load stored parameters
            load_data = load(self.parameters_folder_name);
            fnames = fieldnames(load_data);
            
            measured_parameters = load_data.(fnames{1});
            
            self.m_a = measured_parameters.mass_cart_A;
            self.v_f_a = measured_parameters.visc_friction_cart_A;
            self.m_b = measured_parameters.mass_cart_B;
            self.v_f_b = measured_parameters.visc_friction_cart_B;
            self.m_p = measured_parameters.mass_pendulum;
            self.l_p = measured_parameters.length_pendulum;
            self.b_p = measured_parameters.damping_pendulum;
            self.b_s = measured_parameters.damping_spring;
            self.k_s = measured_parameters.stiffnes_spring;
        end
        
    end
end

