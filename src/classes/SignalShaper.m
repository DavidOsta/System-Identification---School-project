classdef SignalShaper
    %SIGNALSHAPER Class for obtaining signal shapers models based on given
    %parameters
    %   1st arg - natural frequency (rad/sec), type - double
    %   2nd arg - damping ratio (-), type - double

    properties(Access = private)
        natural_freq;
        damping_ratio;
        complex_conjugate_poles;
        beta_re;
        damped_natural_freq;

    end

    methods(Access = public)

        function self = SignalShaper(varargin)

            if isnumeric(varargin)
                error('wca:WrongConstructorsArg',...
                    ['Wrong constructors argument']);
            elseif nargin == 2
                self.natural_freq = varargin{1};
                self.damping_ratio = varargin{2};
                self = self.init_parameters();
                
            elseif nargin == 0
                error('cam:ConstructorArgMiss',...
                    'constructor argument is missing');
            else
                error('wca:WrongConstructorsArg',...
                    ['Wrong constructors argument']);
            end

        end

          function [Omega, beta, complex_poles] = get_init_parameters(self)
            Omega = self.damped_natural_freq;
            beta = self.beta_re;
            complex_poles = self.complex_conjugate_poles;
          end

        function [ZV_shaper, shaper_par] = get_ZV_shaper(self)
            %GET_ZV_SHAPER returns tf of ZV shaper and its parameters
            [a, tau_add, tau_sub] = self.get_ZV_shapers_parameters(1);
            s = tf('s');

            ZV_shaper = a + (1-a) * exp(-tau_sub*s);

            shaper_par.gain_a = a;
            shaper_par.tau_add = tau_add;
            shaper_par.tau_sub = tau_sub;
            shaper_par.alpha = 0;
        end

        function [ZV_inv_shaper, shaper_par] = get_inverse_ZV_Shaper(self)
            %GET_INVERSE_ZV_SHAPER returns tf of inverse ZV shaper and its parameters 
            [ZV_shaper, shaper_par] = self.get_ZV_shaper();
            ZV_inv_shaper = 1 / ZV_shaper;
        end


        function [DZV_alpha_shaper, shaper_par] = get_DZV_alpha_Shaper(self,alpha)
            %GET_DZV_ALPHA_SHAPER returns tf of DZV_alpha shaper and its parameters
            [a, tau_add, tau_sub] = self.get_shapers_basic_parameters(alpha);
%             [a_gain1, tau_add1, tau_sub1] = self.get_shapers_optimal_parameters(alpha)
            s = tf('s');
%
%             DZV_alpha_shaper =...
%                 a + ((1-a) / (1-alpha))*( (exp(-s*alpha*tau)-exp(-s*tau)) / (s*tau) );
            DZV_alpha_shaper =...
                a + ((1-a) / (1-alpha))*( (exp(-s*tau_add)-exp(-s*tau_sub)) / (s*tau_sub) );


            shaper_par.gain_a = a;
            shaper_par.tau_add = tau_add;
            shaper_par.tau_sub = tau_sub;
            shaper_par.alpha = alpha;
        end

        function [DZV_alpha_inv_shaper, shaper_par] = get_inverse_DZV_alpha_Shaper(self,alpha)
            %GET_INVERSE_DZV_ALPHA_SHAPER returns tf of DZV_alpha shaper and its parameters

            [DZV_alpha_shaper, shaper_par] = self.get_DZV_alpha_shaper(alpha);
            DZV_alpha_inv_shaper = 1 / DZV_alpha_shaper;
        end


        function [DZV_shaper, shaper_par] = get_DZV_shaper(self)
            %GET_DZV_SHAPER returns ss model of DZV shaper and its parameters
            s = tf('s');

            [a, tau_add, tau_sub] = self.get_DZV_shapers_parameters();
            
             DZV_shaper =...
                a + ((1-a)/tau_sub)*((1-exp(-s*tau_sub)) / s);
            
            shaper_par.gain_a = a;
            shaper_par.tau_add = tau_add;
            shaper_par.tau_sub = tau_sub;
            shaper_par.alpha = 0;
        end

        function [DZV_inv_shaper, shaper_par] = get_inverse_DZV_shaper(self)
            %GET_INVERSE_DZV_SHAPER returns ss model of inverse DZV shaper
            %and its parameters
            [DZV_shaper, shaper_par] = self.get_DZV_shaper();
            DZV_inv_shaper = 1 / DZV_shaper;
        end

%         function [DZV_inv_opt_shaper] = get_inverse_DZV_optimal_shaper(self,alpha)
%         [a, tau_add, tau_sub] = get_shapers_optimal_parameters(self,alpha);
% 
%         s = tf('s');
%         DZV_opt_shaper =...
%             a + ((1-a) / tau_sub)*( (exp(-s*tau_add)-exp(-s*(tau_sub+tau_add))) / s );
%         DZV_inv_opt_shaper = 1/DZV_opt_shaper;
%         end





    end

    methods(Access = private)

        function self =  init_parameters(self)
            self.beta_re = self.natural_freq * self.damping_ratio;
            self.damped_natural_freq =...
                self.natural_freq * sqrt(1 - self.damping_ratio^2);
            pole1 = complex(-self.beta_re, (self.damped_natural_freq));
            pole2 = complex(-self.beta_re, -(self.damped_natural_freq));
            self.complex_conjugate_poles = [pole1; pole2];

        end

        
        function [a_gain, tau_add, tau_sub] = get_ZV_shapers_parameters(self, alpha)
            % GET_SHAPERS_BASIC_PARAMETERS calculate shapers parameter
            exp_term = exp(self.beta_re*pi / self.damped_natural_freq);
            a_gain = exp_term / (1 + exp_term);
            tau_sub = pi / self.damped_natural_freq;
            tau_add = tau_sub * alpha;
        end
        
        
        function [a_gain, tau_add, tau_sub] = get_DZV_shapers_parameters(self)
            % GET_SHAPERS_BASIC_PARAMETERS calculate shapers parameter
            
            function r = find_root(tau)
                %helper function to find root of non-lin. eq.
                r = wo*exp(-b*tau) + b*sin(wo*tau) - wo*cos(wo*tau);
            end
            
            b = self.beta_re;
            wo = self.damped_natural_freq;
            
            tau0 = pi/wo; % initial guess, tau lies between <pi/wo;2pi/wo>
            l_bound = tau0;
            u_bound = 2*pi/wo;
            tau = fzero(@find_root, tau0);
            
            if(tau > u_bound || tau < l_bound)
                error('found root doesn'' lie in interval <pi/wo;2*pi/wo ');
            end
            
            a_gain = sin(wo * tau) / (sin(wo * tau) - tau * wo * exp(-b*tau));
            
            tau_sub = tau;
            tau_add = 0;
            
        end
        
        
        function [a_gain, tau_add, tau_sub] = get_shapers_parameters(self, alpha)
            % GET_SHAPERS_BASIC_PARAMETERS calculate shapers parameter
            exp_term = exp(self.beta_re*pi / self.damped_natural_freq);
            a_gain = exp_term / (1 + exp_term);
            tau_sub = pi / self.damped_natural_freq;
            tau_add = tau_sub * alpha;
        end
        
        
        function [a_gain, tau_add, tau_sub] = get_shapers_optimal_parameters(self,alpha)
            
            W = self.damped_natural_freq;
            B = self.beta_re;
            tau = pi / W;
            tau_sub = alpha * tau;
            G = (1 - exp(-(B+1j*W)*tau_sub))/(tau_sub*(B+1j*W));
            m = abs(G);
            phi = angle(G);
            
            exp_term = m * exp(-B/W * (phi + pi));
            
            a_gain = exp_term / (1 + exp_term);
            tau_add = (pi + phi)/W;
            
        end


    end

end
