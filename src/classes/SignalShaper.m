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
            s = tf('s');

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


        function [DZV_shaper, shaper_par] = get_DZV_shaper_slow(self)
            %GET_DZV_SHAPER returns ss model of DZV shaper and its parameters
            s = tf('s');

            [a, tau_add, tau_sub] = self.get_DZV_shapers_parameters_analytical();
            
             DZV_shaper =...
                a + ((1-a)/tau_sub)*((1-exp(-s*tau_sub)) / s);
            
            shaper_par.gain_a = a;
            shaper_par.tau_add = tau_add;
            shaper_par.tau_sub = tau_sub;
            shaper_par.alpha = 0;
        end
        function [DZV_inv_shaper, shaper_par] = get_inverse_DZV_shaper_slow(self)
            %GET_INVERSE_DZV_SHAPER returns ss model of inverse DZV shaper
            %and its parameters
            [DZV_shaper, shaper_par] = self.get_DZV_shaper_slow();
            DZV_inv_shaper = 1 / DZV_shaper;
        end
        
        
        
        function [DZV_shaper, shaper_par] = get_DZV_shaper(self)
            %GET_DZV_SHAPER returns ss model of DZV shaper and its parameters
            s = tf('s');
            
            [a, tau_add, tau_sub] = self.get_DZV_shapers_parameters_analytical();
            
                    [a_gain, tau_add1, tau_sub1] = self.get_DZV_shapers_parameters_numerical()
    
            
            
            DZV_shaper =...
                a + ((1-a)/tau_sub)*((1-exp(-s*tau_sub)) / s) * exp(-s*tau_add);
            
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
        
        
        function [a_gain, tau_add, tau_sub] =...
                get_DZV_shapers_parameters_numerical(self)
            % GET_SHAPERS_BASIC_PARAMETERS calculate shapers parameter
            
            function r = find_root(tau)
                %helper function to find root of non-lin. eq.
                r = w0*exp(-b*tau) + b*sin(w0*tau) - w0*cos(w0*tau);
            end
            
            b = self.beta_re;
            w0 = self.damped_natural_freq;
            
            tau0 = pi/w0; % initial guess, tau lies between <pi/wo;2pi/wo>
            l_bound = tau0;
            u_bound = 2*pi/w0;
            tau = fzero(@find_root, tau0);
            
            if(tau > u_bound || tau < l_bound)
                error(['found root doesn'' lie in interval <pi/wo;2*pi/wo\n',...
                ' numerical method has probably failed ']);
            end
            
            a_gain = sin(w0 * tau) / (sin(w0 * tau) - tau * w0 * exp(-b*tau));
            
            tau_sub = tau;
            tau_add = 0;
           
         
        end
        
        function [a_gain, tau_add, tau_sub] =...
                get_DZV_shapers_parameters_analytical(self)
            % GET_SHAPERS_BASIC_PARAMETERS calculate shapers parameter
            
           
            b = self.beta_re;
            w0 = self.damped_natural_freq;
            tau_sub = pi / w0;

            
            G = (1 - exp(-tau_sub * (b + 1j * w0))) / (tau_sub * (b + 1j * w0));
            
            m = abs(G);
            phi = angle(G);
            exp_term = m * exp((b / w0) * (pi + phi));
            a_gain = exp_term / (1 + exp_term);
            tau_add = (pi + phi) / w0;
        end
        
        
        function [a_gain, tau_add, tau_sub] = get_shapers_parameters(self, alpha)
            % GET_SHAPERS_BASIC_PARAMETERS calculate shapers parameter
            exp_term = exp(self.beta_re*pi / self.damped_natural_freq);
            a_gain = exp_term / (1 + exp_term);
            tau_sub = pi / self.damped_natural_freq;
            tau_add = tau_sub * alpha;
        end
        
        
   end

end
