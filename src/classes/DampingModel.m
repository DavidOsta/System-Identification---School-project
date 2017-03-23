classdef DampingModel
    %DAMPINGIDENT Summary of this class goes here
    %   Detailed explanation goes here

    properties(Access = private)
        input_signal;
        output_signal;
        estimated_parameters;
        time_sample;
        estimated_sys;
        parameters;
        mass;
    end

    methods(Access = public)
        function self = DampingModel(varargin)
            if nargin == 3
                self.output_signal = varargin{1};
                self.input_signal = varargin{2};
                self.mass = varargin{3};
                self.time_sample =...
                    self.input_signal.Time(2) - self.input_signal.Time(1);
                self.estimated_sys = self.find_sys;
            elseif nargin == 0
                error('cam:ConstructorArgMiss',...
                    'constructor argument is missing');
            else
                error('wca:WrongConstructorsArg',...
                    ['Wrong constructors argument']);
            end

        end

        function sys_par = get_parameters(self)
            sys_koef = getpvec(self.estimated_sys);
            b0 = sys_koef(1);
            a0 = sys_koef(2);
            a1 = sys_koef(3);
            a2 = sys_koef(4);

            gain = b0;
            natural_freq = sqrt(a1);
            damping_ratio = a0 / (natural_freq * 2);

            k_stiff =  self.mass * natural_freq^2; % [N/m] - Stiffness
            c_crit_damp = 2 * self.mass * natural_freq; % [N.s / m] - Critical damping
            c_damp = damping_ratio * c_crit_damp; % [N.s / m] - damping, z = c / c_c

            sys_par = struct(...
                'gain',gain,...
                'damping_ratio',damping_ratio,...
                'natural_freq',natural_freq,...
                'stiffnes', k_stiff,...
                'damping', c_damp);
        end

        function self = plot_comparison(self)
            % plot compare of measured data and found model
            figure;
            [y_sys,t_sys] =...
                impulse(self.estimated_sys, self.output_signal.Time);

            plot(self.output_signal);grid on; hold on;
            plot(t_sys,y_sys);

            sys_par = self.get_parameters;
            title('Impulse response');
            xlabel('Time seconds');
            ylabel('Position');

            dim = [.2 .5 .4 .4];
            str = strcat('k = ', num2str(sys_par.stiffnes),...
                ' [N/m]', ' c = ', num2str(sys_par.damping), ' [N.s/m]');
            annotation('textbox',dim,'String',str,'FitBoxToText','on');
            legend('Measured Output', 'Simulated Output');



            time_axes = self.output_signal.Time;
            [y_min, indx_min] = min(y_sys);
            envelope_min =...
                y_min * exp(-sys_par.damping_ratio * sys_par.natural_freq*time_axes); %* sqrt(1/(1-sys_par.damping_ratio^2)


            [y_max, indx_max] = max(y_sys);
            envelope_max =...
                y_max * exp(-sys_par.damping_ratio * sys_par.natural_freq*time_axes);

            plot(time_axes + time_axes(indx_min), envelope_min, 'k--');
            plot(time_axes + time_axes(indx_max), envelope_max, 'k--');
        end
    end

    methods(Access = private)

        function tf_sys = find_sys(self)
            % Find transfer function of 2nd order system for given data
            [prepended_out_Data, prepended_inp_Data] = prepend_zeros(self);

            % input and output might have diff dimension
            len_out_data = length(prepended_out_Data);
            len_in_data = length(prepended_inp_Data);
            if(len_out_data > len_in_data)
                prepended_out_Data = prepended_out_Data(1:len_in_data);
            elseif(len_out_data < len_in_data)
                prepended_inp_Data = prepended_inp_Data(1:len_out_data);
            end

            data = iddata(prepended_out_Data, prepended_inp_Data, self.time_sample);

            num_of_poles = 2;
            num_of_zeros = 0;
            tf_sys = tfest(data,num_of_poles,num_of_zeros);
        end

        function [prepended_out_Data, prepended_inp_Data] = prepend_zeros(self)
            % For a better functionality of 'tfest' function
            time_shift = 1;
            num_of_zeros = time_shift / self.time_sample;
            prepend_zeros = zeros(num_of_zeros,1);

            prepended_out_Data = [prepend_zeros; self.output_signal.Data];
            prepended_inp_Data = [prepend_zeros; self.input_signal.Data];
        end
    end
end
