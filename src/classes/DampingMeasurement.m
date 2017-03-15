classdef DampingMeasurement
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = private)
        
        simulink_model = 'damping_measurement';
        simout;
        
        models_errors;
        tested_parameters;
        file_path;
    end
    
    methods(Access = public)
        
        % Constructor
        function self = DampingMeasurement(varargin)
            add_paths();  
        end
        
        function self = measure_response(self)
            % BREAK INTO SMALLER FUNCTIONS
            sim_model = self.simulink_model;
            
            % Unfortunately this has to be done manually
            fprintf(['\n \tRun simulink manually, and excite pulses on system\n',...
                '\twhen finished press enter to continue\n']);
            open_system(sim_model);
            ws_positon = 'position';
            ws_angle = 'angle';
            ws_filtered_angle = 'filtered_angle';
            ws_impulse_signal = 'impulse_signal';
            set_param('damping_measurement/position', 'VariableName', ws_positon);
            set_param('damping_measurement/angle', 'VariableName', ws_angle);
            set_param('damping_measurement/filtered_angle', 'VariableName', ws_filtered_angle);
            set_param('damping_measurement/impulse_signal', 'VariableName', ws_impulse_signal);

            
            step_time = 0.001;
            period = 100;
            pulse_width = (period / 100) * step_time;
            
            set_param(['damping_measurement' '/zoh'],...
                'SampleTime', num2str(step_time));
            set_param(['damping_measurement' '/input_pulse'],...
                'Amplitude', num2str(1/step_time));
            set_param(['damping_measurement' '/input_pulse'],...
                'PulseWidth', num2str(pulse_width));
            set_param(['damping_measurement' '/input_pulse'],...
                'Period', num2str(period));
             set_param(['damping_measurement' '/input_pulse'],...
                'PhaseDelay', num2str(0));
            
            fprintf('\n \tProgram is paused, press enter\n');
            
            pause;
            
            % little hack
            base_ws_position = evalin('base',ws_positon,'0');
            base_ws_angle = evalin('base',ws_angle,'0');
            base_ws_impulse_signal = evalin('base',ws_impulse_signal,'0');
            %             base_ws_filtered_angle = evalin('base',ws_filtered_angle,'0');
            
            if(isa(base_ws_position,'timeseries') &&...
                    isa(base_ws_angle,'timeseries'))                
                measured_position = base_ws_position;
                measured_angle = base_ws_angle;
                impulse_signal = base_ws_impulse_signal;
            else
                fprintf(['\n \tYou have forgotten to run simulink'...
                    ' & execute measurement\n']);
                bdclose(sim_model);
                self = self.measure_response();
            end
            
            self.simout.position = self.adjust_measured_data(measured_position);
            self.simout.angle = self.adjust_measured_data(measured_angle);
            self.simout.impulse_signal = impulse_signal;

            self.file_path = save_data(self.simout, 'data/measured_data',...
                'damping_measurement');
        end
        
        
        function results_path = get_results_path(self)
            % return path to measured and saved data for the object
            if(isempty(self.file_path))
                error('You need to measure damping first');
            else
                results_path = self.file_path;
            end
        end
        
    end
    
    
    methods(Access = private)
        
        function adjusted_data = adjust_measured_data(self, measured_data)
            % change time axes of measured signal
            
            figure;
            plot(measured_data);grid on;
            
            fprintf(['\n\tSelect new beginning and end of measured signal so that'...
                'response starts from cca time 0\n']);
            
            while(true)
                input_start_time =...
                    input('Start time (number) must be greater than 0 >>> ');
                input_end_time =...
                    input(['End time (number) must be smaller than',...
                    num2str(measured_data.Time(end)) ,' >>> ']);
                
                sample_time = measured_data.Time(2) - measured_data.Time(1);
                first_sample = input_start_time / sample_time + 1;
                last_sample = input_end_time / sample_time - 1;
                
                adj_data_time = measured_data.Time(first_sample:last_sample) - input_start_time ;
                adj_data_data = measured_data.Data(first_sample:last_sample);
                
               
                adjusted_data = timeseries(adj_data_data, adj_data_time);
                
                figure;
                plot(adjusted_data);grid on;
                
                try_again = input('Do you want to adjust it again? [y/n] >>> ','s');
                
                if(strcmp(try_again,'n'))
                    close;
                    break;
                end
                close;
                
            end
            close;
        end
    end
end