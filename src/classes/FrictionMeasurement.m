classdef FrictionMeasurement
    %FRICTIONMEASURING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
%         sim_folder = 'sim_models';
        
        simulink_model = 'friction_ident';
        save_folder;
        
        
        switch_pulses = 1;
        switch_sinuses = 2;
        
        pulses_amp = 1:0.5:4;
        
        sin_amp_freqs =  [0.5, 3;
            1, 4;
            1.5, 4;
            2, 4;
            2.5, 4;
            3, 5.5;
            3.5, 6;
            4, 6.5;
            4.5, 7;
            5,8];
        simout;
        
    end
    
    methods(Access = public)
        % Constructor
        function self = FrictionMeasurement(varargin)
            % dataset,num_of_gen,size_of_pop,
            
            if nargin == 2
                self.pulses_amp = varargin{1};
                self.sin_amp_freqs = varargin{2};
                addpath(genpath('models'));
                addpath(genpath('measured_data'));
            end
        end
        
        function self = measure_response(self)
            % BREAK IT INTO SMALLER PARTS
            
            model = self.simulink_model;
            load_system(model);
            % init observer
            [A,B,C,D,L,X,Xhat] =  init_friction_ident_observer();
            
            set_param([model '/observer/A'], 'Gain', mat2str(A));
            set_param([model '/observer/B'], 'Gain', mat2str(B));
            set_param([model '/observer/C'], 'Gain', mat2str(C));
            set_param([model '/observer/L'], 'Gain', mat2str(L));
            set_param([model '/observer/Xhat'], 'InitialCondition', mat2str(Xhat));
            %         set_param([model '/observer/D'], 'Gain', mat2str(D));
            % init lin modeal just for testing
            set_param([model '/linear_system'],...
                'A', mat2str(A),...
                'B', mat2str(B),...
                'C', mat2str(C),...
                'D', mat2str(D),...
                'X0', mat2str(X));
            
  
            set_param([model '/observer/A'], 'Gain', mat2str(A));

            
            puls_simout = self.measure_puls_response();
            sinus_simout = self.measure_sin_response();
            
            self.simout = self.repair_data([puls_simout sinus_simout]);
            self.save_folder = get_path_of_data_folder('measurement');           
        end
        
        function save_data(self,name_arg)
            %Save measured results
            %   save folder - save folder destination
            %   measured_data - data to save
            
            % !! replace with static function, save_data !!
            
            prefix = 'friction_measurement_' ;
            if ~exist('name_arg','var')
                timestamp = strrep(datestr(now), ':', '_');
                filename = [prefix, timestamp];
            else
                filename = [prefix, name_arg];
            end
            
            file_path = [self.save_folder, '/', filename];
            simout = self.simout;
            save(file_path, 'simout');
        end
        
        function plot_results(self, file_name)
            %Plot measured results
            % file_name - for saved data, otherwise function will use data
            % from object's properties
            
            if exist('file_name','var')
                file_path = [self.save_folder, '/', file_name];
                if(exist(file_path, 'file') == 2)
                    loaded_data = load(file_path);
                    fn_data = fieldnames(loaded_data);
                    % use data from specified file
                    sims = loaded_data.(fn_data{1});
                else
                    error(['file: ', file_name, ' does not exist, '...
                        'try to specify file type (.mat, .csv, etc..)']);            
                end
            else
                sims = self.simout;
            end
            
             if isempty(sims) && ~exist('file_name','var')
                error(['Nothing to plot, measure response first or add'...
                    'filename as argument']);
             end
                   
            sims_len = length(sims);            
            n_of_cols = 3;
            n_of_rows = ceil(sims_len / n_of_cols);
            self.plot_subplt_of_responses(sims, n_of_cols, n_of_rows);
            self.plot_subplt_of_inputs(sims, n_of_cols, n_of_rows);       
        end
        
        
    end
    methods(Access = private)
        
        function simout = measure_sin_response(self)
            % Measure sinus response
            model = self.simulink_model;
            
            set_param([model '/signal_switcher/switch_val'],...
                'Value', num2str(self.switch_sinuses)); % pulses switch
            
            sin_par = self.sin_amp_freqs;
            sin_len = length(sin_par(:,1));
            
            simout(sin_len) = Simulink.SimulationOutput;
            % start from last index of simout
            for k = 1:sin_len
                sim_time = 6*pi / (sin_par(k, 2)); % sim_time needed for 6 turns
                
                amp_str = num2str(sin_par(k, 1));
                freq_str = num2str(sin_par(k, 2));
                
                fprintf(['\n\tIdentification paused, next sinus amplitude value is %s [V]',...
                    '\n\tpress enter to continue\n'], amp_str);
                pause;
                
                set_param([model '/sin_input'],...
                    'Frequency', freq_str, 'Amplitude', amp_str);
                
                simout(k) = sim(model, 'StartTime','0',...
                    'StopTime',num2str(sim_time), 'SaveOutput','on');
            end
        end
        
        function simout = measure_puls_response(self)
            
            model = self.simulink_model;
            
            
            set_param([model, '/signal_switcher/switch_val'],...
                'Value', num2str(self.switch_pulses));
            
            pulses_len = length(self.pulses_amp);
            simout(pulses_len) = Simulink.SimulationOutput;
            sim_time = 12;
            
            for k = 1:pulses_len
                amp_str = num2str(self.pulses_amp(k));
                fprintf(['\n\tIdentification paused, next value of pulse amplitude is %s [V]',...
                    '\n\tpress enter to continue\n'], amp_str);
                pause;
                
                set_param([model '/pulse_input'],...
                    'PulseWidth', '3', 'Amplitude', amp_str);
                
                %     simout_pulses(k) = sim(model, 'ReturnWorkspaceOutputs','on' );
                
                simout(k) = sim(model, 'StartTime','0',...
                    'StopTime',num2str(sim_time), 'SaveOutput','on');
            end
        end
        
        function simout = repair_data(self, sim_data)
            
            len_simout = length(sim_data);
            rep_sim = cell(len_simout, 1);
            
            
            % mozna mergnout sim vysledky uz predtimto loopem
            for k=1:len_simout
                
                % repair broken signal
                input_signal = sim_data(k).get('input_signal');
                measured_position = sim_data(k).get('measured_position');
                observed_position = sim_data(k).get('observed_position');
                
                results.input_signal = input_signal;
                results.measured_position = measured_position;
                results.measured_position_repaired =...
                    repair_measured_position(measured_position);
                results.observed_position = observed_position;
                results.observed_position_repaired =...
                    repair_measured_position(observed_position);
                
                rep_sim{k} = results;
                
            end
            
            simout = rep_sim;
            
        end
        
        function plot_subplt_of_responses(self,sims, n_of_cols,n_of_rows)
            % Plot Measured response
            sims_len = length(sims);
            figure;
            sbplt_position = subplot(n_of_rows,n_of_cols,1);
            for k = 1:sims_len
                subplot(n_of_rows,n_of_cols,k);
                plot(sims{k}.measured_position);grid on; hold on
                plot(sims{k}.measured_position_repaired);
                plot(sims{k}.observed_position);
                plot(sims{k}.observed_position_repaired);
                title(['Maximal abs value of input signal : ',...
                    num2str(round(max(sims{k}.input_signal.Data))), ' [V]']);
                
                if(k <= sims_len-n_of_cols) % turn of xlabelx
                    xlabel('');
                else
                    xlabel('Time (sec)');
                end
                if(mod(k-1,n_of_cols) ~= 0)
                    ylabel('');
                else
                    ylabel('Position (m)');
                end
            end
            
            suptitle('Measured response');
            legend(sbplt_position, 'measured', 'measured repaired',...
                'observerd', 'observed repaired');
        end
        
        function  plot_subplt_of_inputs(self,sims, n_of_cols,n_of_rows)
            % Plot Input signals
            sims_len = length(sims);
            
            figure;
            sbplt_input = subplot(n_of_rows,n_of_cols,1);
            for k = 1:sims_len
                subplot(n_of_rows,n_of_cols,k);
                inp_sig = sims{k}.input_signal;
                plot(inp_sig); grid on;
                title(['Maximal abs value : ',...
                    num2str(round(max(inp_sig.Data))), ' [V]']);
                
                if(k <= sims_len-n_of_cols) % turn of xlabelx
                    xlabel('');
                else
                    xlabel('Time (sec)');
                end
                if(mod(k-1,n_of_cols) ~= 0)
                    ylabel('');
                else
                    ylabel('Voltage (V)');
                end
            end
            
            suptitle('input signals');
            legend(sbplt_input,'input signal');
            
        end
        
    end
    
    
end

