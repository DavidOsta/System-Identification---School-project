 classdef Friction_model
    %PLOTER Summary of self class goes here
    %   get friction model
    
    
    properties(Access = public)
    end
    
    properties(Access = private)
        simulink_model = 'sim_friction_model';
        init_s_force_1 = 0.8; % static friction force
        init_s_force_2 = 0.8;
        init_c_force_1a = 0.1; % coloumb friction force
        init_c_force_2a = 0.1;
        init_c_force_1b = 0.4; % coloumb friction force for velocity under threshold
        init_c_force_2b = 0.4;
        init_c_velocity_threshold = 0.4; % threshold for coloumb force
        num_of_gen = 10; % number of genetic algo's populations
        num_of_trials = 50 ; % number of genetic algos's samples in one population
        num_of_best = 10; % number of selected best samples of population
        
        dataset;
        models_errors;
        trials_s_force_1;
        trials_s_force_2;
        trials_c_force_1a;
        trials_c_force_2a;
        trials_c_force_1b;
        trials_c_force_2b;
        trials_v_treshold_range;
        trials_v_treshold;
        test = 0;
    end
    
    
    methods
        % Constructors
        function self = Friction_model(varargin)
            % dataset,num_of_gen,num_of_trials,num_of_best_pop are given
            if nargin > 3
                self.dataset = varargin{1};
                self.num_of_gen = varargin{2};
                self.num_of_trials = varargin{3};
                self.num_of_best = varargin{4};
                self = self.generate_init_parameters();
                % dataset,num_of_gen,num_of_trials are given
            elseif nargin > 2
                self.dataset = varargin{1};
                self.num_of_gen = varargin{2};
                self.num_of_trials = varargin{3};
                self = self.generate_init_parameters();
                % dataset,num_of_gen
            elseif nargin > 1
                self.dataset = varargin{1};
                self.num_of_gen = varargin{2};
                self = self.generate_init_parameters();
                % only dataset is given
            elseif nargin > 0
                self.dataset = varargin{1};
                self = self.generate_init_parameters();
            else
                error('constructor argument is missing')
            end
        end
        
    end
    methods(Access = public)
        function self = find_friction_model(self)
            num_of_samples = length(self.dataset);
            
            for gen = 1:self.num_of_gen % generations
                %parfor
                for trial = 1:self.num_of_trials
                    trial_error =...
                        self.get_trial_error(trial,num_of_samples);
                    
                    self.models_errors(trial,1) = trial_error;
                        
                end
                self.print_current_top_pop(gen);
                self = self.update_population(gen);
            end
        end
        
        function best_models = get_best_models(self,num)
            if nargin < 1
                num_m = self.num_of_best;
            else
                num_m = num;
            end
            
            [sorted_err, indexes_err] = sort(self.models_errors);
            
            best_models = cell(num_m,1);
            for k = 1:num_m
                
                best_models{k} = struct(...
                    'sf1',self.trials_s_force_1(indexes_err(k)),...
                    'sf2',self.trials_s_force_2(indexes_err(k)),...
                    'cf1a',self.trials_c_force_1a(indexes_err(k)),...
                    'cf2a',self.trials_c_force_2a(indexes_err(k)),...
                    'cf1b',self.trials_c_force_1b(indexes_err(k)),...
                    'cf2b',self.trials_c_force_2b(indexes_err(k)),...
                    'vtr',self.trials_v_treshold(indexes_err(k)),...
                    'error',sorted_err(k));
            end
        end
        
        
    end
    
    
    
    % Private methods
    methods(Access = private)
        
        function  self = generate_init_parameters(self)
            % Generate trials for the first population
            self.trials_s_force_1 = self.init_s_force_1 + (-0.5 + (1)*rand(self.num_of_trials,1));
            self.trials_s_force_2 = self.init_s_force_2 + (-0.5 + (1)*rand(self.num_of_trials,1));
            self.trials_c_force_1a = self.init_c_force_1a + (-0.05 + (0.5)*rand(self.num_of_trials,1));
            self.trials_c_force_2a = self.init_c_force_2a + (-0.05 + (0.5)*rand(self.num_of_trials,1));
            self.trials_c_force_1b = self.init_c_force_1b +  (-0.05 + (0.5)*rand(self.num_of_trials,1));
            self.trials_c_force_2b = self.init_c_force_2a + (-0.05 + (0.5)*rand(self.num_of_trials,1));
            self.models_errors = zeros(self.num_of_trials,1);
            
            self.trials_v_treshold_range = transpose(linspace(0, 0.95, self.num_of_trials));
            self.trials_v_treshold = self.shuffle_array(self.trials_v_treshold_range);
            
        end
        
        function shufled_arr = shuffle_array(self, arr)
            % shuffle array 'trials_v_treshold'
            shufled_arr = arr(randperm(numel(arr)));
        end
        
        function errors_sum = get_trial_error(self, trial, num_of_samples)
            
            % because of parallelization store values in the vector instead
            % of simple addition of one variable
            errors = zeros(num_of_samples,1);
            
            % THIS NEED PARALELIZATION
            for k = 1:num_of_samples
                % Load of sim model 
                % because of parallel for loop
                load_system(self.simulink_model); % load simulink model

                dataset_struct = self.dataset{k}; % get measurement results
                
                self.update_base_workspace(trial, dataset_struct.force_input);
                
                measured_position = dataset_struct.measured_position;

                sim_time = measured_position.Time(end);
                sim_position = self.run_simulation(sim_time);
                
                errors(k)= self.calculate_sample_error(...
                    measured_position.Data, sim_position.Data);
            end
            errors_sum = sum(errors);
        end
        
        function update_base_workspace(self,trial, input_signal)
            % update matlab base workspace for simulink model
            % because of parallelization
            friction_parametrs = struct(...
                'sf1', self.trials_s_force_1(trial),...
                'sf2', self.trials_s_force_2(trial),...
                'cf1a', self.trials_c_force_1a(trial),...
                'cf2a', self.trials_c_force_2a(trial),...
                'cf1b', self.trials_c_force_1b(trial),...
                'cf2b', self.trials_c_force_2b(trial),...
                'vtr',self.trials_v_treshold(trial) );
 
            assignin('base','input_signal',input_signal);
            assignin('base','friction_parametrs' ,friction_parametrs);
            
            [A,B,C,D,L,X,Xhat] = init_parametrs();
            assignin('base','A' ,A);
            assignin('base','B' ,B);
            assignin('base','C' ,C);
            assignin('base','D' ,D);
            assignin('base','L' ,L);
            assignin('base','X' ,X);
            assignin('base','Xhat' ,Xhat);


        
        end
        
        function error = calculate_sample_error(self,measured,simulated)
            % Calculate mean square error of measured / simulated timeseries
            error = sqrt(mean((measured - simulated).^2));
            % normalized rmse
            %max_observed = max(max(simulated_position.Data), max(measured_position.Data));
            %min_obserbed = min(min(simulated_position.Data), min(measured_position.Data));
            % nms_err = ms_err / (max_observed-min_obserbed);
        end
        
        function sim_position = run_simulation(self,sim_time)
            %PLOT_IMPULSE_IDENT Summary of this function goes here
            %   Run simulation of system with given friction model,
     
            simout = sim(self.simulink_model,'StartTime','0','StopTime',num2str(sim_time),...
                'SaveOutput','on');
            
            sim_position = simout.get('cart_response');
        end
        
        function self = update_population(self,gen)
            % Evaluate current population of model parameters
            % select top 5 errors and their parameters
            
            best_trials = self.get_best_trials();
            mutated_trials = self.mutate_best_trials(best_trials, gen);
            
            % update generation
            self.trials_s_force_1 = [best_trials.sf1; mutated_trials.sf1];
            self.trials_s_force_2 = [best_trials.sf2; mutated_trials.sf2];
            self.trials_c_force_1a = [best_trials.cf1a; mutated_trials.cf1a];
            self.trials_c_force_2a = [best_trials.cf2a; mutated_trials.cf2a];
            self.trials_c_force_1b = [best_trials.cf1b; mutated_trials.cf1b];
            self.trials_c_force_2b = [best_trials.cf2b; mutated_trials.cf2b];
            self.trials_v_treshold = [best_trials.vtr; mutated_trials.vtr];
            
        end
        
        function best_trials = get_best_trials(self)
            % select best trials from current population
            [sorted_err, indexes_err] = sort(self.models_errors); % sort from small to large
            
            % indexes of best populations
            best_pop_indxs = indexes_err(1:self.num_of_best);
            % pull out best trials
            best_trials = struct(...
                'sf1',self.trials_s_force_1(best_pop_indxs),...
                'sf2',self.trials_s_force_2(best_pop_indxs),...
                'cf1a',self.trials_c_force_1a(best_pop_indxs),...
                'cf2a',self.trials_c_force_2a(best_pop_indxs),...
                'cf1b',self.trials_c_force_1b(best_pop_indxs),...
                'cf2b',self.trials_c_force_2b(best_pop_indxs),...
                'vtr',self.trials_v_treshold(best_pop_indxs));
        end
        
        function mutated_trials = mutate_best_trials(self, best_trials, gen)
            % mutates best population trials
            % gen - is used for scaling, higher generation = smaller
            % changes
            num_of_mut = self.num_of_trials - self.num_of_best;
            len_repmat = num_of_mut/self.num_of_best;
            function mutated_trials = mutate(trials)
                %helper function
                mutated_trials =...
                    repmat(trials,len_repmat,1) + randn(num_of_mut,1) / gen;
            end
            shuffled_vtr =...
                self.shuffle_array(self.trials_v_treshold_range(1:num_of_mut));
            
            mutated_trials = struct(...
                'sf1',mutate(best_trials.sf1),...
                'sf2',mutate(best_trials.sf2),...
                'cf1a',mutate(best_trials.cf1a),...
                'cf2a',mutate(best_trials.cf2a),...
                'cf1b',mutate(best_trials.cf1b),...
                'cf2b',mutate(best_trials.cf2b),...
                'vtr',shuffled_vtr);
            
        end
        
        function print_current_top_pop(self,gen)
            % Print model parameters of best population in current
            % generation
            [sorted_err, indexes_err] = sort(self.models_errors);
            
            fprintf('---- Best population of generation n. %f ----\n', gen);
            fprintf('Static force 1 = %f \n',...
                self.trials_s_force_1(indexes_err(1)));
            fprintf('Static force 2 = %f \n',...
                self.trials_s_force_2(indexes_err(1)));
            fprintf('Coloumb force 1a = %f \n',...
                self.trials_c_force_1a(indexes_err(1)));
            fprintf('Coloumb force 2a = %f \n',...
                self.trials_c_force_2a(indexes_err(1)));
            fprintf('Coloumb force 1b = %f \n',...
                self.trials_c_force_1b(indexes_err(1)));
            fprintf('Coloumb force 2b = %f \n',...
                self.trials_c_force_2b(indexes_err(1)));
            fprintf('Velocity threshold = %f \n',...
                self.trials_v_treshold(indexes_err(1)));
            fprintf('Total Mean square Error %f \n',sorted_err(1));
            fprintf('----------------------------------------\n');
            
        end
        
        
        %         function create_figure(self, sufix)
        %             fig_name = strcat(self.path,'-', sufix);
        %             figure('Name',fig_name,'NumberTitle','off');
        %         end
        %
        %         function file_name = set_file_name(self)
        %             file_struct = dir(fullfile(self.folder_path,'*.mat'));
        %             file_name = file_struct.name;
        %         end
        %
        %         function data = load_data(self, file_name)
        %             file_path = strcat(self.folder_path,'/',file_name);
        %             data = load(file_path); % load to struct
        %         end
        
        
    end

end