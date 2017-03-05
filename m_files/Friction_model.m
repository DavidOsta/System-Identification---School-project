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
        tried_parameters;% rows - set of parameters, col - trials
        
%         trials_s_force_1;
%         trials_s_force_2;
%         trials_c_force_1a;
%         trials_c_force_2a;
%         trials_c_force_1b;
%         trials_c_force_2b;
%         trials_v_treshold;
        
        trials_v_treshold_range;
        test = 0;
        
        lin_system_par_struct;
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

            %parameters of linear model
            lin_sys_A = self.lin_system_par_struct.A;
            lin_sys_B = self.lin_system_par_struct.B;
            lin_sys_C = self.lin_system_par_struct.C;
            lin_sys_D = self.lin_system_par_struct.D;
            lin_sys_X = self.lin_system_par_struct.X;
            
            % Parameters for parfor cycle
            num_of_samples = length(self.dataset);
            par_trials_iters = self.num_of_trials;
            par_simulink_model = self.simulink_model;
            par_dataset = self.dataset;

            % genegeration for loop
            for gen = 1:self.num_of_gen % generations
                % init of parameters for parallel loop
                par_trial_errors = zeros(par_trials_iters, 1);
    
                % friction parameters - no struct because of 'parfor'
                friction_curr_parameters = self.tried_parameters(:,:,gen);
                
                sf1_current_gen = friction_curr_parameters(:,1);
                sf2_current_gen = friction_curr_parameters(:,2);
                cf1a_current_gen = friction_curr_parameters(:,3);
                cf2a_current_gen = friction_curr_parameters(:,4);
                cf1b_current_gen = friction_curr_parameters(:,5);
                cf2b_current_gen = friction_curr_parameters(:,6);
                vtr_current_gen = friction_curr_parameters(:,7);
                
                % parallel loop for simulink simulation of friction parametrs
                % Because of parloop I could not creat functions that could
                % be used inside loop
                parfor trial = 1:par_trials_iters
                    
                    % setup sim model - need to do it here because of parloop
                    load_system(par_simulink_model); % load simulink model
                    
                    current_friction_parameters =...
                        [sf1_current_gen(trial);...
                        sf2_current_gen(trial);...
                        cf1a_current_gen(trial);...
                        cf2a_current_gen(trial);...
                        cf1b_current_gen(trial);...
                        cf2b_current_gen(trial);...
                        vtr_current_gen(trial)];
                    
                    set_param([par_simulink_model '/friction_parameters'],...
                        'Value', mat2str(current_friction_parameters));
                           
                    set_param([par_simulink_model '/linear_system'],...
                        'A', mat2str(lin_sys_A),...
                        'B', mat2str(lin_sys_B),...
                        'C', mat2str(lin_sys_C),...
                        'D', mat2str(lin_sys_D),...
                        'X0', mat2str(lin_sys_X));


                 % the workers in a parallel pool cannot start or access further parallel pools
                 % thus simple loop has to be used
%                  sample_errors = zeros(num_of_samples,1);
                 sample_errors = 0;
                             % for loop - for each measured input signal
                             for k = 1:num_of_samples
                                 % Load of sim model
                                 % because of parallel for loop
                                 
                                 dataset_struct = par_dataset{k}; % get measurement results
                                 input_signal = dataset_struct.force_input;
                                 input_signal_mat = [input_signal.Time input_signal.Data];

                                 % I wonder which of these bad ways how to
                                 % load timeseries is better...
                                 set_param([par_simulink_model '/input_signal'],...
                                     'VariableName', mat2str(input_signal_mat));
                                 
%                                  set_param([par_simulink_model '/input_signal'],...
%                                      'VariableName', 'input_signal');
%                                  assignin('base','input_signal',input_signal);
                                 
                                 measured_position = dataset_struct.measured_position;
                                 sim_time = measured_position.Time(end);
                                 simout = sim(par_simulink_model,'StartTime','0','StopTime',num2str(sim_time),...
                                     'SaveOutput','on');
                                 sim_position = simout.get('cart_response');
                                 
                                 % mean square error
                                 error = sqrt(mean((measured_position.Data - sim_position.Data).^2));
%                                  sample_errors(k) = error;
                                 sample_errors = sample_errors + error;
                              end
                    
                    par_trial_errors(trial) = sample_errors;
                        
                end
                
                 self.models_errors(:, gen) = par_trial_errors;

                self = self.update_population(gen);

                self.print_current_top_pop(gen);
            end
        end
        
        function best_models = get_best_models(self,num_of_best)
            % select friction models with smallest error
            if nargin < 1
                num_m = self.num_of_best;
            else
                num_m = num_of_best;
            end
            
            [sorted_err, indexes_err] = sort(self.models_errors);
            
            best_models = cell(num_m,1);
            for k = 1:num_m
                
%                 best_models{k} = struct(...
%                     'sf1',self.trials_s_force_1(indexes_err(k)),...
%                     'sf2',self.trials_s_force_2(indexes_err(k)),...
%                     'cf1a',self.trials_c_force_1a(indexes_err(k)),...
%                     'cf2a',self.trials_c_force_2a(indexes_err(k)),...
%                     'cf1b',self.trials_c_force_1b(indexes_err(k)),...
%                     'cf2b',self.trials_c_force_2b(indexes_err(k)),...
%                     'vtr',self.trials_v_treshold(indexes_err(k)),...
%                     'error',sorted_err(k));
            end
        end
        
        function print_sim_results(self)
        % Plot simulation results
            gen_range = 1:1:self.num_of_gen;
            
            figure;
            subplot 211;
            plot(gen_range, self.models_errors(1,:),... % 1st top error
                 gen_range, self.models_errors(2,:)) % 2nd top error
            grid on;
            title('Top models errors vs generation');
            xlabel('Generation');
            ylabel('Mean square Error');
            legend('Smallest error during simulation','2nd smallest error');
            set(gca,'xtick', gen_range)
            
            
            subplot 212;
            best_param_of_gens =...
                permute(self.tried_parameters(1,:,:), [3 2 1]);
            
            parameters = { 'sf1','sf2','cf1a','cf2a',...
                'cf1b','cf2b','vtr'};
            par_len = length(parameters);
            
            for p = 1:par_len
                plot(gen_range,best_param_of_gens(:,p));hold on; grid on
            end
            title('Parameters of top model in each generation');
            xlabel('Generation');
            ylabel('Parameters');
            legend(parameters);
            set(gca,'xtick', gen_range);
            
            figure;
            subplot(par_len,1,1);
            concat_errors = self.models_errors(:); % concat generations
          
            sub_plots = ceil(par_len/2); 
            for p = 1:par_len
                subplot(sub_plots, 2, p);
                % get parameters for each trial and generation
                p_par = permute(self.tried_parameters(:,p,:), [1 3 2]);
                % concat & sort parameters 
                [sorted_parameters, sorted_indexes] = sort(p_par(:));
                
                plot(sorted_parameters, concat_errors(sorted_indexes)); grid on;
                title(['Error vs Parameter', ' ', parameters(p)]);
                xlabel('Parameters value');
                ylabel('Models error');
            end
        end
    end
    
    
    
    % Private methods
    methods(Access = private)
        
        function  self = generate_init_parameters(self)
            % Generate trials for the first population
         
            sf1 = self.init_s_force_1 + (-0.5 + (1)*rand(self.num_of_trials,1));
            sf2 = self.init_s_force_2 + (-0.5 + (1)*rand(self.num_of_trials,1));
            cf1a = self.init_c_force_1a + (-0.05 + (0.5)*rand(self.num_of_trials,1));
            cf2a = self.init_c_force_2a + (-0.05 + (0.5)*rand(self.num_of_trials,1));
            cf1b = self.init_c_force_1b +  (-0.05 + (0.5)*rand(self.num_of_trials,1));
            cf2b = self.init_c_force_2a + (-0.05 + (0.5)*rand(self.num_of_trials,1));
            
            self.trials_v_treshold_range = transpose(linspace(0, 0.95, self.num_of_trials));
            vtr = self.shuffle_array(self.trials_v_treshold_range);
            
       
            self.models_errors = zeros(self.num_of_trials, self.num_of_gen);
            % row - friction parameters for each trial, column - array of
            % parameters
            self.tried_parameters = [sf1, sf2,...
                cf1a, cf2a,...
                cf1b, cf2b,...
                vtr];
            
            self.lin_system_par_struct = init_parametrs();

        end
        
        function shufled_arr = shuffle_array(self, arr)
            % shuffle array 'trials_v_treshold'
            shufled_arr = arr(randperm(numel(arr)));
        end
               
         
       
        function self = update_population(self,gen)
            % Evaluate current population of model parameters
            % select top 5 errors and their parameters

            
            if(gen < self.num_of_gen) % no need to update population after last generation
                best_trials = self.get_best_trials(gen);
                mutated_trials = self.mutate_best_trials(best_trials, gen);
                
                
                % !!!! this will be class propertie created
                parameters = { 'sf1','sf2','cf1a','cf2a',...
                    'cf1b','cf2b','vtr'};
                % !!!
                len_par = length(parameters);
                
                par_matrix =...
                    zeros(self.num_of_trials, len_par);
                for p = 1:len_par
                    par_matrix(:,p) =...
                        [best_trials.(parameters{p});
                        mutated_trials.(parameters{p})];
                end
                
                self.tried_parameters(:,:, gen + 1) = par_matrix;
              
            end
            
        end
        
        function best_trials = get_best_trials(self,gen)
            % select best trials from current population           
            [sorted_err, indexes_err] = sort(self.models_errors(:,gen)); 
            % indexes of best populations
            best_pop_indxs = indexes_err(1:self.num_of_best);
     
            % !!!! this will be class propertie created
            parameters = { 'sf1','sf2','cf1a','cf2a',...
                'cf1b','cf2b','vtr'};
            % !!! 
            len_par = length(parameters);
            
            % Dynamically create a struct
            for par = 1:len_par
                best_par = self.tried_parameters(best_pop_indxs,par, gen);
                best_trials.(parameters{par}) =...
                    self.replace_negative_elements(best_par);
            end

        end
        
        function mutated_trials = mutate_best_trials(self, best_trials, gen)
            % mutates best population trials
            % gen - is used for scaling, higher generation = smaller
            % variance
         
            function altered_arr= mutate(arr, mod)
                %helper function
%                 mutated_tr =...
%                     repmat(trials,len_repmat,1) + randn(num_of_mut,1) / gen;

                % normal distr w/ mean = 0, variance = 0.1 / gen (scaling)
                
                if(strcmp(mod, 'mutate'))
                    arr =...
                        repmat(arr,len_repmat,1) + normrnd(0, 0.1 / gen,num_of_mut,1);
                    altered_arr = self.replace_negative_elements(arr);
                elseif(strcmp(mod, 'shuffle'))
                    arr = self.trials_v_treshold_range(1:num_of_mut);
                    altered_arr = self.shuffle_array(arr);
                else
                    altered_arr = arr;
                end
            end
 
            num_of_mut = self.num_of_trials - self.num_of_best;
            len_repmat = num_of_mut/self.num_of_best;

            f_names = fieldnames(best_trials);
            % create struct
            for p = 1:length(f_names)
                par = f_names{p};
                if(strcmp(par,'vtr'))
                    mutated_trial = mutate(best_trials.(par), 'shuffle');
                else
                    mutated_trial = mutate(best_trials.(par), 'mutate');
                end
                mutated_trials.(par) = mutated_trial;
            end
        
        end
        
        function positive_arr = replace_negative_elements(self,arr)
            % change negative elements to positive of same absolute size
            negative_indexes = arr < 0;
            if(any(negative_indexes) == 1 ) % has neg value
                arr(negative_indexes) = arr(negative_indexes) * -1;
            end
            positive_arr = arr;
        end
        
        
        function print_current_top_pop(self,gen)
            % Print model parameters of best population in current
            % generation
            self.tried_parameters(:,:, gen);
            [min_err, min_indx] = min(self.models_errors(:, gen));
            
            fprintf('---- Best population of generation n. %f ----\n', gen);
            fprintf('Static force 1 = %f \n',...
                self.tried_parameters(min_indx, 1 , gen));
            fprintf('Static force 2 = %f \n',...
                self.tried_parameters(min_indx, 2 , gen));
            fprintf('Coloumb force 1a = %f \n',...
                self.tried_parameters(min_indx, 3 , gen));
            fprintf('Coloumb force 2a = %f \n',...
                self.tried_parameters(min_indx, 4 , gen));
            fprintf('Coloumb force 1b = %f \n',...
                self.tried_parameters(min_indx, 5 , gen));
            fprintf('Coloumb force 2b = %f \n',...
                self.tried_parameters(min_indx, 6 , gen));
            fprintf('Velocity threshold = %f \n',...
                self.tried_parameters(min_indx, 7 , gen));
            fprintf('Total Mean square Error %f \n',min_err(1));
            fprintf('----------------------------------------\n');
        
        end
        
       
        
    end

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
 
                     % Load variables to worker's base workspace
                    % cannot set up it via for loop because of PARFOR !!!
%                     set_param([par_simulink_model '/friction_parameters' '/sf1'],...
%                         'Value', num2str(sf1_current_gen(trial)));
%                     set_param([par_simulink_model '/friction_parameters' '/sf2'],...
%                         'Value', num2str(sf2_current_gen(trial)));
%                     set_param([par_simulink_model '/friction_parameters' '/cf1a'],...
%                         'Value', num2str(cf1a_current_gen(trial)));
%                     set_param([par_simulink_model '/friction_parameters' '/cf2a'],...
%                         'Value', num2str(cf2a_current_gen(trial)));
%                     set_param([par_simulink_model '/friction_parameters' '/cf1b'],...
%                         'Value', num2str(cf1b_current_gen(trial)));
%                     set_param([par_simulink_model '/friction_parameters' '/cf2b'],...
%                         'Value', num2str(cf2b_current_gen(trial)));
%                     set_param([par_simulink_model '/friction_parameters' '/vtr'],...
%                         'Value', num2str(vtr_current_gen(trial)));
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 