 classdef FrictionModel
    %PLOTER Summary of self class goes here
    %   get friction model
    
     
    properties(Access = private)
        simulink_model = 'sim_friction_model';
        
        num_of_gen = 10; % number of genetic algo's populations
        size_of_pop = 10 ;% 50, 100 % number of genetic algos's samples in one population
        
        %TODO - this won't be here
        coloumb_model_init_vals = struct(...
            'fc1a', 0.8,...
            'fc2a', 0.8,...
            'fc1b', 0.4,...
            'fc2b', 0.4,...
            'vtr', 0.2);
        
        stiction_model_init_vals = struct(...
            'fs1', 0.8,...
            'fs2', 0.8,...
            'fc1a', 0.2,...
            'fc2a', 0.2,...
            'fc1b', 0.4,...
            'fc2b', 0.4,...
            'vtr', 0.2);
            
        striberck_model_init_vals = struct('faf', 12, 'vafa', 3);


        dataset;
        models_errors;
        tested_parameters;% rows - set of parameters, col - trials
        
        ga_engine; % !!111 REDOOOO
                  
        lin_system_par_struct;
    end
    
    
    methods(Access = public)
        % Constructor
        function self = FrictionModel(varargin)
            % dataset,num_of_gen,size_of_pop,
            
            %TODO - TEMPORARY CHANGE THIS
            init_vals = [0.8; 0.8; 0.1;  0.1;  0.4; 0.4; 0.4];
            l_bounds =  [0;   0;   0;    0;    0;   0;   0];
            u_bounds =  [3;   3;   2;    2;    2;   2;   0.8];

            
            self.ga_engine = GeneticAlgoEngine(self.size_of_pop,init_vals,...
                l_bounds,u_bounds);

            if nargin > 1
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

            %TODO - Break it to smaller parts + optimize simulink model
            
            %parameters of linear model
            lin_sys_A = self.lin_system_par_struct.A;
            lin_sys_B = self.lin_system_par_struct.B;
            lin_sys_C = self.lin_system_par_struct.C;
            lin_sys_D = self.lin_system_par_struct.D;
            lin_sys_X = self.lin_system_par_struct.X;
            
            % Parameters for parfor cycle
            num_of_samples = length(self.dataset);
            par_size_of_pop = self.size_of_pop;
            par_simulink_model = self.simulink_model;
            par_dataset = self.dataset;
            

            % genegeration for loop
            for gen = 1:self.num_of_gen % generations
                % init of parameters for parallel loop
                par_pop_errors = zeros(par_size_of_pop, 1);
    
                % friction parameters - no struct because of 'parfor'
                curr_gen_pop_parameters = self.tested_parameters(:,:,gen);
                
                % parallel loop for simulink simulation of friction parametrs
                % Because of parloop I could not creat functions that could
                % be used inside loop
                parfor pop = 1:par_size_of_pop
                    
                    % setup sim model - need to do it here because of parloop
                    load_system(par_simulink_model); % load simulink model
                    
                    current_friction_parameters = ...
                        curr_gen_pop_parameters(pop,:)';
                    
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

                 sample_errors = 0;
                             % for loop - for each measured input signal
                             for k = 1:num_of_samples
                                 
                                 dataset_struct = par_dataset{k}; % get measurement results
                                 input_signal = dataset_struct.force_input;
                                 input_signal_mat = [input_signal.Time input_signal.Data];

                                 %TODO - there must be a better way
                                  set_param([par_simulink_model '/input_signal'],...
                                     'VariableName', mat2str(input_signal_mat));
                                                              
                                 measured_position = dataset_struct.measured_position;
                                 sim_time = measured_position.Time(end);
                                 simout = sim(par_simulink_model,'StartTime','0','StopTime',num2str(sim_time),...
                                     'SaveOutput','on');
                                 sim_position = simout.get('cart_response');
                                 
                                 % mean square error
                                 error = sqrt(mean((measured_position.Data - sim_position.Data).^2));
                                 sample_errors = sample_errors + error;
                              end
                    
                    par_pop_errors(pop) = sample_errors;
                    
                    % clean up
                    bdclose('all');
                end
                           
                 self.models_errors(:, gen) = par_pop_errors;
             
                if(gen < self.num_of_gen) % no need to update after last gen
                    self.tested_parameters(:,:, gen + 1) =...
                        self.ga_engine.get_new_population(par_pop_errors,...
                        curr_gen_pop_parameters, gen);
                end
                
                self.print_current_top_pop(gen);
            end
        end
        

        
        function save_friction_model(self)
            
          %TODO - Implement

        % save simulation results
%             model_erros = self.models_errors;
%             num_of_gen = self.num_of_gen;
%             model_par_names = self. ;
%             all_tried_par = self.tested_parameters;
%         
        end
        
        function print_sim_results(self)
        % Plot simulation results
            gen_range = 1:1:self.num_of_gen;
            gen_len = self.num_of_gen;
            
            [sorted_errors, err_indexes] = sort(self.models_errors);

            figure;
            subplot 211;
            plot(gen_range, sorted_errors(1,:),... % 1st top error
                 gen_range, sorted_errors(2,:)) % 2nd top error
            grid on;
            title('Top models errors vs generation');
            xlabel('Generation');
            ylabel('Mean square Error');
            legend('Smallest error during simulation','2nd smallest error');
            set(gca,'xtick', gen_range)
            
            
            subplot 212;
            
            
            parameters = { 'sf1','sf2','cf1a','cf2a',...
                'cf1b','cf2b','vtr'};
            par_len = length(parameters);
            
            top_parameters = zeros(gen_len,par_len);
            
            for gen = 1:gen_len
                  
                index_of_best_par = err_indexes(1,gen);
                top_par_in_gen =...
                    self.tested_parameters(index_of_best_par,:,gen);
                top_parameters(gen,:) = top_par_in_gen;
            end

            
            for p = 1:par_len
                plot(gen_range,top_parameters(:,p));hold on; grid on
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
                p_par = permute(self.tested_parameters(:,p,:), [1 3 2]);
                % concat & sort parameters 
                [sorted_parameters, sorted_indexes] = sort(p_par(:));
                % plot models error for a given parameter
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
                
            self.models_errors = zeros(self.size_of_pop, self.num_of_gen);

            self.tested_parameters =...
                self.ga_engine.ceate_first_population();

            self.lin_system_par_struct = init_parametrs();

        end
        

 
        function print_current_top_pop(self,gen)
            % Print model parameters of best population in current
            % generation
            self.tested_parameters(:,:, gen);
            [min_err, min_indx] = min(self.models_errors(:, gen));
            
            fprintf('---- Best population of generation n. %f ----\n', gen);
            fprintf('Static force 1 = %f \n',...
                self.tested_parameters(min_indx, 1 , gen));
            fprintf('Static force 2 = %f \n',...
                self.tested_parameters(min_indx, 2 , gen));
            fprintf('Coloumb force 1a = %f \n',...
                self.tested_parameters(min_indx, 3 , gen));
            fprintf('Coloumb force 2a = %f \n',...
                self.tested_parameters(min_indx, 4 , gen));
            fprintf('Coloumb force 1b = %f \n',...
                self.tested_parameters(min_indx, 5 , gen));
            fprintf('Coloumb force 2b = %f \n',...
                self.tested_parameters(min_indx, 6 , gen));
            fprintf('Velocity threshold = %f \n',...
                self.tested_parameters(min_indx, 7 , gen));
            fprintf('Total Mean square Error %f \n',min_err(1));
            fprintf('----------------------------------------\n');
        
        end
        
       
        
    end

 end

 
 
 
 
 
 
 
 
 