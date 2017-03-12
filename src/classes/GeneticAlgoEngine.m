classdef GeneticAlgoEngine
    %BASICFRICTIONMODEL Summary of this class goes here
    %   Detailed explanation goes here

    properties(Access = public)
                generation_counter = 0; % for parameter scaling   
    end
    
    properties(Access = private)

        pop_distr_dec = struct(... % distribution of a new population for GA
        'best', 0.1, ...
        'best_mut', 0.5, ...
        'notbest', 0.1, ...
        'notbest_mut', 0.20, ...
        'new_rnd', 0.1); 
         
        pop_distr_size;
    
 
        pop_size;
        pop_lower_bounds;
        pop_upper_bounds
        pop_init_values;
        num_of_parameters;
        rnd_distr_objects;


    end
    

    
    methods(Access = public)
        
        % Constructors
        function self = GeneticAlgoEngine(varargin)
            
            if(nargin == 0)
                error('cam:ConstructorArgMiss',...
                    'constructor argument is missing');
            elseif(nargin == 4)
                pop_size = varargin{1};
                pop_init_values = varargin{2};
                low_bounds = varargin{3};
                up_bounds = varargin{4};
                % check args
                if(~isnumeric(pop_size) || pop_size <= 0 || mod(pop_size,10) ~= 0)
                    error('wca1:WrongConstructorsArg',...
                        ['Wrong constructors argument, mod(pop_size,'...
                        '10) has to be equal to 0 and greater than 0']);
                elseif((length(pop_init_values) ~= length(low_bounds)) ||...
                       (length(pop_init_values) ~= length(up_bounds)))
                        error('wca2:WrongConstructorsArg',...
                            ['Wrong constructors argument, length of init values'...
                            ' and lower,upper bounds has to match']);
                elseif(~all(low_bounds < pop_init_values) ||...
                        ~all(pop_init_values < up_bounds))
                    error('wca3:WrongConstructorsArg',...
                        ['Wrong constructors argument, lower bounds have to be'...
                        ' smaller than init values and init values have to be' ...
                        ' smaller than upper bounds']);
                else % args are ok
                    self.pop_size = pop_size;
                    self.pop_init_values = pop_init_values;
                    self.pop_lower_bounds = low_bounds;
                    self.pop_upper_bounds = up_bounds;
                    self.num_of_parameters = length(pop_init_values);
                    
                    self = self.init_population_size();
                    self.rnd_distr_objects = self.create_prob_distr_objects();
                end
            else
                error('wca4:WrongConstructorsArg',...
                    ['Wrong constructors argument']);
            end     
        end
        
        function  first_pop = ceate_first_population(self)
            % Generate trials for the first population
             first_pop = create_random_pop(self, self.pop_size);
        end
        
        function new_population = get_new_population(self,...
                curr_gen_errors, curr_gen_parameters, generation)
            % Generate a new population for GAlg
                        
            [best_pops, notbest_pops] =...
                self.select_pop(curr_gen_errors, curr_gen_parameters);
            
            best_mut_pops = self.mutate_pop(best_pops,...
                self.pop_distr_size.best_mut, generation);
            
            notbest_mut_pops = self.mutate_pop(notbest_pops,...
                self.pop_distr_size.notbest_mut, generation);
            
            new_rnd_pop = self.create_random_pop(self.pop_distr_size.new_rnd);
            
           
            new_population = [best_pops; notbest_pops; ...
                best_mut_pops; notbest_mut_pops; new_rnd_pop];
        end
               
        function pop_size = get_pop_size(self)
            pop_size = self.pop_size;
        end
        
        function [pop_distr_size, pop_distr_dec] = get_pop_distribution(self)
            % get distribution and size of populations
            pop_distr_size = self.pop_distr_size;   
            pop_distr_dec = self.pop_distr_dec;
        end
    end
    
    methods(Access = private)
        
        function [best_pop, notbest_pops] = select_pop(self,curr_gen_errors,...
                curr_gen_parameters)
            % select best and randomly not best pop from current
            % generation
            [~, indexes_err] = sort(curr_gen_errors);
            
            best_end_indx = self.pop_distr_size.best;
            % indexes of best populations
            best_pop_indexes = indexes_err(1:best_end_indx, :);    
            best_pop = curr_gen_parameters(best_pop_indexes, :);
            
            % randomly selected non best population
            notbest_pops = datasample(curr_gen_parameters(best_end_indx+1:end, :),...
                self.pop_distr_size.notbest);
        end
        
        function [mut_pop] =...
                mutate_pop(self, curr_population, size_of_mut_pop, gen)
            % mutates best population trials
            % gen - is used for scaling, higher generation = smaller
            % variance

            size_of_curr_pop = length(curr_population(:,1));
            mut_len = size_of_mut_pop/size_of_curr_pop;
                  
            num_of_pars = self.num_of_parameters;
            prob_distr_objs = self.rnd_distr_objects;
            mut_pop = zeros(size_of_mut_pop, num_of_pars);
                        
            for k = 1:num_of_pars
                d_obj = prob_distr_objs{k};
                slid_window = 1:mut_len;
                for p = 1:size_of_curr_pop
                    d_obj.mu = curr_population(p,k);
                    d_obj.sigma = 0.1; % / gen; % or gen const or slower decay
                    
                    mut_pop(slid_window,k) = random(d_obj, mut_len,1); 
                    slid_window = slid_window + mut_len;
                end
            end
        end
        
        function rnd_pops = create_random_pop(self, rnd_pop_size)
            % create new random populations
            
            num_of_pars = self.num_of_parameters;
            
            rnd_pops = zeros(rnd_pop_size, num_of_pars);
            rnd_objs = self.rnd_distr_objects;

            for k = 1:num_of_pars
        
                 rnd_par_values = random(rnd_objs{k}, rnd_pop_size,1);
                 rnd_pops(:, k) = rnd_par_values;           
            end          
        end
        
        function [rnd_distr_objects] = create_prob_distr_objects(self)
        % creates a probability distribution object for each parameter
        
            num_of_pars = self.num_of_parameters;
            rnd_distr_objects = cell(num_of_pars, 1);
            init_vals = self.pop_init_values;
            u_bnd = self.pop_upper_bounds;
            l_bnd = self.pop_lower_bounds;
            
            for k = 1:num_of_pars
                mean = init_vals(k);
                sig = init_vals(k); % or init_vals(k) / 2
                pd = makedist('normal','mu',mean,'sigma',sig);
                pd = pd.truncate(l_bnd(k), u_bnd(k));
                rnd_distr_objects{k} = pd;
            end       
        end
        
       
        function self = init_population_size(self)  
            % init population size
            f_names = fieldnames(self.pop_distr_dec);
            for k = 1:length(f_names)
                field = f_names{k};
                self.pop_distr_size.(field) =...
                    self.pop_distr_dec.(field) * self.pop_size;
            end
        end  
    end
end