classdef GeneticAlgoEngineTest < matlab.unittest.TestCase
    
    properties (TestParameter)
        
        
    end
    
    
    methods (Test)
        function test_constructor(testCase)
            
            exp_pop_size = 100;
            init_vals = [0.8; 0.8; 0.1; 0.1; 0.4; 0.4; 0.4];
            l_bounds = zeros(length(init_vals),1);
            u_bounds = ones(length(init_vals),1);
            
            fm_obj = GeneticAlgoEngine(exp_pop_size,init_vals,l_bounds,u_bounds);
            
            act_num_of_trials = fm_obj.get_pop_size();
            testCase.verifyEqual(act_num_of_trials,exp_pop_size);
            
            testCase.verifyError(@()GeneticAlgoEngine(-100, ones(3,1), ones(3,1), ones(3,1)),...
                'wca1:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(0, ones(3,1), ones(3,1), ones(3,1)),...
                'wca1:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine('kek', ones(3,1), ones(3,1), ones(3,1)),...
                'wca1:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(100, ones(20,1), ones(21,1), ones(20,1)),...
                'wca2:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(100, ones(20,1), ones(20,1), ones(22,1)),...
                'wca2:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(100, ones(3,1), [0;0;0], [2;2;1]),...
                'wca3:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(100, ones(3,1), [1;0;0], [2;2;2]),...
                'wca3:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(100, ones(3,1), [1;1;1], [1;1;1]),...
                'wca3:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(122,255),...
                'wca4:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(-100, [1; 1]),...
                'wca4:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(-100, [1; 1], [2; 2]),...
                'wca4:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine('hello'),...
                'wca4:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(95),...
                'wca4:WrongConstructorsArg');
            testCase.verifyError(@()GeneticAlgoEngine(),...
                'cam:ConstructorArgMiss');
        end
        
        function test_init_of_attributes(testCase)
            % test of class atributes
            function struct_sum = sum_fields_of_struct(strct)
                % Helper function
                f_names = fieldnames(strct);
                f_len = length(f_names);
                f_vals = zeros(f_len,1);
                for k = 1:f_len
                    f_vals(k) = strct.(f_names{k});
                end
                struct_sum = sum(f_vals);
            end
            
            
            
            size_of_pop = 1000;
            init_vals = [0.8; 0.8; 0.1; 0.1; 0.4; 0.4; 0.4];
            l_bounds = zeros(length(init_vals),1);
            u_bounds = ones(length(init_vals),1);
            
            fm_obj = GeneticAlgoEngine(size_of_pop,init_vals,...
                l_bounds,u_bounds);
            
            [pop_distr_size, pop_distr_dec] =...
                fm_obj.get_pop_distribution();
            
            
            act_pop_distr_dec_sum = sum_fields_of_struct(pop_distr_dec);
            exp_pop_distr_dec_sum = 1.0;
            
            testCase.verifyEqual(act_pop_distr_dec_sum,exp_pop_distr_dec_sum,...
                'AbsTol', 0.0001);
            
            
            act_pop_distr_size_sum = sum_fields_of_struct(pop_distr_size);
            exp_pop_distr_size_sum = size_of_pop;
            
            testCase.verifyEqual(act_pop_distr_size_sum,exp_pop_distr_size_sum,...
                'AbsTol', 0.0001);
            
            
            
        end
        
        function test_ceate_first_population(testCase)
            % Should return a matrix of random all elements must be
            % greater or equall to zero
            
            
            size_of_pop = 100;
            init_vals = [0.8; 0.8; 0.1; 0.1; 0.4; 0.4; 0.4];
            l_bounds = zeros(length(init_vals),1);
            u_bounds = ones(length(init_vals),1);
            
            fm_obj = GeneticAlgoEngine(size_of_pop,init_vals,...
                l_bounds,u_bounds);
            
            par_len = length(init_vals);
            
            act_result1 =...
                fm_obj.ceate_first_population();
            
            act_result2 =...
                fm_obj.ceate_first_population();
            
            exp_result = zeros(size_of_pop, par_len);
            
            % test size
            testCase.verifyEqual(size(act_result1),size(exp_result));
            % test randomness
            testCase.verifyNotEqual(act_result1, act_result2);
            % test upper & lower bounds
            for k = 1:size_of_pop
                pop = fm_obj.ceate_first_population();
                for p = 1:par_len
                    testCase.verifyTrue(all(pop(:, p) > l_bounds(p)));
                    testCase.verifyTrue(all(pop(:, p) < u_bounds(p)));
                end
            end
            
            init_vals = [0.8; 0.8; 0.1;  0.1;  0.4; 0.4; 0.4];
            l_bounds =  [0.5;   -1;  0.05; 0.01; 0;   0;   0.01];
            u_bounds =  [1;   1;   0.5;  0.2;    0.6;   1;   6];
            
            fm_obj2 = GeneticAlgoEngine(size_of_pop,init_vals,...
                l_bounds,u_bounds);
            
            for k = 1:size_of_pop
                pop = fm_obj2.ceate_first_population();
                for p = 1:par_len
                    testCase.verifyTrue(all(pop(:, p) > l_bounds(p)));
                    testCase.verifyTrue(all(pop(:, p) < u_bounds(p)));
                end
            end
            
            
            
        end
        
        function test_get_new_population(testCase)
            % Should return a parameter matrix of new population
            % best are kept and rest is randomize, all elements must be
            % greater or equall to zero
            
                       
            size_of_pop = 100;
            init_vals = [0.8; 0.8; 0.1;  0.1;  0.4; 0.4; 0.4];
            l_bounds =  [0.5;   -1;  0.05; 0.01; 0;   0;   0.01];
            u_bounds =  [1;   1;   0.5;  0.2;    0.6;   1;   6];
            
            fm_obj = GeneticAlgoEngine(size_of_pop,init_vals,...
                l_bounds,u_bounds);
            
            
            par_len = length(init_vals);
            
            
            pop_rnd_par = fm_obj.ceate_first_population();
            pop_errors = rand(size_of_pop,1); % random 

     
            act_result1 =...
                fm_obj.get_new_population(pop_errors,...
                pop_rnd_par, 1);
                       
            act_result2 =...
                fm_obj.get_new_population(pop_errors,...
                pop_rnd_par, 2);
            
            [~, indx] = sort(pop_errors);
            act_result_first_10 = act_result1(1:10,:);
            exp_result_first_10 = pop_rnd_par(indx(1:10),:);
            
            testCase.verifyEqual(act_result_first_10,exp_result_first_10);
            testCase.verifyNotEqual(act_result1,act_result2);
            
            for k = 1:size_of_pop*10
                for p = 1:par_len
                    testCase.verifyTrue(all(act_result1(:, p) > l_bounds(p)));
                    testCase.verifyTrue(all(act_result1(:, p) < u_bounds(p)));
                end
            end
            
        end
        
    end
end


