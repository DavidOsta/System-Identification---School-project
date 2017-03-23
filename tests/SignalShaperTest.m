classdef SignalShaperTest < matlab.unittest.TestCase

    properties (TestParameter)

    end


    methods (Test)
        function test_constructor(testCase)
            % test initialized parameters

            [shaper_obj, sys_obj] = testCase.create_shaper_obj();

            [act_damped_natural_freq, act_beta, act_complex_poles] =...
                shaper_obj.get_init_parameters();

            [Wn,zeta] = damp(sys_obj);
            natural_freq = Wn(1);
            damping_ratio = zeta(1);
            
            exp_beta = natural_freq * damping_ratio;
            exp_damped_natural_freq = natural_freq * sqrt(1 - damping_ratio^2);
            pole1 = complex(-act_beta, (act_damped_natural_freq));
            pole2 = complex(-act_beta, -(act_damped_natural_freq));
            exp_complex_poles = [pole1; pole2];

            testCase.verifyEqual(act_damped_natural_freq,exp_damped_natural_freq,'AbsTol', 0.0001);
            testCase.verifyEqual(act_beta,exp_beta, 'AbsTol', 0.0001);
            testCase.verifyEqual(act_complex_poles,exp_complex_poles, 'AbsTol', 0.0001);


            testCase.verifyError(@()SignalShaper(1),...
                'wca:WrongConstructorsArg');
            testCase.verifyError(@()SignalShaper(1,3,'faf'),...
                'wca:WrongConstructorsArg');
            testCase.verifyError(@()SignalShaper(1,2,'hy'),...
                'wca:WrongConstructorsArg');
            testCase.verifyError(@()SignalShaper(),...
                'cam:ConstructorArgMiss');
        end


        function test_implementation_of_feedforward_ZV_shaper(testCase)
        % compare implementation of ZV shaper in state space form
        % with ZV shaper implemented in Simulink time domain

        [shaper_obj, sys_obj] = testCase.create_shaper_obj();

        [ZV_shaper_ss_model, shaper_par] = shaper_obj.get_ZV_shaper();
                
        [control_force, system_response] = ...
                testCase.run_sim(ZV_shaper_ss_model,...
                'tests_feedforward_shapers',...
                'ZV',...
                shaper_par,...
                sys_obj);

       testCase.verifyEqual(system_response.Data(:,1), system_response.Data(:,2),...
           'AbsTol', 0.01);

       testCase.verifyEqual(control_force.Data(:,1), control_force.Data(:,2),...
           'AbsTol', 0.01);

        end
        
        function test_implementation_of_feedforward_DZV_shaper(testCase)
        [shaper_obj, sys_obj] = testCase.create_shaper_obj();

        [DZV_shaper_ss_model, shaper_par] = shaper_obj.get_DZV_shaper();
                
        [control_force, system_response] = ...
                testCase.run_sim(DZV_shaper_ss_model,...
                'tests_feedforward_shapers',...
                'DZV',...
                shaper_par,...
                sys_obj);

       testCase.verifyEqual(system_response.Data(:,1), system_response.Data(:,2),...
           'AbsTol', 0.01);

       testCase.verifyEqual(control_force.Data(:,1), control_force.Data(:,2),...
           'AbsTol', 0.01);
        end


        function test_implementation_of_inverse_DZV_shaper(testCase)
        % compare implementation of ZV shaper in state space form
        % with ZV shaper implemented in Simulink time domain

        [shaper_obj, sys_obj] = testCase.create_shaper_obj();

        [ZV_shaper_ss_model, shaper_par] = shaper_obj.get_DZV_shaper();
                
        [control_force, system_response] = ...
                testCase.run_sim(ZV_shaper_ss_model,...
                'tests_inverse_shapers',...
                'DZV',...
                shaper_par,...
                sys_obj);

       testCase.verifyEqual(system_response.Data(:,1), system_response.Data(:,2),...
           'AbsTol', 0.01);

       testCase.verifyEqual(control_force.Data(:,1), control_force.Data(:,2),...
           'AbsTol', 0.01);

        end


        
        
    end
    methods(Access = private)

        function [shaper_obj, sys_obj] =...
                create_shaper_obj(testCase, tf_sys)
          %create_shaper_obj create shaper object for test cases
      
          if(exist('tf_sys', 'var'))
            sys_obj = tf_sys;
          else
            sys_obj= tf(1, [1 0.2 1]);
          end;
  
          [Wn,zeta] = damp(sys_obj);
          nat_freq = Wn(1);
          damping_rat = zeta(1);
          
          shaper_obj = SignalShaper(nat_freq, damping_rat);
          sys_parameters.natural_freq = nat_freq;
          sys_parameters.damping_ratio = damping_rat;
        end




        function [control_force, system_response] = ...
                run_sim(testCase, shaper_model, sim_model_name,...
                shaper_type, shaper_parameters, sys_obj)

            switch shaper_type
                case 'ZV'
                    shaper_switch = '1';
                case 'DZV'
                    shaper_switch = '2';
                case 'DZValpha'
                    shaper_switch = '3';
                otherwise
                    error('unknown type of shaper')
            end
            
            set_param([sim_model_name,'/shaper_switch'],...
                'Value', shaper_switch);
            set_param([sim_model_name,'/shaper_switch'],...
                'Value', shaper_switch);
            
            assignin('base', 'system_tf', sys_obj);
            assignin('base', 'shaper_model', shaper_model);
            assignin('base', 'a', shaper_parameters.gain_a);
            assignin('base', 'tau_add', shaper_parameters.tau_add);
            assignin('base', 'tau_sub', shaper_parameters.tau_sub);
            assignin('base', 'alpha', shaper_parameters.alpha);
            
            warning('off','all');
            simout = sim(sim_model_name, 'ReturnWorkspaceOutputs', 'on');
            
            control_force = simout.get('control_force');
            system_response = simout.get('system_response');

        end




    end



%         function test_init_of_attributes(testCase)
%             % test of class atributes
%             function struct_sum = sum_fields_of_struct(strct)
%                 % Helper function
%                 f_names = fieldnames(strct);
%                 f_len = length(f_names);
%                 f_vals = zeros(f_len,1);
%                 for k = 1:f_len
%                     f_vals(k) = strct.(f_names{k});
%                 end
%                 struct_sum = sum(f_vals);
%             end
%
%
%
%             size_of_pop = 1000;
%             init_vals = [0.8; 0.8; 0.1; 0.1; 0.4; 0.4; 0.4];
%             l_bounds = zeros(length(init_vals),1);
%             u_bounds = ones(length(init_vals),1);
%
%             fm_obj = GeneticAlgoEngine(size_of_pop,init_vals,...
%                 l_bounds,u_bounds);
%
%             [pop_distr_size, pop_distr_dec] =...
%                 fm_obj.get_pop_distribution();
%
%
%             act_pop_distr_dec_sum = sum_fields_of_struct(pop_distr_dec);
%             exp_pop_distr_dec_sum = 1.0;
%
%             testCase.verifyEqual(act_pop_distr_dec_sum,exp_pop_distr_dec_sum,...
%                 'AbsTol', 0.0001);
%
%
%             act_pop_distr_size_sum = sum_fields_of_struct(pop_distr_size);
%             exp_pop_distr_size_sum = size_of_pop;
%
%             testCase.verifyEqual(act_pop_distr_size_sum,exp_pop_distr_size_sum,...
%                 'AbsTol', 0.0001);
%
%
%
%         end
%
%         function test_ceate_first_population(testCase)
%             % Should return a matrix of random all elements must be
%             % greater or equall to zero
%
%
%             size_of_pop = 100;
%             init_vals = [0.8; 0.8; 0.1; 0.1; 0.4; 0.4; 0.4];
%             l_bounds = zeros(length(init_vals),1);
%             u_bounds = ones(length(init_vals),1);
%
%             fm_obj = GeneticAlgoEngine(size_of_pop,init_vals,...
%                 l_bounds,u_bounds);
%
%             par_len = length(init_vals);
%
%             act_result1 =...
%                 fm_obj.ceate_first_population();
%
%             act_result2 =...
%                 fm_obj.ceate_first_population();
%
%             exp_result = zeros(size_of_pop, par_len);
%
%             % test size
%             testCase.verifyEqual(size(act_result1),size(exp_result));
%             % test randomness
%             testCase.verifyNotEqual(act_result1, act_result2);
%             % test upper & lower bounds
%             for k = 1:size_of_pop
%                 pop = fm_obj.ceate_first_population();
%                 for p = 1:par_len
%                     testCase.verifyTrue(all(pop(:, p) > l_bounds(p)));
%                     testCase.verifyTrue(all(pop(:, p) < u_bounds(p)));
%                 end
%             end
%
%             init_vals = [0.8; 0.8; 0.1;  0.1;  0.4; 0.4; 0.4];
%             l_bounds =  [0.5;   -1;  0.05; 0.01; 0;   0;   0.01];
%             u_bounds =  [1;   1;   0.5;  0.2;    0.6;   1;   6];
%
%             fm_obj2 = GeneticAlgoEngine(size_of_pop,init_vals,...
%                 l_bounds,u_bounds);
%
%             for k = 1:size_of_pop
%                 pop = fm_obj2.ceate_first_population();
%                 for p = 1:par_len
%                     testCase.verifyTrue(all(pop(:, p) > l_bounds(p)));
%                     testCase.verifyTrue(all(pop(:, p) < u_bounds(p)));
%                 end
%             end
%
%
%
%         end
%
%         function test_get_new_population(testCase)
%             % Should return a parameter matrix of new population
%             % best are kept and rest is randomize, all elements must be
%             % greater or equall to zero
%
%
%             size_of_pop = 100;
%             init_vals = [0.8; 0.8; 0.1;  0.1;  0.4; 0.4; 0.4];
%             l_bounds =  [0.5;   -1;  0.05; 0.01; 0;   0;   0.01];
%             u_bounds =  [1;   1;   0.5;  0.2;    0.6;   1;   6];
%
%             fm_obj = GeneticAlgoEngine(size_of_pop,init_vals,...
%                 l_bounds,u_bounds);
%
%
%             par_len = length(init_vals);
%
%
%             pop_rnd_par = fm_obj.ceate_first_population();
%             pop_errors = rand(size_of_pop,1); % random
%
%
%             act_result1 =...
%                 fm_obj.get_new_population(pop_errors,...
%                 pop_rnd_par, 1);
%
%             act_result2 =...
%                 fm_obj.get_new_population(pop_errors,...
%                 pop_rnd_par, 2);
%
%             [~, indx] = sort(pop_errors);
%             act_result_first_10 = act_result1(1:10,:);
%             exp_result_first_10 = pop_rnd_par(indx(1:10),:);
%
%             testCase.verifyEqual(act_result_first_10,exp_result_first_10);
%             testCase.verifyNotEqual(act_result1,act_result2);
%
%             for k = 1:size_of_pop*10
%                 for p = 1:par_len
%                     testCase.verifyTrue(all(act_result1(:, p) > l_bounds(p)));
%                     testCase.verifyTrue(all(act_result1(:, p) < u_bounds(p)));
%                 end
%             end
%
%         end

end
