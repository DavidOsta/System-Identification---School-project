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

        [DZV_shaper_ss_model, shaper_par] = shaper_obj.get_DZV_shaper_slow();
                
        [control_force, system_response] = ...
                testCase.run_sim(DZV_shaper_ss_model,...
                'tests_feedforward_shapers',...
                'DZVslow',...
                shaper_par,...
                sys_obj);

       testCase.verifyEqual(system_response.Data(:,1), system_response.Data(:,2),...
           'AbsTol', 0.01);

       testCase.verifyEqual(control_force.Data(:,1), control_force.Data(:,2),...
           'AbsTol', 0.01);
        end


        function test_implementation_of_inverse_DZV_shaper_slow(testCase)
        % compare implementation of ZV shaper in state space form
        % with ZV shaper implemented in Simulink time domain

        [shaper_obj, sys_obj] = testCase.create_shaper_obj();

        [ZV_shaper_ss_model, shaper_par] = shaper_obj.get_inverse_DZV_shaper_slow();
                
        [control_force, system_response] = ...
                testCase.run_sim(ZV_shaper_ss_model,...
                'tests_inverse_shapers',...
                'DZVslow',...
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

        [ZV_shaper_ss_model, shaper_par] = shaper_obj.get_inverse_DZV_shaper();
                
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
                case 'DZVslow'
                    shaper_switch = '2';
                case 'DZV'
                    shaper_switch = '3';
                case 'DZValpha'
                    shaper_switch = '4';
                otherwise
                    error('unknown type of shaper')
            end
            
            load_system(sim_model_name);

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
            bdclose(sim_model_name);
        end




    end




end
