function setup_statespace_block( model_name, block_name, struct_parameters )
%SETUP_LTI_BLOCK setup simulink's state space block
%   Detailed explanation goes here

try
    find_system(model_name);
catch % model is not opened
    open_system(model_name);
end

try
    find_system([model_name,'/', block_name]);
    
    set_param([model_name, '/', block_name],...
    'A', mat2str(struct_parameters.A),...
    'B', mat2str(struct_parameters.B),...
    'C', mat2str(struct_parameters.C),...
    'D', mat2str(struct_parameters.D),...
    'X0', mat2str(struct_parameters.X));

catch ME
    switch ME.identifier
        case {'Simulink:Commands:FindSystemValidBlockDiagramOrParameter',...
                'Simulink:Commands:FindSystemNoBlock'}
            fprintf(['\n\tBlock "', block_name, '" does not exist']);
        otherwise
            rethrow(ME)
    end
end


end

