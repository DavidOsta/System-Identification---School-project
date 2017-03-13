function [ path ] = get_path_of_data_folder(class_type)
%GET_PATH_OF_MEASURED_DATA_FOLDER Static function returns path for
%measurement data
%
switch class_type
    case 'ident'
        path = 'data/identified_models';
    case 'measurement' %{'pie','pie3'}
        path = 'data/measured_data';
    otherwise
        path = 'data';
end
end

