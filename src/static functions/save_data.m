function file_path = save_data(save_folder, measured_data, filename)
%Save measured results
%   save folder - save folder destination
%   measured_data - data to save

% check if an argument exists
if ~exist('filename','var')
      filename = strrep(datestr(now), '_', '');
end

file_path = [save_folder, '/', filename];
save(file_path, 'measured_data');
end

