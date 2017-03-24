function file_path = save_data(measured_data,folder_path, filename)
%Save measured results
%   save folder - save folder destination
%   measured_data - data to save

% check if an argument exists
time_stamp = strrep(datestr(now), ':', '_');
if ~exist('filename','var')
    filename = time_stamp;
else
    filename = [filename, '-', time_stamp];
end

if ~exist('folder_path','var')
    folder_path = 'data';
end


file_path = [folder_path, '/', filename];
save(file_path, 'measured_data');
end
