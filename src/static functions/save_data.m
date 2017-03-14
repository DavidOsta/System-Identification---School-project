function file_path = save_data(measured_data,folder_path, filename)
%Save measured results
%   save folder - save folder destination
%   measured_data - data to save

% check if an argument exists
time_stamp = strrep(datestr(now), '_', '');
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


%         function save_data(self,name_arg)
%             %Save measured results
%             %   save folder - save folder destination
%             %   measured_data - data to save
%             
%             prefix = 'friction_measurement_' ;
%             if ~exist('name_arg','var')
%                 timestamp = strrep(datestr(now), ':', '_');
%                 filename = [prefix, timestamp];
%             else
%                 filename = [prefix, name_arg];
%             end
%             
%             file_path = [self.save_folder, '/', filename];
%             simout = self.simout;
%             save(file_path, 'simout');
%         end