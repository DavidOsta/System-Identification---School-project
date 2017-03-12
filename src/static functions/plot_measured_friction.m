function plot_measured_friction( file_name )
% PLOT_MEASURED_FRICTION Plot results from friction measurement
% file_path - filename of measured data


% plot_friction_ident
obj = FrictionMeasurement();
obj.plot_results(file_name);
end

