function [ repaired_ts ] = repair_measured_position( measured_ts )
%REPAIR_DATA There was a overflow problem in our sensors
%   Use this to correct your data until sensors are fixed

data = measured_ts.Data;
max_data = max(data);
min_data = min(data);
offset = max_data - min_data;

for k = 1:length(data)
    if(data(k) < -max_data/2)
        data(k) = data(k) + offset;
    end;
end;

repaired_ts = timeseries(data, measured_ts.Time);
end

