clear; clc; close all;

dart = readcell('ZephyrMk2Flight1_Dart.csv');
booster = readcell('Zephyrmk2_flight1_booster.csv');

dart_time = cell2mat(dart(2:end,5)); T = find(dart_time==0,1,'last');
dart_time = dart_time(T:end);
dart_height = cell2mat(dart((T+1):end,11));
dart_accelx = cell2mat(dart((T+1):end,17));
dart_accely = cell2mat(dart((T+1):end,18));
dart_accelz = cell2mat(dart((T+1):end,19));
dart_accel = sign(dart_accelx).*sqrt(dart_accelx.^2 + dart_accely.^2 + dart_accelz.^2);
dart_veloc = 0;
for i = 2:length(dart_time)
    dart_veloc(i,:) = dart_veloc(i-1) + 0.5*(dart_accel(i) + dart_accel(i-1))*(dart_time(i) - dart_time(i-1));
end

plot(dart_time,dart_accel);
xlim([0 6]);