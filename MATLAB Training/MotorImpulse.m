clear; close all; clc; warning('off');

MP_L22 = readcell('MotorProperties_L2200.csv');
MP_M37 = readcell('MotorProperties_M3700.csv');

MP_L22 = cell2mat(MP_L22(:,2:end));
MP_M37 = cell2mat(MP_M37(:,2:end));

tBurn_L22 = MP_L22(1,end);
tBurn_M37 = MP_M37(1,end);

Tavg_L22 = mean(MP_L22(2,:));
Tavg_M37 = mean(MP_M37(2,:));

Tmax_L22 = max(MP_L22(2,:));
Tmax_M37 = max(MP_M37(2,:));

JT_L22 = trapz(MP_L22(1,:),MP_L22(2,:));
JT_M37 = trapz(MP_M37(1,:),MP_M37(2,:));

fprintf('Aerotech L2200G Properties:\n');
fprintf('Burn Time: %.2f seconds\n',tBurn_L22);
fprintf('Average Thrust: %.2f N\n',Tavg_L22);
fprintf('Maximum Thrust: %.2f N\n',Tmax_L22);
fprintf('Total Impulse: %.1f Ns\n',JT_L22);