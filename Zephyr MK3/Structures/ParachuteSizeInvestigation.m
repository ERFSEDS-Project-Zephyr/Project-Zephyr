clear; close all; clc; warning('off'); addpath('../..');
set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',24,'defaultAxesFontName','Times New Roman');

% Givens
CK = 1.3; % Opening Force Coefficient
rho = 1.225; % Air Density (kg m^-3)
g = 9.81; % Gravity (m s^-2)
y0 = 712; % Height at line stretch (m)
v0 = 362; % Velocity at line stretch (m s^-1)
m_ball = 0; % Ballast Mass (kg)
m_dart = 3.767 + m_ball; % Dart Mass (kg)
m_nose = 0.088 + m_ball; % Nose Mass (kg)

% Simulation Parameters
force_nose_max = 6*convert(75,'lbf','N'); % Max. Nose Cone Force (N)
apogee_max = convert(14500,'ft','m'); % Max Apogee Height (m)
diameter_min = 0*0.0254; % Min. Diameter (m)
diameter_max = 60*0.0254; % Max. Diameter (m)
coeffDrag_min = 0.4; % Min. Drag Coefficient
coeffDrag_max = 0.6; % Max. Drag Coefficient
n_diameter = 1e5; % Resolution of Diameter
n_coeffDrag = 1e2; % Resolution of Drag Coeffcient

[D,CD] = meshgrid(linspace(diameter_min,diameter_max,n_diameter),...
                  linspace(coeffDrag_min,coeffDrag_max,n_coeffDrag));
CD_D_range = zeros(n_coeffDrag,3);

AD = (0.25*pi*D.^2).*CD;
vT = sqrt(2*m_dart*g./(rho*AD));
tf = 0;
% tf = 0.65*0.5*D/v0;
% alpha = (1+CK)*g*tf./(vT.^2);
vf = v0;
% vf = (sqrt(1 + 2*(v0 - g*tf).*alpha) - 1)./alpha;
force_nose = 0.5*rho*(vf.^2)*CK.*AD*(m_nose/m_dart);
apogee = y0 + v0*tf - 0.5*g*(tf.^2).*(1 + (1/3)*(1+2*CK)*((vf./vT).^2)) + ((vT.^2)/g).*log(sqrt(1 + (vf./vT).^2));

S_fmax = (1+sign(force_nose_max - force_nose))/2;
S_apogee = (1+sign(apogee_max - apogee))/2;
S = S_fmax.*S_apogee;
for k = 1:n_coeffDrag
    Imin = find(S(k,:)==1,1,'first'); Imax = find(S(k,:)==1,1,'last');
    CD_D_range(k,:) = [CD(k,1), D(1,Imin)/0.0254, D(1,Imax)/0.0254];
end
save('ParachuteSizeInvestigation.mat','CD_D_range');

figure(1);
subplot(1,2,1); hold on;
surf(D/0.0254,CD,S_fmax,255*cat(3,1-S,S_fmax,S_fmax*0)); shading interp;
hold off; view(0,90); pbaspect([1 1 1]); xticks(0:6:(diameter_max/0.0254));
xlim([diameter_min diameter_max]/0.0254); ylim([coeffDrag_min coeffDrag_max]);
title('Nose Cone Force Bifurcation'); xlabel('Parachute Diameter (inches)'); ylabel('Drag Coefficient');

subplot(1,2,2); hold on;
surf(D/0.0254,CD,S_apogee,255*cat(3,1-S,S_apogee,S_apogee*0)); shading interp; 
hold off; view(0,90); pbaspect([1 1 1]); xticks(0:6:(diameter_max/0.0254));
xlim([diameter_min diameter_max]/0.0254); ylim([coeffDrag_min coeffDrag_max]);
title('Apogee Height Bifurcation'); xlabel('Parachute Diameter (inches)'); ylabel('Drag Coefficient');
saveas(gca,'ParachuteSizeInvestigation','png');