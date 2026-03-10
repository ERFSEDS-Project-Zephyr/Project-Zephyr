function [Setup,BodyDims,FinDims,MassDistribution] = MasterFileReader(filename)
% % % % Function takes the name of an excel file specialized for Zephyr data and
% outputs MATLAB structures that hold data referenced by their names.
%   Detailed explanation goes here

Datasheet = readcell(filename);

%% Read Setup Section
SetupCell = string(Datasheet(3:end,1:4));
SetupCell(ismissing(SetupCell)) = [];
SetupCell = reshape(SetupCell,[],4);

for k = 1:size(SetupCell,1)
    eval(strcat("Setup.",SetupCell(k,3),"=str2double(SetupCell(k,4));"));
end

%% Read Body Dims Section
BodyCell = string(Datasheet(3:end,5:8));
BodyCell(ismissing(BodyCell)) = [];
BodyCell = reshape(BodyCell,[],4);

for k = 1:size(BodyCell,1)
    eval(strcat("BodyDims.",BodyCell(k,3),"=str2double(BodyCell(k,4));"));
end

%% Read Fin Dims Section
FinCell = string(Datasheet(3:end,9:(11+Setup.numStages)));
FinCell(ismissing(FinCell)) = [];
FinCell = reshape(FinCell,[],3+Setup.numStages);

for k = 1:size(FinCell,1)
    eval(strcat("FinDims.",FinCell(k,3),"=str2double(FinCell(k,4:(3+Setup.numStages)));"));
end

%% Read Mass Distribution
MassCell = string(Datasheet(2:end,(12+Setup.numStages):end));
MassCell(ismissing(MassCell)) = [];
MassCell = reshape(MassCell,[],4);

eval(strcat("MassDistribution.Component=MassCell(2:end,1);"));
for k = 2:size(MassCell,2)
    eval(strcat("MassDistribution.",extractBefore(MassCell(1,k)," "),"=str2double(MassCell(2:end,k));"));
end

%% Display to Command Window
fprintf('Setup \n');
fprintf('_____________________________________________________________________________________________\n');
for t = 1:size(SetupCell,1)
    if strcmpi(SetupCell(t,2),"1")
        fprintf('%s: %.2f \n',SetupCell(t,1),str2double(SetupCell(t,4)));
    else
        fprintf('%s(%s): %.2f \n',SetupCell(t,1),SetupCell(t,2),str2double(SetupCell(t,4)));
    end
end

fprintf('\nBody Dimensions \n')
fprintf('_____________________________________________________________________________________________\n');
for t = 1:size(BodyCell,1)
    if strcmpi(BodyCell(t,2),"1")
        fprintf('%s: %.2f \n',BodyCell(t,1),str2double(BodyCell(t,4)));
    else
        fprintf('%s(%s): %.2f \n',BodyCell(t,1),BodyCell(t,2),str2double(BodyCell(t,4)));
    end
end

fprintf('\nFin Dimensions \n')
fprintf('_____________________________________________________________________________________________\n');
for t = 1:size(FinCell,1)
    if strcmpi(FinCell(t,2),"1")
        fprintf('%s: [',FinCell(t,1));
        fprintf('%5.2f, ',str2double(FinCell(t,4:(3+Setup.numStages))));
        fprintf('\b\b]\n');
    else
        fprintf('%s(%s): [',FinCell(t,1),FinCell(t,2));
        fprintf('%5.2f, ',str2double(FinCell(t,4:(3+Setup.numStages))));
        fprintf('\b\b]\n');
    end
end

fprintf('\nPoint Mass Distributions \n')
fprintf('_____________________________________________________________________________________________\n');
fprintf('%25s\t%s\t%s\t%s\n',MassCell(1,1),MassCell(1,2),MassCell(1,3),MassCell(1,4));
for t = 2:(size(MassCell,1)-1)
    fprintf('%25s\t%.3f\t\t%.3f\t\t%.3f\n',MassCell(t,1),str2double(MassCell(t,2)),str2double(MassCell(t,3)),str2double(MassCell(t,4)));
end