function readThrustCurve(filename,motorName)
% Function takes in .txt file with RockSim data of motor and saves data as 
% .mat file in the local file folder

% Read data file
motorData = readcell(filename);
motorData = motorData(:,2:5);

% Eliminate unwanted characters from RockSim format
text2remove = ['"',"t","f","m","cg","/",">","="];
for k = 1:size(motorData,1)
    for j = 1:size(motorData,2)
        motorData(k,j) = erase(motorData(k,j),text2remove);
    end
end
motorData = str2double(motorData)';

% Convert matrix to struct and save
eval([motorName,'.time = motorData(1,:)']);
eval([motorName,'.thrust = motorData(2,:)']);
eval([motorName,'.mass = motorData(3,:)']);
eval([motorName,'.cg = motorData(4,:)']);
save([motorName,'_Properties.mat'],motorName);