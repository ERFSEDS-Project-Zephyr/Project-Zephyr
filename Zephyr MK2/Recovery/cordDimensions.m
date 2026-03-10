function [foldCount,foldLength] = cordDimensions(cordLength,cordThickness)

% Note: Assumes internal diameter of rocket is 2.75 in
% Only inches work 

foldCount = 2.75/cordThickness ;
foldLength = (cordLength-3) / foldCount ;

fprintf('Fold Count=%0.2f ; Fold Length=%0.2f \n',foldCount,foldLength)