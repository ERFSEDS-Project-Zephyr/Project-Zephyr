function [noFold,halfFold,zFold] = parachuteDimensions(parachuteRadius)

% Note: Untested for Experimental Parachute
% Any units work, inches recommended 

noFold = parachuteRadius ;
halfFold = parachuteRadius/2 ;
zFold = parachuteRadius/3 ; 

fprintf('No Fold=%0.2f ; Half Fold=%0.2f ; Z Fold=%0.2f \n',noFold,halfFold,zFold)

