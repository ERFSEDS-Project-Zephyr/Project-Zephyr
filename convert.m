function x = convert(x,u1,u2,AOD)
% Units: https://www.mathworks.com/help/symbolic/units-list.html

x = sym(x); u1 = str2symunit(u1); u2 = str2symunit(u2);

if ~exist('AOD','var')
    x = double(separateUnits(unitConvert(x*u1,u2)));
else
    if AOD==0
        x = double(separateUnits(unitConvert(x*u1,u2,'Temperature','Absolute')));
    else
        x = double(separateUnits(unitConvert(x*u1,u2,'Temperature','Difference')));
    end
end