% This was originally a function named obliquerelations downloaded from the
% MATLAB file exchange. The original credits are 
% David Padgett (2023). Oblique Shock Relations Solver 
% (https://www.mathworks.com/matlabcentral/fileexchange/
% 28242-oblique-shock-relations-solver), MATLAB Central File Exchange. 
% Retrieved November 9, 2023.
% The comments were subsequently modified for AE 319
%
% if delta and mach number are given then the format is
% theta=flowoblique('mach',M,'delta',delta,gamma);
% if theta and mach number are given then the format is
% theta=flowoblique('mach',M,'theta',theta,gamma);
% if theta and delta are given then the format is
% mach=flowoblique('delta',delta,'theta',theta,gamma);
%
% NOTE -  ALL ANGLES ARE IN RADIANS
%
% Function to solve the mach, theta, delta oblique shock relations for
% supersonic flow incident on a wedge.
%
% Currently works for solving weak shock angles
%
% Mach - the mach number of the flow (dimensionless)
% theta - the shock angle (radians)
% delta - the wedge half-angle (in radians).  In some cases called the turning angle of
%         the flow, but this code is not intended to solve Prandtle-Meyer expansion
%         problems.
%
% See, e.g., http://en.wikipedia.org/wiki/Oblique_shock for a description
% of the geometric relationship between theta and delta.
%
% Given a mach number and theta, delta is solved using 
% 
%  tan(delta) = 2cot(theta)*(M^2sin^2(theta) - 1)/(M^2(gamma + cos(2theta) + 1)
%
% Given a delta and a mach number, the above equation is solved numerically
% using Newton's method. 
%
% Given a delta and a theta, the above equation is solved for M using
% Newton's method.
%
% The two input vectors must be of equal length.  If they are not, the
% shorter input vector will be lengthened and the last value of the shorter
% vector will be repeated to expand the vector.

function retval = flowoblique(label1, in1, label2, in2, gamma)
retval = -1;

% Expand the shortest vector to be the same length as the longest vector.
% Hold the last value of the shorter vector.
if length(in1) < length(in2)
    newin1 = zeros(size(in2));
    newin1(1:length(in1)) = in1;
    newin1(length(in1):length(newin1)) = ones(size(in2))*in1(length(in1));
    in1 = newin1;
elseif length(in1) > length(in2)
    newin2 = zeros(size(in1));
    newin2(1:length(in2)) = in2;
    newin2(length(in2):length(newin2)) = ones(size(in1))*in2(length(in2));
    in2 = newin2;
end


if length(in1) ~= length(in2)
    error('Inputs must be of the same length');
else

    % Given a mach and a theta, solve for delta directly
    if strcmpi(label1, 'mach') && strcmpi(label2, 'theta')
        retval = atan2(2*cot(in2).*(in1.^2.*(sin(in2)).^2.-1), in1.^2.*(gamma + cos(2*in2))+2);
    elseif strcmpi(label2, 'mach') && strcmpi(label1, 'theta')
        retval = atan2(2*cot(in1).*(in2.^2.*(sin(in1)).^2.-1), in1.^2.*(gamma + cos(2*in1))+2);
        
    % Given a delta and a theta, solve for M using Newton's method
    elseif strcmpi(label1, 'theta') && strcmpi(label2, 'delta')
        b = in1;
        delta = in2;
        f = @(M) (2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2)-tan(delta);
        fp = @(M) -1*(4*M.*cot(b).*(M.^2.*sin(b).^2-1).*(cos(2*b)+gamma))./((M.^2.*(cos(2*b)+gamma)+2).^2)...                
                + (4*M.*cos(b).*sin(b))./(M.^2.*(cos(2*b)+gamma) +2);
            
        xold = ones(size(in1))*1.5;
        xnew = xold - f(xold)./fp(xold);
        while abs(xold - xnew) > 10^-8*ones(size(in1))
            
            xold = xnew;
            xnew = xold - f(xold)./fp(xold);
        end
        retval = xnew;
    elseif strcmpi(label1, 'delta') && strcmpi(label2, 'theta')
        b = in2;
        delta = in1;
        f = @(M) (2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2)-tan(delta);
        fp = @(M) -1*(4*M.*cot(b).*(M.^2.*sin(b).^2-1).*(cos(2*b)+gamma))./((M.^2.*(cos(2*b)+gamma)+2).^2)...                
                + (4*M.*cos(b).*sin(b))./(M.^2.*(cos(2*b)+gamma) +2);
             
        xold = ones(size(in1));
        xnew = xold - f(xold)./fp(xold);
        while abs(xold - xnew) > 10^-8*ones(size(in1))
            
            xold = xnew;
            xnew = xold - f(xold)./fp(xold);
        end
        retval = xnew;
        
    % Given a mach and a delta, solve for theta using Newton's method
    % (return NaN for out of bounds values)
    elseif strcmpi(label2, 'mach') && strcmpi(label1, 'delta')
        M = in2;
        delta = in1;
        f = @(b) (2*cot(b).*(M.^2.*(sin(b)).^2-1))/(M.^2.*(gamma + cos(2*b))+2)-tan(delta);
        fp = @(b) (4*M.^2.*sin(2*b).*cot(b).*(M.^2.*sin(b).^2-1))./((M.^2.*(cos(2*b)+gamma)+2).^2)...                
                + (4*M.^2.*cos(b).^2 - 2*csc(b).^2.*(M.^2.*(sin(b)).^2-1))./(M.^2.*(cos(2*b)+gamma) +2);
        overmax = ones(size(delta));
        for i = 1:length(in1)
           M = in1(i);
           f1 = @(b) -1*(2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2);
           x = fminbnd(f1,0,pi/2);
           if(delta(i) >= atan(-1*f1(x)))
               overmax(i) = NaN;
           end
        end  
        xold = ones(size(in1))*0.1;
        xnew = xold - f(xold)/fp(xold);
        while abs(xold - xnew) > 10^-8*ones(size(in1))
            
            xold = xnew;
            xnew = xold - f(xold)/fp(xold);
        end
        retval = xnew;
        retval = retval.*overmax;
    elseif strcmpi(label1, 'mach') && strcmpi(label2, 'delta')
        M = in1;
        delta = in2;
        f = @(b) (2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2)-tan(delta);
        fp = @(b) (4*M.^2.*sin(2*b).*cot(b).*(M.^2.*sin(b).^2-1))./((M.^2.*(cos(2*b)+gamma)+2).^2)...                
                + (4*M.^2.*cos(b).^2 - 2*csc(b).^2.*(M.^2.*(sin(b)).^2-1))./(M.^2.*(cos(2*b)+gamma) +2);
        overmax = ones(size(delta));
        for i = 1:length(in1)
           M = in1(i);
           f1 = @(b) -1*(2*cot(b).*(M.^2.*(sin(b)).^2-1))./(M.^2.*(gamma + cos(2*b))+2);
           x = fminbnd(f1,0,pi/2);
           if(delta(i) >= atan(-1*f1(x)))
               overmax(i) = NaN;
           end
        end

        xold = ones(size(in1))*0.1;
        xnew = xold - f(xold)./fp(xold);
        while abs(xold - xnew) > 10^-8*ones(size(in1))
            xold = xnew;
            xnew = xold - f(xold)./fp(xold);
        end
        retval = xnew;
        retval = retval.*overmax;
    end
end
