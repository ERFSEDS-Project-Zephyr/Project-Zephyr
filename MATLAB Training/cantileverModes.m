function [l,N,V] = cantileverModes(n,res)
% cantileverModes: Determines the natural frequencies and natural modes of
% a uniform cantilever beam
%       l = cantileverModes(n) returns the nth natural frequency
%
%   [l,N] = cantileverModes(n) also returns the zeros of the natural mode shape.
%
% [l,N,V] = cantileverModes(n) also returns the displacements of the mode 
%           shape, scaled to have a maximum displacement of 1 at the far end.
%
% [l,N,V] = cantileverModes(__,res) allows for finer control over the resolution 
%           of the x-axis of the mode shape.
%
% Returns accurate mode shapes for n < 13
% This function requires the Optimization Toolbox.  

if nargin<2
    res = 100;
end
options = optimset('TolFun',1e-13,'TolX',1e-13,'Display','off');

l = fsolve(@(x) 1+cos(x).*cosh(x),(n - 0.5)*pi,optimset('Display','off'));
kappa = (cosh(l)+cos(l))/(sinh(l)+sin(l));
if n > 1
    N = fsolve(@(x) (cosh(l*x)-cos(l*x))-kappa*(sinh(l*x)-sin(l*x)),...
               ((1:(n-1))*pi + atan(1/kappa))/l,options);
    N = [0,N];
else
    N = 0;
end

x = linspace(0,1,res);
V = ((cosh(x*l) - cos(x*l))-kappa*(sinh(x*l) - sin(x*l)))/...
     ((cosh(l) - cos(l))-kappa*(sinh(l) - sin(l)));