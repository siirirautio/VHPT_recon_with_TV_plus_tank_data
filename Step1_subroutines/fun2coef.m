% Assume given a complex-valued function on the unit circle represented 
% by point values at the evaluation points 2*pi*[1:Nx]/Nx.
%
% This function computes the coefficients of f in the basis
%   phi_0(t) = (2pi)^{-1/2},
%   phi_1(t) = pi^{-1/2} cos t,
%   phi_2(t) = pi^{-1/2} sin t,
%   phi_3(t) = pi^{-1/2} cos 2t,
%   phi_4(t) = pi^{-1/2} sin 2t,
%   ...
%   phi_{2N-1}(t) = pi^{-1/2} cos Ntrig*t,
%   phi_{2N}(t)   = pi^{-1/2} sin Ntrig*t.
%
% The real and imaginary parts of f are treated separately: 
% the coefficient vector returned has 2*(2*Ntrig+1) elements.
% The first 2*Ntrig+1 elements are the (real-valued) coefficients
% of the real part of f, and the remaining 2*Ntrig+1 elements 
% are the coefficients of the imaginary part of f.
%
% Arguments:
% f       Point values of f, vertical vector of length Nx
% Nx      Number of evaluation points
% Ntrig   Order of trigonometric approximation
%
% Returns:
% coef    Vertical vector of length 2*(2*Ntrig+1)
%
% Samuli Siltanen April 2009

function coef = fun2coef(f,Nx,Ntrig)

% Check that the argument function f has correct size
if ~((size(f,1)==Nx) && (size(f,2)==1))
    error('Argument vector f has wrong size')
end

% Initialize the result
coef = zeros(2*(2*Ntrig+1),1);

% Create evaluation points
fii  = 2*pi*[1:Nx]/Nx;
fii  = fii(:);
Dfii = fii(2)-fii(1);

% Constant component of real part
coef(1) = 1/sqrt(2*pi)*Dfii*sum(real(f));
% Frequency components of real part
for iii = 1:Ntrig
    coef(2*iii)   = 1/sqrt(pi)*Dfii*sum(real(f).*cos(iii*fii));
    coef(2*iii+1) = 1/sqrt(pi)*Dfii*sum(real(f).*sin(iii*fii));
end

% Constant component of imaginary part
coef(2*Ntrig+1+1) = 1/sqrt(2*pi)*Dfii*sum(imag(f));
% Frequency components of imaginary part
for iii = 1:Ntrig
    coef(2*Ntrig+1+2*iii)   = 1/sqrt(pi)*Dfii*sum(imag(f).*cos(iii*fii));
    coef(2*Ntrig+1+2*iii+1) = 1/sqrt(pi)*Dfii*sum(imag(f).*sin(iii*fii));
end

