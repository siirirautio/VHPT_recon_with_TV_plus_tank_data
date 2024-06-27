% Assume given the coefficients of a complex-valued function f in the basis
%
%   phi_0(t) = (2pi)^{-1/2},
%   phi_1(t) = pi^{-1/2} cos t,
%   phi_2(t) = pi^{-1/2} sin t,
%   phi_3(t) = pi^{-1/2} cos 2t,
%   phi_4(t) = pi^{-1/2} sin 2t,
%   ...
%   phi_{2N-1}(t) = pi^{-1/2} cos Ntrig*t,
%   phi_{2N}(t)   = pi^{-1/2} sin Ntrig*t.
%
% The coefficient vector has 2*(2*Ntrig+1) elements.
% The first 2*Ntrig+1 elements are the (real-valued) coefficients
% of the real part of f, and the remaining 2*Ntrig+1 elements 
% are the coefficients of the imaginary part of f.
%
% This function returns the (complex) values of f on the unit circle 
% at the evaluation points 2*pi*[1:Nx]/Nx.
%
% Arguments:
% coef    Vertical vector of length 2*(2*Ntrig+1)
% Ntrig   Order of trigonometric approximation
% Nx      Number of evaluation points
%
% Returns:
% f       Point values of f, vertical vector of length Nx
% fii     Evaluation angles (in radians)
%
% Samuli Siltanen May 2009

function [f,fii] = coef2fun(coef,Ntrig,Nx)

% Check that the argument function coef has correct size
if ~((size(coef,1)==2*(2*Ntrig+1)) && (size(coef,2)==1))
    error('Argument vector f has wrong size')
end

% Create evaluation points
fii  = 2*pi*[1:Nx]/Nx;
fii  = fii(:);
Dfii = fii(2)-fii(1);

% Reconstruct real part approximately from the coefficients
f_re = zeros(size(fii));
% Constant component
f_re = f_re + 1/sqrt(2*pi)*coef(1);
% Frequency components
for iii = 1:Ntrig
    f_re = f_re + 1/sqrt(pi)*coef(2*iii)*cos(iii*fii);
    f_re = f_re + 1/sqrt(pi)*coef(2*iii+1)*sin(iii*fii);
end

% Reconstruct imaginary part approximately from the coefficients
f_im = zeros(size(fii));
% Constant component
f_im = f_im + 1/sqrt(2*pi)*coef(2*Ntrig+1+1);
% Frequency components
for iii = 1:Ntrig
    f_im = f_im + 1/sqrt(pi)*coef(2*Ntrig+1+2*iii)*cos(iii*fii);
    f_im = f_im + 1/sqrt(pi)*coef(2*Ntrig+1+2*iii+1)*sin(iii*fii);
end

% Combine the results
f = f_re + i*f_im;
