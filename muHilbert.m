% Computes and saves the mu-Hilbert transforms H_\mu and H_{-\mu} and H0
% acting on _real_valued_ functions. Applying the operators to
% complex-valued functions needs an extra step which is implemented outside
% this file.
%
% Hmu is solved from the equation
%
%   TanDer(Hmu(g)) = DN(g).
%
% Furthermore, we compute the matrix Hmmu of H_{-mu} as follows:
% (i)   remove first row and first column from Hmu, call the result Htmp,
% (ii)  write Hmmu = -inv(Htmp),
% (iii) add zero row and column into the middle of Hmmu.
%
% Routine DN_comp.m must be computed before this file.
%
% Arguments:
% DN      Dirichlet-to-Neumann matrix of conductivity sigma in the trig basis
% DN1     Dirichlet-to-Neumann matrix of conductivity 1 in the trig basis
% Ntrig   order of trigonometric approximation
%
% Returns:
% H0      Matrix of (0)-Hilbert transform
% Hmu     Matrix of (mu)-Hilbert transform
% Hmmu    Matrix of (-mu)-Hilbert transform
%
% Samuli Siltanen April 2019
% MODIFIED: the WHY comments removed and the sign problem fixed by
%           correcting the TD matrix. (Allan Perämäki)

function [H0,Hmu,Hmmu] = muHilbert(DN,DN1,Ntrig)

% Construct tangential derivative matrix analytically in the subspace of
% functions with mean value zero
TD = zeros(2*Ntrig,2*Ntrig);
for iii = 1:Ntrig
    TD(2*(iii-1)+1,2*(iii-1)+2) = iii;
    TD(2*(iii-1)+2,2*(iii-1)+1) = -iii;
end

% Construct a matrix that acts as inverse of TD in the subspace of
% functions with mean value zero
TDinv = inv(TD);

% Add appropriate zero row and zero column to complete the computation.
% Note that we choose the TDinv matrix to map constant functions to zero.
TD = [zeros(2*Ntrig,1),TD];
TD = [zeros(1,2*Ntrig+1);TD];
TDinv = [zeros(2*Ntrig,1),TDinv];
TDinv = [zeros(1,2*Ntrig+1);TDinv];

% Construct Hmu and H0 operators
Hmu  = TDinv*DN; 
H0   = TDinv*DN1; 

% Construct H_{-\mu} operator as matrix Hmmu
Htmp = Hmu(2:end,2:end);
Hmmu = -inv(Htmp);

% Add appropriate zero row and zero column to complete the computation
Hmmu = [zeros(2*Ntrig,1),Hmmu];
Hmmu = [zeros(1,2*Ntrig+1);Hmmu];


