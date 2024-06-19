% Compute the matrix representation of the real-linear operator P_0
% for compex-valued functions.
%
% Argument:
% Ntrig    order of trigonometric approximation
%
% Returns:
% res      matrix of operator, size (2*(2*Ntrig+1))x(2*(2*Ntrig+1))
%
% Samuli Siltanen May 2009

function res = P0_matrix(Ntrig,H0)

% Construct complex averaging operator
aver = Aver_matrix(Ntrig);

% Construct complex identity operator
iden = eye(2*(2*Ntrig+1));

% Construct Hmu
%load  data/Hmu H0
HH = [[H0,zeros(2*Ntrig+1)];[zeros(2*Ntrig+1),H0]];

% Construct multiplication operator i*
imult = [[zeros(2*Ntrig+1),-eye(2*Ntrig+1)];[eye(2*Ntrig+1),zeros(2*Ntrig+1)]];

% Put all of the above together
res = .5*(iden + imult*HH) + .5*aver;
