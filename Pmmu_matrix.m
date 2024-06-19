% Compute the matrix representation of the real-linear operator P_{-mu}
% for compex-valued functions.
%
% Argument:
% Ntrig    order of trigonometric approximation
% Hmu      Matrix of (mu)-Hilbert transform Hmu in the trig basis
% Hmmu     Matrix of (-mu)-Hilbert transform Hmmu in the trig basis
%
% Returns:
% res      matrix of operator, size (2*(2*Ntrig+1))x(2*(2*Ntrig+1))
%
% Samuli Siltanen May 2009


function res = Pmmu_matrix(Ntrig,Hmu,Hmmu)

% Construct complex averaging operator
aver = Aver_matrix(Ntrig);

% Construct complex identity operator
iden = eye(2*(2*Ntrig+1));

% Construct Hilbert transform matrix
HH = [[Hmmu,zeros(2*Ntrig+1)];[zeros(2*Ntrig+1),Hmu]];

% Construct multiplication operator i*
imult = [[zeros(2*Ntrig+1),-eye(2*Ntrig+1)];[eye(2*Ntrig+1),zeros(2*Ntrig+1)]];

% Put all of the above together
res = .5*(iden + imult*HH) + .5*aver;
