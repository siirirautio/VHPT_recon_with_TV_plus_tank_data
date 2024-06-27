% Compute the matrix representation of the averaging operator Aver.
% We construct the Aver operator for real-valued functions only;
% then we construct the complex version block by block.
%
% Argument:
% Ntrig    order of trigonometric approximation
%
% Returns:
% res      matrix of averaging operator, size (2*2*Ntrig+1)x(2*2*Ntrig+1)
%
% Samuli Siltanen May 2009


function res = Aver_matrix(Ntrig)

% Construct the averaging operator matrix for real-valued functions 
Aver      = zeros(2*Ntrig+1);
Aver(1,1) = 1;

% Build the complex-valued version block by block
res = [[Aver,zeros(size(Aver))];[zeros(size(Aver)),Aver]];
