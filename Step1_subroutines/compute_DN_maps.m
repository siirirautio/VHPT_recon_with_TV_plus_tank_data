%===================================================================================================
%===================================================================================================
% This function computes the DN maps necessary for VHPT imaging from ACT5 EIT electrode data. It's
% assumed here that there are an even number L of electrodes.

% INPUTS:
% datadir = directory where data is stored, including forward or backslashes as appropriate
% V_clb_fname   = name of .mat file containing calibration voltages from ACT5
% V_trg_fname   = name of .mat file containing target voltages from ACT5
% V_1cem_fname  = name of .mat file containing voltages simulated on a homogeneous domain using CEM

% OUTPUTS:
% DN            = a symmetric, noise-robust DN map for input in the VHPT imaging chain
% DN1           = the DN map for the continuum model, corresponding to sigma = 1
% Ntrig         = the number of sines / cosines in the trigonometric current pattern basis

% Copyright (c) 2024 Melody Alsaker
% Date last modified: June 2024
%===================================================================================================
%===================================================================================================

function [DN, DN1, Ntrig] = compute_DN_maps(datadir, V_clb_fname, V_trg_fname, V_1cem_fname)


% Format of Current Patterns: Each column of the voltage matrices corresponds to a trigonometric
% current pattern, which initially use the convention of the ACT5 machine:
% [ cos(theta), cos(2*theta), ... ,cos(L*theta/2), sin(theta), sin(2*theta), ..., sin(L*theta/2) ];
% The sin(L*theta/2) column will be dropped because we need linear independence.
% The cos(L*theta/2) column will be dropped to ensure there are the same number of sines & cosines.
% We will then permute the columns to correspond to the ordering
% [ cos(theta), sin(theta), ... ,cos(15*theta), sin(15*theta)];


% Load and extract calibration and target voltage data, size LxL.
% Average the voltages over all frames of each dataset
V_clb_multi = load([datadir, V_clb_fname]);
V_clb = real(V_clb_multi.frame_voltage);
V_clb = mean(V_clb,3)';

V_trg_multi = load([datadir, V_trg_fname]);
V_trg = real(V_trg_multi.frame_voltage);
V_trg = mean(V_trg,3)';

% Current pattern matrix, size L x L. Each column gives one current pattern. This is read in
% along with the voltage data from ACT5.
L = length(V_trg_multi.cur_pattern);  % Number of electrodes
J = V_trg_multi.cur_pattern;          % Current pattern matrix

% Load (homogeneous) simulated CEM data, make it size LxL (just to match the others)
V_1cem_multi = load([datadir, V_1cem_fname]);
V_1cem = full(V_1cem_multi.voltages);
V_1cem = [V_1cem, zeros(L,1)];

Ntrig = (L-2) / 2;          % Number of sines / cosines in our basis

% Create column permutation matrix that will be applied to reorder current patterns
CPM = zeros(2*Ntrig);
CPM(:,1:2:end) = [eye(Ntrig); zeros(Ntrig)];
CPM(:,2:2:end) = [zeros(Ntrig); eye(Ntrig)];

% Columns corresponding to cos(16*theta) and sin(16*theta) will get dropped
dropCols = [Ntrig+1, 2*Ntrig+2];

% Arrange Current Pattern Matrix
J(:,dropCols) = [];         % Delete columns for cos(16*theta) and sin(16*theta)
J =  J * CPM;               % Permute columns to match CP convention

% Prepare all 3 voltage matrices. For each matrix we need to:
% -Drop columns for cos(16*theta) and sin(16*theta),
% -Apply the column permutation matrix,
% -Ensure voltages sum to zero in each column to satisfy Kirchoff's law.
V_1cem(:,dropCols) = [];
V_1cem = V_1cem * CPM;
adj = sum(V_1cem)/L;
V_1cem = V_1cem - adj(ones(L,1),:);

V_clb(:,dropCols) = [];
V_clb = V_clb * CPM;
adj = sum(V_clb)/L;
V_clb = V_clb - adj(ones(L,1),:);

V_trg(:,dropCols) = [];
V_trg = V_trg * CPM;
adj = sum(V_trg)/L;
V_trg = V_trg - adj(ones(L,1),:);


% Raw ND maps, size 2*Ntrig x 2xNtrig
ND1_cem  =   V_1cem' * J;
ND_clb   =   V_clb' * J;
ND_trg   =   V_trg' * J;

% Continuum-model ND & DN maps corresponding to sigma = 1
ND1 = diag(reshape(cat(1,1./(1:Ntrig),1./(1:Ntrig)),1,[]),0);
DN1 = diag(reshape(cat(1,1:Ntrig,1:Ntrig),1,[]),0);

% Compute a symmetric, noise-robust DN map using calibration formula
DN_right = inv(ND1_cem / ND_clb * ND_trg - ND1_cem + ND1);
DN = (DN_right + DN_right')/2;

% Add appropriate zero row and zero column corresponding to constant basis functions
DN = [zeros(2*Ntrig,1), DN];
DN = [zeros(1,2*Ntrig+1); DN];
DN1 = [zeros(2*Ntrig,1), DN1];
DN1 = [zeros(1,2*Ntrig+1); DN1];

end

