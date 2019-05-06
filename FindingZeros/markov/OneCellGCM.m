function [I, C, P] = OneCellGCM(ncell, N, lb, ub, div, ds, S)
% Construct GCM for a given cell by sampling Method
% proposed in section 11.12
% Note that ds is called as [xnew, ~, ~, fold] = ds(xold, max(h))
%
% -------------------------------------------------------------------------
% I transplated the GCM code directly to here that helps build the Markov
% transitional matrix, sampling technique proposed in the Cell-to-Cell
% Mapping book is applied here for probability calculation.
%
% Input Arguments:
%      ncell:   cell index
%   N, lb, ub:  cell space info
%       div:    sampling division
%        ds:    dynamic system for point to point mapping
%        S:     cell set of our interest
%
% Output Arguments:
%      I:   number of image cells
%      C:   image cell array
%      P:   probability of each image cell
%
% Free Xiong: 2013/08/05 modification
% -------------------------------------------------------------------------
%
% image cells over refined space
[subcells, ~] = refine(ncell, lb, ub , N, div);
M = length(subcells);
%
div = 2*div + 1;
h = (ub - lb)./N;
h_div = h./div;
N_div = N.*div;
%
C_origin = zeros(M,1); % image cells of all sampled cells within ncell
for i = 1:M
    z = celltoz(subcells(i),N_div);
    xSubCell = ztox(z,h_div,lb);
    xImgCell = ds(xSubCell);
    z = xtoz(xImgCell,h,lb);
    C_origin(i) = ztocell(z,N);
    %
    if ~ismember(C_origin(i),S)
        C_origin(i) = 0; % sink cell in GCM
    end
end
%
% sampling
C = unique(C_origin); % image cell array of ncell
I = length(C); % No. of image cells of ncell
%
P = zeros(I,1);
for i = 1:I
    indices = find(C(i)==C_origin);
    Mi = length(indices);
    P(i) = Mi/M;
end