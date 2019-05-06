function [I, C, P] = OneCellGCM(ncell, N, lb, ub, div, ds, S)
% Construct GCM for a given cell by sampling Method
% proposed in section 11.12.
%
Nc = prod(N);
if nargin < 7
    S = 1:Nc;
end
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
C_origin = zeros(M,1); % image cell array
for i = 1:M
    z = celltoz(subcells(i),N_div);
    xSubCell = ztox(z,h_div,lb);
    xImgCell = ds(xSubCell);
    z = xtoz(xImgCell,h,lb);
    C_origin(i) = ztocell(z,N);
    %
    if ~ismember(C_origin(i),S)
        C_origin(i) = 0; % sink cell, we take them
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