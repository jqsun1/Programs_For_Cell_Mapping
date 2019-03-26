function [C, S] = graph(dysy, ub, lb, N, S)
% -------------------------------------------------------------------------
% Construct the simple cell mapping from point mapping
%
% By: Furui Xiong; 10/24/2015
% -------------------------------------------------------------------------
h = (ub - lb) ./ N;
nz = prod(N);                 
if nargin < 5
    S = (1:nz);
end
C = zeros(length(S), 1);

for i = 1:length(S)
    cs = S(i);
    
    z = celltoz(cs, N); % Compute coordinates of cs
    xa = ztox(z,h,lb); % Compute central point of cs
    xd = dysy(xa);      % Point evolution
    z = xtoz(xd,h,lb); %Compute cell coordinates
    ncell = ztocell(z, N);
    
    if ~isempty(find(S == ncell, 1))
        C(i) = ncell; % cs is not a sink cell
    else
%         C(i) = 0;
        C(i) = nz+1;
    end
end