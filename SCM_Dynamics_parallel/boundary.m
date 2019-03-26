function mcells = boundary(gr,N,cell,S)
% -------------------------------------------------------------------------
% Boundary cells extracted from different group numbers, the group number
% is acquired from scm, with global searching, note the input cell should
% contain unstable saddle point
%
% By: Free Xiong; 2015-02-27
% -------------------------------------------------------------------------
if nargin < 4
    S = (1:length(gr))';
end
tgt = zeros(length(gr),1);
bn_old = zeros(length(gr),1);
bn_new = zeros(length(gr),1);

% continuation like extraction
tgt(cell==S) = 1;
bn_old(cell==S) = 1;
bn_new(cell==S) = 1;
while true
    tgt_cells = S(tgt==1);
    for i = 1:length(tgt_cells)
        cs = tgt_cells(i);
        ncells = neighbour_finder(N,cs);
        for j = 1:length(ncells)
            nc = ncells(j);
            if ismember(nc,S)
                if (gr(nc==S)~=gr(cs==S)) && bn_old(nc==S)==0
                    % new boundary cell brought in
                    tgt(nc==S) = 1;
                    bn_new(nc==S) = 1;
                end
            end
        end
        tgt(cs==S) = 0; % cs is already been processed
    end
    
    if norm(bn_new-bn_old)==0
        break
    end
    
    bn_old = bn_new;
end

% manifold cells
mcells = S(bn_old==1);