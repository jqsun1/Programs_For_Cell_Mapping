function Dm = domicile(pg, gr, g, P)
% -------------------------------------------------------------------------
% Determine domiciles of each transient cell, the Dm is with the matrix
% structure of N-by-max(gr). Clearly, those with all zeros in Dm are the
% persistent cells or the cells that cannot reach to any persistent groups.
% Those persistent cells are labelled with all -1 at their domicile Dm
% entries.
%
% Input arguments:
%      pg:    pg of each cell, logical array with size N
%      gr:    pg group number of each cell, transient cells are with 0
%       g:    total number of pg, i.e., max(gr)
%       P:    sparse matrix representation of gcm
%
% Output argument:
%      Dm:    N-by-g matrix with each row as the domiciles of a cell, all
%             zeros of one row indicates that cell cannot reach any pg, all
%             -1 of one row represents that cell is a pg cell.
%
% By: Free Xiong; 2014-08-20
% -------------------------------------------------------------------------
Dm = zeros(length(pg),g); % g=max(gr)
for i = 1:g
    Pr_old = zeros(length(pg),1);
    Ti_old = zeros(length(pg),1);
    
    % Domi_init, set the original target set
    for z = 1:length(pg)
        if(gr(z)==i)
            Ti_old(z) = 1;
        end
    end
        
    % backward expansion
    n = 1;
    while(true)
        [Pr_new, Dm, Ti_new] = Domi_tran_i(Pr_old, Ti_old, i, P, Dm);
        if(norm(Pr_new-Pr_old)==0)
            break
        end
        Pr_old = Pr_new;
        Ti_old = Ti_new;
        n = n+1;
    end
    disp(['Backward step for pg ',num2str(i),' is ',num2str(n)]);
end
Dm(pg~=0,:) = -1;