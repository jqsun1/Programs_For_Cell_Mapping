function S_final = Forward_shoot_eliminate(S_old, C, Nm)
% -------------------------------------------------------------------------
% For a given cell set S_old and a simple cell mapping C of S_old, if an
% element of S_old can generate a sequence with Nm times forward mapping,
% we retain that element in S_old, otherwise, we kill it. This code aims at
% implimenting the Nm times forward shooting by using only one time cell
% mapping. It avoids the previous approach calling the graph subroutine and
% conduct point-to-point iteration every time.
%
% Input arguments:
%     S_old:   Old set remain to be trimed
%       C:     SCM of S_old
%      Nm:     Times of forward shooting
%
% Output argument:
%     S_final: Final subset of S_old survives the Nm times of shooting
% -------------------------------------------------------------------------
%
mRmv = [];
length(S_old)
for m = 1:length(S_old)
    i = m;
    % generate Nm times forward shooting sequence
    for n = 1:length(Nm)
        j = find(S_old == C(i));
        if isempty(j)
            mRmv = [mRmv; m];
            break
        else
            i = j;
        end
    end
end
%
S_old(mRmv) = [];
S_final = S_old;   