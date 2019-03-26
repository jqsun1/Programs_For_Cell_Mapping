function [Pr_new, Dm, Ti_new] = Domi_tran_i(Pr_old, Ti_old, i, P, Dm)
% -------------------------------------------------------------------------
% Determine the transient cells with domcile at the ith pg. Backward
% searching with expanding traget set is applied here.
%
% Input arguments:
%      Pr_old:   processed cell indicator
%      Ti_old:   target cell indicator
%        i:      persistent group number
%        P:       gcm transitional matrix
%       Dm:     domicile indicator with (i,j) entry as the jth pg of cell i
%
% By: Free Xiong; 2014-08-20
% -------------------------------------------------------------------------
Pr_new = Pr_old;
Ti_new = Ti_old;
for z = 1:size(P,1)
    if(Pr_old(z)==0 && Ti_old(z)==1)
        Pr_new(z) = 1;
        pimg = find(P(:,z)~=0);
        for j = 1:length(pimg)
            w = pimg(j);
            Ti_new(w) = 1;
            Dm(w,i) = i;
        end
    end
end