function pcells = search_Pgroup(scm, max_period, s)
% -------------------------------------------------------------------------
% Simplified unravelling algorithm that determines the P-K cells remain to
% be refined further. This simplified version aviods the calculation of
% group and step numbers, which implies the interactions among different
% cells can be canelled. Hence, the code is very much suitable for parallel
% implementation.
%
% By: Free Xiong; 2014/02/14
% -------------------------------------------------------------------------
if nargin < 3
    s = 1:length(scm);
end
p = zeros(length(s),1); % period array
for i = 1:length(s)
    cell = s(i);
    % shoot forward with maximum expected period
    loop = zeros(max_period+1,1);
    for j = 1:max_period+1
        loop(j) = cell;
        cs = scm(s==cell); % image cell
        if cs(1) == s(i)
            p(ismember(s,loop(1:j))) = j; % period of the loop
            break
        else
            cell = cs(1);
        end
    end
end
% periodic cells are candidates
pcells = s(p~=0);