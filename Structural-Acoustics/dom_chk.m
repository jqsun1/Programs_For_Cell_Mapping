function [pf,ps,iRemv] = dom_chk(rel_pf, rel_ps)
%dom_chk does the dominance check on the inputs
%   The inputs include the 
% rel_pf: The set to be checked for dominancy
% rel_ps: The corresponding points in variable space for pareto front
% rel_sol: The cell numbers correspoding to those points in pareto front
%   The outputs
% pf,ps,sol_cells: The pareto set, front and corresponding cell numbers
% after the dominated solutions are removed.
% Author: Yousef Naranjani 3/3/2014
nObj = size(rel_pf,2);
[pf, ind] = sortrows(rel_pf); %Order the pareto front like dictionary
len = size(rel_pf,1);
iRemv = [];     % index of the entires to be removed
ikeep = 1;      % index of the entries to be kept
chk_pop = 1;    % the population of cells that current cell needs to be compared to
for i = 2 : len
    flag = zeros(chk_pop,1);
    for k = chk_pop : -1 : 1
        if all(pf(i,:) >= pf(ikeep(k),:)) && ~all(pf(i,:)==pf(ikeep(k),:)) 
            flag(k) = 1;
            break
        end
    end
    if any(flag)
        iRemv = [iRemv;i];
    else
        chk_pop = chk_pop+1;
        ikeep = [ikeep;i];
    end
end
pf(iRemv,:) = [];   % Removing the dominated cells
ps = rel_ps(ind,:);
ps(iRemv,:) = [];