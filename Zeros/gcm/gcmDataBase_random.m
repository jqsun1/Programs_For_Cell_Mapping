function GCM = gcmDataBase_random(N, lb, ub, k, ds, iter)
% GCM database with kFF sampling during each cell mapping construction.
%
if nargin < 6    
    Nc = prod(N);
    GCM = cell(Nc,3);
    S = 1:Nc;
else
    GCM = cell(length(iter),3);
    S = iter;
end
%
for i = 1:length(S)
    [I, C, P] = OneCellGCM_random(S(i), N, lb, ub, k, ds, S);
    GCM{i,1} = I;
    GCM{i,2} = C;
    GCM{i,3} = P;
end