function GCM = gcmDataBase(N, lb, ub, div, ds, iter)
% Return the cell array of GCM database, the database
% consists three arrays that are defined in section 11.1
% Sink cells are specially treated as the (Nc+1)th row of
% GCM cell array
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
    [I, C, P] = OneCellGCM(S(i), N, lb, ub, div, ds, S);
    GCM{i,1} = I;
    GCM{i,2} = C;
    GCM{i,3} = P;
end