function GCM = gcmDataBase_aug(N, lb, ub, div, ds, iter)
% Return the cell array of GCM database, the database
% consists three arrays that are defined in section 11.1
% Sink cells are specially treated as the (Nc+1)th row of
% GCM cell array. A pre-image array is added in GCM, which
% turns this code as the augement of traditional GCM.
%
if nargin < 6    
    Nc = prod(N);
    GCM = cell(Nc,4);
    S = 1:Nc;
else
    GCM = cell(length(iter),4);
    S = iter;
end
%
for i = 1:length(S)
    [I, C, P] = OneCellGCM(S(i), N, lb, ub, div, ds, S);
    GCM{i,1} = I;
    GCM{i,2} = C;
    GCM{i,3} = P;
end
%
% use directed graph to locate all pre-image cells in GCM
DG = gcm2graph(GCM, S); % directed graph
iDG = DG'; % inversed graph
for i = 1:length(S)
    GCM{i,4} = S(iDG(i,:));
end