function SCM = SelectedSCM(N, GCM, iter, choice)
% -------------------------------------------------------------------------
% Select a SCM for persistent group identification which is compatiable 
% with the generated GCM, the SCM is selected based on section 11.3. Note
% GCM is input in the form as the cell array.
% -------------------------------------------------------------------------
if nargin < 3
    Nc = prod(N);
    SCM = zeros(Nc,1);
    S = 1:Nc;
    choice = 'max';
else
    SCM = zeros(length(iter),1);
    S = iter;
end
%
for i = 1:length(S)
    % define image cell choosing model
    if strcmp(choice,'max')
        [~, index] = max(GCM{i,3});
    elseif strcmp(choice,'min')
        [~, index] = min(GCM{i,3});
    else
        I = GCM{i,1};
        index = Roulette_wheel(1:I,1/I*ones(1,I));
    end
    %
    if isempty(GCM{i,2})
        SCM(i) = 0;
    else
        SCM(i) = GCM{i,2}(index);
    end
end