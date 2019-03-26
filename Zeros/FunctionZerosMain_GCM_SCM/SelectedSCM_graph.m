function SCM = SelectedSCM_graph(N, DG, P, iter, choice)
% -------------------------------------------------------------------------
% Select a SCM for persistent group identification which is compatiable 
% with the generated GCM, the SCM is selected based on section 11.3. Note
% GCM is input in the form of spare matrix graph.
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
    % define image cell choosing mode
    img = S(DG(i,:)); % image cells
    temp = P(i,:);
    pro = temp(DG(i,:)); % transitional probability
    if isempty(img)
        SCM(i) = 0;
    else
        if strcmp(choice, 'max')
            [~, index] = max(pro);
        elseif strcmp(choice, 'min')
            [~, index] = min(pro);
        elseif strcmp(choice, 'rand')
            index = randperm(length(img),1);
        else
            I = length(img);
            index = Roulette_wheel(1:I,1/I*ones(1,I));
        end
        SCM(i) = img(index);
    end
end