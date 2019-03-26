function [d, class, Q, P]=period(G, P)
% -------------------------------------------------------------------------
% Classify the persistent group (irreducible markov chain) into 'd'
% subgroups whose period is d.
%
% Input arguments:
%      G:     graph represented GCM with logical sparse matrix
%      P:     associated GCM probability matrix (optional)
%
% Output arguments:
%      d:     period of the irreducible markov chain
%    class:   subgroup No. of each state in the markov chain
%      Q:     canoical cyclic pattern of graph matrix
%      P:     canoical cyclic pattern of probability matrix (optional)
%
% Written by Free: Nov/12/2013
% -------------------------------------------------------------------------
if nargin < 2
    P = [];
end
%
n = length(G);
step = zeros(n,1); % shooting steps from node 1
%
% Note for irreducible markov chain, there is only one period. Thus,
% without loss of generality, we choose the inital shoot from node 1.
states = (1:n)';
current_states = 1;
processed_states = [];
level = 0; 
d = 0;
while true
    processed_states = [processed_states; current_states];
    img_states = [];
    for i = current_states'
        img_states = [img_states; states(G(i,:)~=0)];
    end
    %
    % level up, assign new step for image states
    level = level + 1;
    step(img_states) = level; % constant update, might cover old values
    %
    % check whether image states reach node 1 again
    if ismember(1,img_states)
        d = gcd(level, d);
        if ~all(ismember(img_states,processed_states))
            % delete node 1 to aviod infinite loop
            img_states(img_states==1) = [];
        else
            break
        end
    end
    %
    % swap states for new level of forward shooting
    current_states = img_states;
end
class = rem(step,d) + 1;
%
% permute positions of states for canoical form
[~, ix] = sort(class);
Q = G(ix,ix);
if ~isempty(P)
    P = P(ix,ix);
end