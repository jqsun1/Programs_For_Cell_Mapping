clear all
close all
clc

% point mapping for attractor finding
% by: furui 10/28/2015

tic;
addpath(genpath('./')) % Add all the sub-folders of the current directory

% define problem, the short description of each problem is listed here:
% 1, duffing oscillator with external excitation
% 2, van-de-pol equation
% 3, henon map, https://en.wikipedia.org/wiki/H%C3%A9non_map
% 4, duffing oscillator without excitation
% 5, discrete duffing map, https://en.wikipedia.org/wiki/Duffing_map
% 6, discrete gingerbreadman map, https://en.wikipedia.org/wiki/Gingerbreadman_map
% 7, discrete ikeda map, https://en.wikipedia.org/wiki/Ikeda_map
% 8, discrete sun's map, from "Effects of small random uncertainties on non-linear systems studied by the generalized cell mapping method"
% 9, another example from the same reference above
tic;
prob_No = 9;
[lb, ub, N, dsys] = ProblemDef(prob_No);

% traverse every central point for point mapping over a trajectory
xc = PlotCells(1:prod(N), lb, ub, N, 'red', 'EdgeOff');
close

% choose which type of sampling in state space
type = 'random';
% type = 'uniform';

if strcmp(type,'uniform')
    %% uniform sampling
    iter = 1000; % iteration steps
    ikeep = 20; % keep the last few points in one trajectory
    xa = zeros(prod(N)*ikeep,length(lb));
    for i = 1:prod(N)/2
        x_old = xc(i,:);
        xtemp = zeros(ikeep,length(lb));
        id = 1;
        for j = 1:iter
            % store the last few points
            if j>iter-ikeep
                xtemp(id,:) = x_old;
                id = id+1;
            end
            x_new = dsys(x_old);
            x_old = x_new;
        end
        xa((i-1)*ikeep+1:i*ikeep,:) = xtemp;
    end
elseif strcmp(type,'random')
    %% random sampling
    samNo = 10000;
    iter = 1000;
    ikeep = 200;
    xs = repmat(lb',samNo,1) + rand(samNo,length(lb)).*repmat(ub'-lb',samNo,1);
    xa = zeros(samNo*ikeep,length(lb));
    for i = 1:samNo
        x_old = xs(i,:);
        xtemp = zeros(ikeep,length(lb));
        id = 1;
        for j = 1:iter
            % store the last few points
            if j>iter-ikeep
                xtemp(id,:) = x_old;
                id = id+1;
            end
            x_new = dsys(x_old);
            x_old = x_new;
        end
        xa((i-1)*ikeep+1:i*ikeep,:) = xtemp;
    end
end

time = toc

xa = unique(xa,'rows');
plot(xa(:,1),xa(:,2),'r.')
axis([lb(1) ub(1) lb(2) ub(2)]);
box on