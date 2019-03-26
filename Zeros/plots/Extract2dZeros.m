% clear all
% close all
% clc

% extract zeros by using clutersing analysis for the 2d example

% load('FunctionZeros_2I2O_N_13x13_iter_4.mat')
% load('FunctionZeros_Fractional_N_13x13_iter_4.mat')

xc = PlotCells(S_new,lb,ub,N_new,'b','EdgeOff');
close
point_sol = []; % individual representative point for each cluster
sol_num = 1; % cluster number

options = optimoptions(@fsolve,'Display','off',...
    'MaxIter',1000,'MaxFunEvals',1000,'TolFun',1e-10,'TolX',1e-12);
[xf,~, flag] = fsolve(f,xc(1,:),options);
point_sol = [point_sol; xf];
delta = 0.001; % radius of the neighbourhood ball
for i = 2:size(xc,1)
    [xf,~, flag] = fsolve(f,xc(i,:),options);
    if any(xf'<lb) || any(xf'>ub) || flag~=1
        continue;
    end
    disjoint = true;
    for j = 1:sol_num
        if sqrt(sum((xf-point_sol(j,:)).^2)) <= delta
            disjoint = false;
            break;
        end
    end
    if disjoint
        point_sol = [point_sol; xf];
        sol_num = sol_num+1;
    end
end
