clear all
close all
clc

% Test of the quantitative analysis without subdivision.
% By: Free Xiong; 2014-08-14

N = [151;151];
No = 5;
[lb, ub, dsys] = ProblemDef(No);
samNo = 15;

% construct gcm
tic;
[P, scm] = gcm(dsys, lb, ub, N, samNo);
time_gcm = toc

% determine persistent group using scc
tic;
pg = pg_cells_scc2(P);
time_pg = toc

% iterating the limiting probability distribution
tic;
p_old = zeros(prod(N),1);
p_old(randperm(prod(N),1)) = 1; % initial condition on a cell
iter = 1:100;
err = zeros(size(iter));
for i = 1:length(iter)
    p_new = P'*p_old;
    err(i) = norm(p_new-p_old);
    p_old = p_new;
end
time_prob = toc

figure(1)
plot(iter,err)
xlabel('Iteration Time')
ylabel('$\epsilon$')

figure(2)
xc = PlotCells(pg, lb, ub, N, 'blue', 'EdgeOff');
scatter(xc(:,1),xc(:,2),8,p_new(pg),'square','filled')
axis([lb(1) ub(1) lb(2) ub(2)])