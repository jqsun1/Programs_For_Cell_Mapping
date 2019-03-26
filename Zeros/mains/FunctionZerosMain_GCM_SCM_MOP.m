clear all
close all
clc
% -------------------------------------------------------------------------
% GCM-SCM-MOP combined algorithm for higher dimensional finding zero
% problems. The code implements all ideas from SCM/GCM/MOP and is extended
% to a more general level. The basic structure is still two-stage:
%
% Stage one:
% SCM-GCM hybrid algorithm with gradient-base search for covering set, this
% step remains the same as I delveloped previously.
%
% Stage two:
% GCM-MOP direct search in an iterative manner for refinement.
%
% Stage three:
% If necassary, back to pointwise space with the help of fsolve command to
% pinpoint the final accurate solution.
% -------------------------------------------------------------------------
wq = cputime;
%
% up to now the most promising combination is (dopt = 1, sopt = 1 or 7)
sopt = 7; % stepsize control (stepsize.m)
dopt = 1; % local iteration control (dynamic_systems.m)
opt = false; % normalization control
%
N_old = [20; 20]; % cell space partition can be larger here since DS is used
div = ones(length(N_old),1); % (2n+1)
No = 8;
[lb, ub, dsys, g, f, fmop] = ProblemDef(No, N_old, dopt, sopt, opt);
%
max_iter = 3; % iteration times for SCM-MOP gradient-based

%% Stage one: GCM-SCM hybrid gradient-based search for covering set
S = covering_set_gcm(N_old, lb, ub, div, dsys, 1:prod(N_old), 'stable');

whos S
time_cover = cputime - wq

%% Stage two: GCM-MOP direct search for 'max_iter' times
mode = 3; % refined scm (inherit MOP)
N_gcm = N_old;
S_new = S;
for i = 1:max_iter
    S_new = reshape(S_new,length(S_new),1);
    [S_new, N_new] = refine(S_new, lb, ub , N_gcm, div);
    S_new = DS_covering_set_gcm_mop(fmop, N_new, lb, ub, g, S_new, 'stable');
    N_gcm = N_new;
end

whos S_new
N_new'
time = cputime - wq

%% Stage three: Back to pointwise space (optional)
% [xa, ~] = Stability_check(S_new, N_new, lb, ub, f, dsys);
% Q = block_classification(S_new, N_new);
% 
% Q = unique(Q);
% whos xa
% length(Q)
% time_Pointwise_final = cputime - wq