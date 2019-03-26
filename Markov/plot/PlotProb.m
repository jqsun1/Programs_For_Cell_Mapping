clear all
close all
clc

% plot limiting probabilty and persistent groups
% load plasma_prob_dist_N_45x45x45.mat

% calculate limiting probability of each PG using iteration
for j = 1:max(gr)
    Pi = P(pg(gr==j),pg(gr==j));
    xc = PlotCells(pg(gr==j), lb, ub, N, 'red', 'EdgeOff');
    p_old = zeros(length(pg(gr==j)),1); 
    p_old(1)=1;
    iter = 0;
    while(iter<150)
        p_new = Pi'*p_old;
        err = norm(p_new-p_old);
        if(err<1e-4)
            break
        end
        p_old = p_new;
    end
    if(length(lb)==2)
        scatter(xc(:,1),xc(:,2),10,p_new,'square','filled')
        axis([lb(1) ub(1) lb(2) ub(2)])
        hold on
    elseif(length(lb)==3)
        scatter3(xc(:,1),xc(:,2),xc(:,3),10,p_new,'square','filled')
        axis([lb(1) ub(1) lb(2) ub(2)])
        hold on
    end
end
view(160,30)

% print -depsc plasma_prob_dist_N_45x45x45.eps

clear all
load duffing_tran_N_100x100.mat
figure(2)
PlotCells(find(gr==1), lb, ub, N, 'red', 'EdgeOff');
hold on
PlotCells(find(gr>1), lb, ub, N, 'blue', 'EdgeOff');
hold on
PlotCells(find(Dm(:,1)~=0 & Dm(:,2)~=0 & Dm(:,3)~=0), lb, ub, N, 'green', 'EdgeOff');
% print -depsc duffing_trans_N_100x100.eps