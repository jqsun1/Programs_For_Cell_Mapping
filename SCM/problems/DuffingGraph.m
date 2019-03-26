function py = DuffingGraph(y, k, alpha, B, tf)
% -------------------------------------------------------------------------
% Build the point to point map of duffing example by integrating at every
% period of external excitation. This is the most straightforward method in
% building discrete mapping in pointwise space.
%
% It is found that the fixed period sampling is quiet subjective. Different
% integration time would return totally different local mapping and thus
% lead to different attraction basin. Hence, a more local scheme is carried
% out by just examing the adjacent cell of each initial condition cell.
% -------------------------------------------------------------------------
% [~, x] = ode45(@(t,x)Duffing(t,x,k,alpha,B), [0 tf], y);
[~, x] = ode45solve(@(t,x)Duffing(t,x,k,alpha,B), 0:tf/50:tf, y);
py = x(end,:)';