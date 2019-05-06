function py = DuffingGraph(y, k, alpha, B, ub, lb, N)
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
[~, x] = ode45(@(t,x)Duffing(t,x,k,alpha,B), [0 2*pi], y);
% [~, x] = ode45solve(@(t,x)Duffing(t,x,k,alpha,B), 0:2*pi/25:2*pi, y);
%
% final point approach
py = x(end,:)';
%
% % adjacent cell approach
% seq = zeros(size(x,1),1);
% h = (ub - lb)./N;
% for i = 1:length(seq)
%     z = xtoz(x(i,:)',h,lb);
%     seq(i) = ztocell(z,N);
% end
% %
% seq = unique(seq,'stable'); % shooting sequence in cell space
% if length(seq) > 1
%     img = seq(2); % image cell is 2nd along the sequence (adjacent)
% else
%     img = seq; % attractor
% end
% z = celltoz(img,N);
% py = ztox(z,h,lb);