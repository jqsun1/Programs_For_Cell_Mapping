function py = FiniteTimeLorenz(y, f, ts)
[~, x] = ode45(f,[0 ts],y);
py = x(end,:)';