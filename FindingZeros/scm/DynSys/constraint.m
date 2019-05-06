function [c, ceq] = constraint(v, delta)
c = norm(v) - delta;
ceq = [];