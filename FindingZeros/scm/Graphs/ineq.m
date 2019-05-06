function cell = ineq(g, x, cell)
% default goal of inequality constraint is to set it less then zero
    if ~isempty(g)
        if any((g(x) <= 0) == 0)
            cell = 0;
        end
    end
end