function z = xtoz(x, h, lb)
%Compute cell coordinates
z = round((x - lb) ./ h + 1/2);