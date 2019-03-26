function [c, N]= refine(apg, lb, ub , N, div)
%Subdivision
%
h = (ub-lb)./N;
z = arrayfun(@(x)celltoz(x, N), apg,'UniformOutput',false);
z = cell2mat(z')'; %store coordinates of every cell (nz*3)
%
n = div;
div = 2*div + 1;
%
zn = [];
rang = {};
%
h = h ./ div;      %Divide h by div
N = N .* div;      %Compute dimension of N
for i=1:length(div)
    rang{i} = -n(i):n(i);
end
comb = allcomb(rang{:});
%
for i = 1:size(z, 1)
    center = z(i,:)'.*(2*n + 1) - n;
    ncell = bsxfun(@plus, comb, center');
    zn = [zn; ncell];
end
z = zn;
%
cells = zeros(length(z), 1);
for i=1:length(z)
    cells(i) = ztocell(z(i,:)', N);
end
c = cells;