function ncell = ztocell(z, N)
% compute cell index according to integer coordinate
% only those within the computational region are
% designated with a integer cell index.
%
if all(z-N<=0) && all(z > 0)  % JQS
    z = z - 1;
    ncell = z(1);
    b = N(1);
    for i = 2:length(N)
        ncell = ncell + z(i)*b;
        b = b * N(i);
    end
    ncell = ncell + 1;
else
%     ncell = 0;
    ncell = prod(N) + 1;  % The sink cell is the last one, JQS
end