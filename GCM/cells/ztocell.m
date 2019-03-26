function ncell = ztocell(z, N)
% compute cell index according to integer coordinate
% only those within the computational region are
% designated with a integer cell index.
%
if all(z-N<=0)
    opt = true; 
else
    opt = false;
end
%
if opt % within computational region
    z = z - 1;
    ncell = z(1);
    b = N(1);
    for i = 2:length(N)
        ncell = ncell + z(i)*b;
        b = b * N(i);
    end
    ncell = ncell + 1;
else
    % SCM sink cell representation
    ncell = 0; % cell locates outta the computational region
    %
%     % GCM sink cell representation
%     Nc = prod(N);
%     ncell = Nc + 1;
end