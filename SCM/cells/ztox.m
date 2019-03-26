function xd = ztox(z, h, lb)
    xd = lb + h .* z - 1/2*h;       %Compute center point
end