clear all
close all
clc

!nvcc -ptx scm_mop.cu
k = parallel.gpu.CUDAKernel('scm_mop.ptx','scm_mop.cu','cellEvaluation');
k.ThreadBlockSize = [1024 1 1];
k.GridSize = [512 1 1];

dim = 3;
obj = 2;

N = [10;10;10];
tot_cells = prod(N);
cells = 1:prod(N);
lb = [0;0;0];
ub = [40;40;40];
xc = zeros(tot_cells*dim,1);
fe = zeros(tot_cells*obj,1);

N = gpuArray(int32(N));
cells = gpuArray(int32(cells));
lb = gpuArray(single(lb));
ub = gpuArray(single(ub));
xc = gpuArray(single(xc));
fe = gpuArray(single(fe));

[xc,fe] = feval(k,cells,lb,ub,N,tot_cells,xc,fe);