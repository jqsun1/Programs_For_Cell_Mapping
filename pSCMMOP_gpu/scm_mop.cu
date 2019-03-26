/*
CUDA kernels for hybrid scm-mop on gpu
Free Xiong
2014/02/27
*/
#include "cuda_runtime.h"
#include <cstdlib>
#include "cell.cuh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>

using namespace std;

#define NumberOfBlock 2048
#define ThreadPerBlock 1024

//evaluate cells of central points and function values
__global__ void cellEvaluation(const int * cells, const float * lb, const float * ub, 
									 const int * N, const int tot_cells, float * xc, float * fe){
	/*
	Input arguments:
	     cells:    cells waiting to be evaluated
		 lb,ub:    lower and upper bound of searching domain
		   N:      celluar space partition
	   tot_cells:  total number of cells in "cells" array

	Output arguments:
	      xc:      1D array with cell central coordinates
		  fe:      1D array with central function values
	*/
	int id = threadIdx.x + blockIdx.x*blockDim.x; //unique id of thread
	int i,j;
	int z[D];
	int cs;
	float h[D];
	float xd[D];
	float fd[O];
	
	for(i=0;i!=D;i++)
		h[i] = (ub[i]-lb[i])/N[i];

	i = id;
	while(i<tot_cells){
		cs = cells[i];
		celltoz(z,cs,N);
		ztox(xd,z,h,lb);
		f(fd,xd);
		for(j=0;j!=D;j++)
			xc[i*D+j]=xd[j];
		for(j=0;j!=O;j++)
			fe[i*O+j]=fd[j];
		i+=blockDim.x*gridDim.x;
	}
}

//creat simple cell mapping
__global__ void scm_ds(const int* cells, const int tot_cells, const float * xc, 
					   const float * fe, const float* lb, const float* ub, const int* N, 
					   int* img){
	//scm with directed search
	/*
	Input arguments:
	cells:     input cell set
	tot_cells: number of cells in the set
	lb,ub:     lower and upper bound
	N:         cell space partition
	xc,fe:     cell info evalauted outside

	Output argument:
	img:      image cells
	*/

	//scm
	int id = threadIdx.x + blockIdx.x*blockDim.x;
	int i,j,k,cs;
	bool flag,better,dominated,constraints;
	int index;
	float xcs[D],xn[D],fc[O],fn[O],cst[C];
	float cmp,df_old,df_new;
	int neighbour[pop];
	//int z_c[D], z_fwd[D], z_bwk[D];

	//for(i=0;i<D;i++)
	//	pop*=3;

	// loop over the entire cell set
	i = id;
	while(i<tot_cells){
		cs = cells[i];
		for(j=0;j<D;j++)
			xcs[j] = xc[i*D+j];// central cell coordinates
		for(j=0;j<O;j++)
			fc[j] = fe[i*O+j];// central cell function values
		constraints = true;

		//check central cell constaints violation
		g(cst,xcs);
		for(j=0;j<C;j++){
			if(cst[j]>0){
				constraints = false;
				break;
			}
		}
		if(!constraints){
			img[i] = 0; //sink cell
			i+=blockDim.x*gridDim.x;
			continue;
		}

		neighbour_finder(neighbour,N,cs); //be aware of the visitation of neighbour array

		////search among neighbours that are along the coordinate directions (2*dim)
		//celltoz(z_c,cs,N);
		//for(j=0;j<D;j++){
		//	z_bwk[j] = z_c[j];
		//	z_fwd[j] = z_c[j];
		//}
		//for(j=0;j<D;j++){
		//	z_fwd[j]+=1; //increase along the jth direction
		//	z_bwk[j]-=1;
		//	if(z_fwd[j]>N[j])
		//		neighbour[j] = 0; //out of upper boundary
		//	else
		//		neighbour[j] = ztocell(z_fwd,N);
		//	if(z_bwk[j]<1)
		//		neighbour[j+D] = 0; //out of lower boundary
		//	else
		//		neighbour[j+D] = ztocell(z_bwk,N);
		//	z_fwd[j] = z_c[j]; //recover the coordinate for one direction increasing only
		//	z_bwk[j] = z_c[j];
		//}
		
		df_old = -1.0f;
		dominated = false;
		for(j=0;j<pop;j++){
			//check for valid neighour and compare
			flag = false;
			for(k=0;k<tot_cells;k++){
				if(cells[k]==neighbour[j]){
					flag = true;
					index = k;
					break;
				}
			}
			if(!flag)
				index = -1;
			
			//dominance optimality when picking the image cell
			if(neighbour[j]==0 || !flag)
				continue;
			else{
				for(k=0;k<O;k++)
					fn[k] = fe[index*O+k]; //neighbour function values
				
				//all(fn<fc)
				cmp = 10.0f;
				better = true;
				for(k=0;k<O;k++){
					cmp = fn[k]-fc[k];
					if(cmp>0.0f){
						better = false;
						break;
					}
				}

				//steepest descent
				if(better){
					dominated = true; //there is better neighbour
					df_new = 0.0f;
					for(k=0;k<O;k++)
						df_new+=fc[k]-fn[k];
					if(df_new>df_old){
						img[i] = neighbour[j];
						df_old = df_new; //swap for next df comparison
					}
					//check violation of constraints of the neighbours
					for(k=0;k<D;k++)
						xn[k] = xc[index*D+k]; //neighbour center
					g(cst,xn);
					for(k=0;k<C;k++){
						if(cst[k]>0){
							img[i] = 0; //constraints violation, sink cell
							break;
						}
					}
				}
			}
		}

		if(!dominated)
			img[i] = cs;

		i+=blockDim.x*gridDim.x;
	}
}

//scm unravelling for p-k cells
__global__ void scm_search(const int * cells, const int * img, const int tot_cells, 
						   const int max_period, int * period){
	//scm search for periodic cells
	/*
	Input arguments:
	cells:        cell set
	img:          image of cells
	tot_cells:    total number of cell set
	max_period:   maximum period extracted, mostly p-2

	Ouput argument:
	period:       period of each cell, transients cell has -1 period
	*/
	int id = threadIdx.x + blockIdx.x*blockDim.x;
	int i,j,k,cs,cell_old,cell_new,index;
	bool pkcell;
	
	//loop over entire cell set
	i = id;
	while(i<tot_cells){
		cs = cells[i];
		cell_old = cs;
		pkcell = false;

		//shoot forward to extract cyclic structures
		for(j=0;j<max_period+1;j++){
			//find index
			for(k=0;k<tot_cells;k++){
				if(cells[k]==cell_old){
					index = k;
					break;
				}
			}
			//shoot forward
			cell_new = img[index];
			if(cell_new==0){
				period[i] = -1; // meet sink cell
				break;
			}
			else if(cell_new==cs){
				period[i] = j+1;
				pkcell = true;
				break;
			}
			cell_old = cell_new;
		}

		//mark transient cells
		if(!pkcell)
			period[i] = -1;

		i+=blockDim.x*gridDim.x;
	}
}

//parallel refinement
__global__ void refine(int * rcells, const int * cells, const int tot_cells, int * N_new,
					   const int * N, const int * div, const int * allcomb, const int dim){
	int id = threadIdx.x + blockIdx.x*blockDim.x;
	int i,cs;
	i = id;

	//loop over candidate cell set
	while(i<tot_cells){
		cs = cells[i];
		subdivision(&rcells[dim*i],cs,N,N_new,dim,div,allcomb);
		i+=blockDim.x*gridDim.x;
	}
}

//parallel dominance check
__global__ void dominance(const float * fe, const int * period, const int total_cells, int * flag){
	//dominance check over the candidate cells
	/*
	Input arguments:
	fe:       function values of all cells
	period:   periodic number of all cells
	tot_cells:  total number of cell set

	Output argument:
	flag:     indicator of whether a cell is a Pareto cell (0,1)
	*/
	int id = threadIdx.x + blockIdx.x*blockDim.x;
	int i,j,k;
	float fc[O],cmp[O];
	int temp[O],gz;
	bool dominated, allequal;

	//loop over entire cell set
	i = id;
	while(i<total_cells){
		dominated = false;
		if(period[i]==-1){
			flag[i] = 0;
			i+=blockDim.x*gridDim.x;
			continue;
		}
		else{
			for(j=0;j<O;j++)
				fc[j] = fe[i*O+j];
			//compare function values with other candidates
			for(j=0;j<total_cells;j++){
				if(period[j]==-1)
					continue;
				else{
					allequal = true;
					//get the function value difference and compare (temp=1 is bad)
					for(k=0;k<O;k++){
						cmp[k] = fc[k]-fe[j*O+k];
						if(cmp[k]>=0)
							temp[k]=1;
						else
							temp[k]=0;
					}
					//check if all objectives are equal happens
					for(k=0;k<O;k++){
						if(cmp[k]!=0){
							allequal = false;
							break;
						}
					}
					//check whether (all(cmp>=0) && !all(cmp==0))  stands
					gz=1;
					for(k=0;k<O;k++)
						gz*=temp[k];
					if(gz==1 && !allequal){
						dominated = true;
						flag[i] = 0;
						break;
					}
				}
			}
		}
		if(!dominated)
			flag[i] = 1;

		i+=blockDim.x*gridDim.x;
	}
}

void saveData(float * xc, float * fe, int write_cells, int * status){
	//save final Pareto set and front into different files
	int i,j;
	ofstream outData;
	int counter=0;

	outData.open("ps.dat");
	for(i=0;i<write_cells;i++){
		if(status[i]==1){
			counter++;
			for(j=0;j<D;j++)
				outData<<scientific<<setw(15)<<xc[i*D+j]<<setprecision(4);
			outData<<endl;
		}
	}
	outData.close();

	outData.open("pf.dat");
	for(i=0;i<write_cells;i++){
		if(status[i]==1){
			for(j=0;j<O;j++)
				outData<<scientific<<setw(15)<<fe[i*O+j]<<setprecision(4);
			outData<<endl;
		}
	}
	outData.close();
	
	cout<<endl;
	cout<<"Number of Pareto cells is: "<<counter<<endl;
}

void scm_mop(int * N, float * lb, float * ub, int * div, int max_period, int max_iter){
	/*
	main structure of scm-mop hybrid algorithm, this code execute kernels
	from host (CPU) to device (GPU)
	*/
	int i,j,iter,tot_cells=1,dim=1;
	int write_cells;
	int * cells;
	float * device_xc, * device_fe;
	int * device_N, * device_div;
	float * device_lb, * device_ub;
	int * device_cells, * device_img, * device_period;
	int * device_candcells;
	int candNo;
	int * device_Nnew, * device_rcells;
	int * allcomb, * device_allcomb;
	int * device_flag;
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start);

	for(i=0;i<D;i++){
		tot_cells *= N[i];
		dim *= 2*div[i]+1;
	}

	allcomb = new int [dim*D];
	cartprod(allcomb,div); //cartesin product of cell coord increasing
	
	cudaMalloc(&device_N, D*sizeof(int));
	cudaMalloc(&device_div, D*sizeof(int));
	cudaMalloc(&device_lb, D*sizeof(float));
	cudaMalloc(&device_ub, D*sizeof(float));
	cudaMalloc(&device_allcomb, dim*D*sizeof(int));
	cudaMalloc(&device_Nnew, D*sizeof(int));

	cudaMemcpy(device_N, N, D*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(device_div, div, D*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(device_lb, lb, D*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(device_ub, ub, D*sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(device_allcomb, allcomb, D*dim*sizeof(int), cudaMemcpyHostToDevice);

	iter = 0;

//--start of the loop
	while(iter<max_iter){
		//intitial cell set
		if(iter==0){
			cells = new int[tot_cells];
			for(i=0;i<tot_cells;i++)
				cells[i] = i+1; //1-base cell number counting
			cudaMalloc(&device_cells, tot_cells*sizeof(int));
			cudaMemcpy(device_cells, cells, tot_cells*sizeof(int), cudaMemcpyHostToDevice);
			delete[] cells;
		}

		write_cells = tot_cells;

		//evaluate cells
		cudaMalloc(&device_xc, tot_cells*D*sizeof(float));
		cudaMalloc(&device_fe, tot_cells*O*sizeof(float));
		cellEvaluation<<<NumberOfBlock,ThreadPerBlock>>>(device_cells,device_lb,device_ub,device_N,tot_cells,device_xc,device_fe);
		cudaDeviceSynchronize();

		//build scm
		cudaMalloc(&device_img, tot_cells*sizeof(int));
		scm_ds<<<NumberOfBlock,ThreadPerBlock>>>(device_cells,tot_cells,device_xc,device_fe,device_lb,device_ub,device_N,device_img);
		cudaDeviceSynchronize();
		//if(iter==0){
		//	int * temp = new int [tot_cells];
		//	cudaMemcpy(temp, device_cells, tot_cells*sizeof(int), cudaMemcpyDeviceToHost);
		//	for(i=0;i<100;i++){
		//			cout<<setw(5)<<temp[i];
		//		cout<<endl;
		//	}
		//	delete[] temp;
		//	//system("pause");
		//}

		//unravel scm
		cudaMalloc(&device_period, tot_cells*sizeof(int));
		scm_search<<<NumberOfBlock,ThreadPerBlock>>>(device_cells,device_img,tot_cells,max_period,device_period);
		cudaDeviceSynchronize();

		//collect candidate cells on CPU
		int * period = new int [tot_cells];
		cudaMemcpy(period, device_period, tot_cells*sizeof(int), cudaMemcpyDeviceToHost);
		candNo = 0;
		for(i=0;i<tot_cells;i++){
			if(period[i]!=-1)
				candNo++;
		}
		cout<<"Iteration "<<iter+1<<" complete. ";
		cout<<candNo<<" cells found..."<<endl;

		int * candcells = new int [candNo];
		j=0;
		cells = new int [tot_cells];
		cudaMemcpy(cells, device_cells, tot_cells*sizeof(int), cudaMemcpyDeviceToHost);
		for(i=0;i<tot_cells;i++){
			if(period[i]!=-1){
				candcells[j] = cells[i];
				j++;
			}
		}
		delete[] cells;

		cudaMalloc(&device_candcells, candNo*sizeof(int));
		cudaMemcpy(device_candcells, candcells, candNo*sizeof(int), cudaMemcpyHostToDevice);
		delete[] candcells;
		delete[] period;
	
		//refine
		if(iter<max_iter-1){
			cudaMalloc(&device_rcells, candNo*dim*sizeof(int));
			refine<<<NumberOfBlock,ThreadPerBlock>>>(device_rcells,device_candcells,candNo,device_Nnew,device_N,device_div,device_allcomb,dim);
			cudaDeviceSynchronize();

			//ready for next iteration
			tot_cells = candNo*dim;
			cudaFree(device_cells);
			cudaMalloc(&device_cells, candNo*dim*sizeof(int));
			cudaMemcpy(device_cells, device_rcells, candNo*dim*sizeof(int), cudaMemcpyDeviceToDevice);
			cudaMemcpy(device_N, device_Nnew, D*sizeof(int), cudaMemcpyDeviceToDevice);

			cudaFree(device_rcells);
			cudaFree(device_period);
			cudaFree(device_xc);
			cudaFree(device_fe);
		}

		iter++;

		cudaFree(device_candcells);
		cudaFree(device_img);
	}
//--end of the loop

	//dominance check
	cudaMalloc(&device_flag, write_cells*sizeof(int));
	dominance<<<NumberOfBlock,ThreadPerBlock>>>(device_fe,device_period,write_cells,device_flag);
	cudaDeviceSynchronize();
	cudaEventRecord(stop);

	//write Pareto results to file
	float * xc = new float [write_cells*D];
	float * fe = new float [write_cells*O];
	int * flag = new int [write_cells];
	cudaMemcpy(flag, device_flag, write_cells*sizeof(int), cudaMemcpyDeviceToHost);
	cudaMemcpy(xc, device_xc, write_cells*D*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(fe, device_fe, write_cells*O*sizeof(float), cudaMemcpyDeviceToHost);

	cudaEventSynchronize(stop);
	float gpu_time;
	cudaEventElapsedTime(&gpu_time, start, stop);
	cout<<endl;
	cout<<"Device (GPU) runtime is: "<<setw(8)<<gpu_time<<setprecision(6)<<" ms"<<endl;

	saveData(xc,fe,write_cells,flag);

	//free up GPU memory
	cudaFree(device_N);
	cudaFree(device_div);
	cudaFree(device_lb);
	cudaFree(device_ub);
	cudaFree(device_allcomb);
	cudaFree(device_Nnew);
}

int main(){
	int N[D] = {7,7,7,7,7};
	float lb[D] = {0,0,0,0,0};
	float ub[D] = {40,40,40,40,40};
	int div[D] = {1,1,1,1,1};
	int max_iter = 2;
	int max_period = 2;

	scm_mop(N,lb,ub,div,max_period,max_iter); //GPU implementation

	return 0;
}