#include "DTLZ.h"

/*DTLZ1-MOP*/
vector<double> DTLZ1(vector<double> xreal,vector<double> obj)
{
	int i,j;
	int n=xreal.size();
	int m=obj.size();
	int k=n-m+1;
	double f,g;

	g=0.0;
	for(i=n-k;i<n;i++)
	{
		g+=pow(xreal[i]-0.5,2.0)-cos(20.0*M_PI*(xreal[i]-0.5));
	}
	g=100.0*(k+g);

	for(i=1;i<=m;i++)
	{
		f=0.5*(1.0+g);
		for(j=m-i;j>=1;j--)
		{
			f*=xreal[j-1];
		}
		if(i>1)
		{
			f*=1-xreal[(m-i+1)-1];
		}
		obj[i-1]=f;
	}
	return obj;
}

/*DTLZ2-MOP*/
vector<double> DTLZ2(vector<double> xreal,vector<double> obj)
{
	int i,j;
	int n=xreal.size();
	int m=obj.size();
	int k=n-m+1;

	double f,g;

	g=0.0;
	for(i=n-k;i<n;i++)
	{
		g+=pow(xreal[i]-0.5,2.0);
	}

	for(i=1;i<=m;i++)
	{
		f=(1.0+g);
		for(j=m-i;j>=1;j--)
		{
			f*=cos(xreal[j-1]*M_PI/2.0);
		}
		if(i>1)
		{
			f*=sin(xreal[(m-i+1)-1]*M_PI/2.0);
		}

		obj[i-1]=f;
	}
	return obj;
}

/*DTLZ3-MOP*/
vector<double> DTLZ3(vector<double> xreal,vector<double> obj)
{
	int i,j;
	int n=xreal.size();
	int m=obj.size();
	int k=n-m+1;
	double f,g;

	g=0.0;
	for(i=n-k;i<n;i++)
	{
		g+=pow(xreal[i]-0.5,2.0)-cos(20.0*M_PI*(xreal[i]-0.5));
	}
	g=100.0*(k+g);

	for(i=1;i<=m;i++)
	{
		f=(1+g);
		for(j=m-i;j>=1;j--)
		{
			f*=cos(xreal[j-1]*M_PI/2.0);
		}
		if(i>1)
		{
			f*=sin(xreal[(m-i+1)-1]*M_PI/2.0);
		}

		obj[i-1]=f;
	}
	return obj;
}

/*DTLZ4-MOP*/
vector<double> DTLZ4(vector<double> xreal,vector<double> obj)
{
	int i,j;
	double alpha=100.0;
	int n=xreal.size();
	int m=obj.size();
	int k=n-m+1;
	double f,g;

	g=0.0;
	for(i=n-k;i<n;i++)
	{
		g+=pow(xreal[i]-0.5,2.0);
	}

	for(i=1;i<=m;i++)
	{
		f=(1+g);
		for(j=m-i;j>=1;j--)
		{
			f*=cos(pow(xreal[j-1],alpha)*M_PI/2.0);
		}
		if(i>1)
		{
			f*=sin(pow(xreal[(m-i+1)-1],alpha)*M_PI/2.0);
		}

		obj[i-1]=f;
	}
	return obj;
}

/*DTLZ5-MOP*/
vector<double> DTLZ5(vector<double> xreal,vector<double> obj)
{
	int i,j;
	int n=xreal.size();
	int m=obj.size();
	int k=n-m+1;
	double *theta,t;
	double f,g;

	g=0.0;
	for(i=n-k;i<n;i++)
	{
		g+=pow(xreal[i]-0.5,2.0);
	}

	theta=(double*)malloc(sizeof(double)*m);

	t=M_PI/(4.0*(1.0+g));
	theta[0]=xreal[0]*M_PI/2.0;

	for(i=1;i<m;i++)
	{
		theta[i]=t*(1.0+2.0*g*xreal[i]);
	}

	for(i=1;i<=m;i++)
	{
		f=(1+g);
		for(j=m-i;j>=1;j--)
		{
			f*=cos(theta[j-1]);
		}
		if(i>1)
		{
			f*=sin(theta[(m-i+1)-1]);
		}

		obj[i-1]=f;
	}

	free(theta);
	return obj;
}

/*DTLZ6-MOP*/
vector<double> DTLZ6(vector<double> xreal,vector<double> obj)
{
	int i,j;
	int n=xreal.size();
	int m=obj.size();
	int k=n-m+1;
	double *theta,t;
	double f,g;

	f=0.0;
	g=0.0;
	for(i=n-k;i<n;i++)
	{
		g+=pow(xreal[i],0.1);
	}

	theta=(double*)malloc(sizeof(double)*m);

	t=M_PI/(4.0*(1.0+g));
	theta[0]=xreal[0]*M_PI/2.0;
	for(i=1;i<m;i++)
	{
		theta[i]=t*(1.0+2.0*g*xreal[i]);
	}

	for(i=1;i<=m;i++)
	{
		f=(1+g);
		for(j=m-i;j>=1;j--)
		{
			f*=cos(theta[j-1]);
		}
		if(i>1)
		{
			f*=sin(theta[(m-i+1)-1]);
		}

		obj[i-1]=f;
	}

	free(theta);
	return obj;
}

/*DTLZ7-MOP*/
vector<double> DTLZ7(vector<double> xreal,vector<double> obj)
{
	int i,j;
	int n=xreal.size();
	int m=obj.size();
	int k=n-m+1;
	double g,h;

	g=0.0;
	for(i=n-k;i<n;i++)
	{
		g+=xreal[i];
	}
	g=1.0+9.0*g/(double)k;

	for(i=0;i<m-1;i++)
	{
		obj[i]=xreal[i];
	}

	h=0.0;
	for(j=0;j<m-1;j++)
	{
		h+=xreal[j]/(1.0+g)*(1.0+sin(3.0*M_PI*xreal[j]));
	}
	h=(double)m-h;

	obj[m-1]=(1.0+g)*h;
	return obj;
}
