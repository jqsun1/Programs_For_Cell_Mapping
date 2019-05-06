#ifndef solver_wing_h
#define solver_wing_h

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	static const std::string slash="\\";
#else
	static const std::string slash="/";
#endif
#define pps 80 //points per side
#define inf 1e4 //infinity value
// MOP setup
#define Dim 12 //Number of dimensions
#define O 3 //Number of objective functions
#define maxArrLen (long long int)  5000000 //Max length of arrays in calculations
#define maxTempLen (long long int) 5000000 //Max length of arrays in calculations
//   		 Rleu   Rlel   Xup   Yup   Yxxup  Xlo   Ylo    Yxxlo Tte  Toff Ate Bte

extern double coeff[];
extern double XX[];


//extern"C" void xfoilfun_(bool * CONV, double * CL, double * CD, double * CM, double *XX, double *YY, double * Re, double * Alfa, long int * npointslong, long int * samples);
//extern"C" {
	// everything is passed by pointer
//	void xfoilfun_(bool * CONV, float * CL, float * CD, float * CM, float *XX, float *YY, float * Re, float * Alfa, long int * npointslong, long int * samples);
//}


void matLU(double *A, double *L, double *U, int dim);
void matsolve(double *A, double *b, double *x, double* L, double* U, double *y, int n);
void parsec(double * a, double * param);
void parsec2pointsPrep(double * XX, double * coeff, int npoints);
void parsec2points(double * X, double * Y, double * a, int & npoints, double * XX, double * coeff);
bool fileRead(double * xtmp, double * ytmp, int & i, std::string & airfoil, std::string filename);
void f(double * y, double * x);

#endif /* solver_wing_h */
