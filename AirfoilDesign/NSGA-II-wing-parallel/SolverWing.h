#ifndef solver_wing_h
#define solver_wing_h

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32) && !defined(__CYGWIN__)
	static const std::string slash="\\";
#else
	static const std::string slash="/";
#endif
#define pps 80 //points per side
#define inf 1e4 //infinity value
extern float coeff[];
extern float XX[];


extern"C" void xfoilfun_(bool * CONV, float * CL, float * CD, float * CM, float *XX, float *YY, float * Re, float * Alfa, long int * npointslong, long int * samples);

void matLU(float *A, float *L, float *U, int dim);
void matsolve(float *A, float *b, float *x, float* L, float* U, float *y, int n);
void parsec(float * a, float * param);
void parsec2pointsPrep(float * XX, float * coeff, int npoints);
void parsec2points(float * X, float * Y, float * a, int & npoints, float * XX, float * coeff);
bool fileRead(float * xtmp, float * ytmp, int & i, std::string & airfoil, std::string filename);
void f(float * y, float * x);

#endif /* solver_wing_h */
