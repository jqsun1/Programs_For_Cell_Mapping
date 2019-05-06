#include <iostream>

using namespace std;

extern"C" {
void fortfunc_(int *ii, float *ff);
}

main()
{

   int ii=5;
   float ff=5.5;

   fortfunc_(&ii, &ff);

   return 0;
}
