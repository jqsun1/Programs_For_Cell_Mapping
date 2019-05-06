/* Definition of random number generation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"
#include <iostream>
using namespace std;
double seed;
double oldrand[55];
int jrand;

/* Get seed number for random and start it up */
void randomize()
{
    int j1;
    for(j1=0; j1<=54; j1++)
    {
        oldrand[j1] = 0.0;
    }
    jrand=0;
    warmup_random (seed);
    return;
}

/* Get randomize off and running */
void warmup_random (double seed)
{
	//cout << "warmup_random started" << endl;
    int j1, ii;
    double new_random, prev_random;
    oldrand[54] = seed;
    new_random = 0.000000001;
    prev_random = seed;
    for(j1=1; j1<=54; j1++)
    {
        ii = (21*j1)%54;
        oldrand[ii] = new_random;
        new_random = prev_random-new_random;
        if(new_random<0.0)
        {
            new_random += 1.0;
        }
        prev_random = oldrand[ii];
    }
    advance_random ();
    advance_random ();
    advance_random ();
    jrand = 0;
    return;
}

/* Create next batch of 55 random numbers */
void advance_random ()
{
    int j1;
    double new_random;
    for(j1=0; j1<24; j1++)
    {
        new_random = oldrand[j1]-oldrand[j1+31];
        if(new_random<0.0)
        {
            new_random = new_random+1.0;
        }
        //cout << "advance_random, j="<< j1 <<", rand= " << new_random<< endl;
        oldrand[j1] = new_random;
    }
    for(j1=24; j1<55; j1++)
    {
        new_random = oldrand[j1]-oldrand[j1-24];
        if(new_random<0.0)
        {
            new_random = new_random+1.0;
        }
        //cout << "advance_random, j="<< j1 <<", rand= " << new_random<< endl;
        oldrand[j1] = new_random;
    }
}

/* Fetch a single random number between 0.0 and 1.0 */
/*
double randomperc()
{
    jrand++;
    if(jrand>=55)
    {
        jrand = 1;
        advance_random();
    }
    return((double)oldrand[jrand]);
}
*/
/* Fetch a single random integer between low and high including the bounds */
/*
int rnd (int low, int high)
{
	cout << "rnd started" << endl;
    int res;
    if (low >= high)
    {
        res = low;
    }
    else
    {
        //res = low + (randomperc()*(high-low+1));
        cout << "randomperc()= " <<randomperc() << endl;
        res = (int) ((double)low + (randomperc()*(double)(high-low+1)));
        if (res > high)
        {
            res = high;
        }
    }
    return (res);
}
*/
/* Fetch a single random real number between low and high including the bounds */
/*
double rndreal (double low, double high)
{
    return (low + (high-low)*randomperc());
}
*/
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
 /* Random number generator */
int rnd(int min, int max) {
	float out;
	out = ((float)rand()/(float)RAND_MAX*((float)max-(float)min + 1))+(float)min;
	return (int)out;
}
////////////////////////////////////////////////////////////////////////////////
 /* Random number generator */
double rndreal (double min, double max) {
	double out;
	out = ((double)rand()/(double)RAND_MAX*((double)max-(double)min))+(double)min;
	return out;
}
////////////////////////////////////////////////////////////////////////////////
 /* Random number generator */
double randomperc () {
	double out;
	out = (double)rand()/(double)RAND_MAX;
	return out;
}
