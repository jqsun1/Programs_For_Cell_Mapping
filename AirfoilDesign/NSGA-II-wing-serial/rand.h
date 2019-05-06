/* Declaration for random number related variables and routines */

# ifndef _RAND_H_
# define _RAND_H_

/* Variable declarations for the random number generator */
extern double seed;
extern double oldrand[55];
extern int jrand;

/* Function declarations for the random number generator */
void randomize(void);
void warmup_random (double seed);
void advance_random (void);
double randomperc(void);
int rnd (int low, int high);
double rndreal (double low, double high);

# endif /* _RAND_H_ */
