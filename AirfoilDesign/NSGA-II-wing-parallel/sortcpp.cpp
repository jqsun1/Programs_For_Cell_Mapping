/* Routines for randomized recursive quick-sort */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"
#include <iostream>
using namespace std;
/* Randomized quick sort routine to sort a population based on a particular objective chosen */
void quicksort_front_obj(population *pop, int objcount, int obj_array[], int obj_array_size)
{
    q_sort_front_obj (pop, objcount, obj_array, 0, obj_array_size-1);
    return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on a particular objective chosen */
void q_sort_front_obj(population *pop, int objcount, int obj_array[], int left, int right)
{
	//cout << "q_sort_front_obj started" << endl;
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right)
    {
    	//cout << "q_sort_front_obj middle 1" << endl;
        index = rnd (left, right);
        //cout << "q_sort_front_obj middle 1.1, obj_array[index]=" <<  obj_array[right]<< endl;
        //cout << "right=" << right << endl;
        //cout << "left=" << left << endl;
        //cout << "index=" << index << endl;
        temp = obj_array[right];
        obj_array[right] = obj_array[index];
        obj_array[index] = temp;
        //cout << "q_sort_front_obj middle 1.2" << endl;
        pivot = pop->ind[obj_array[right]].obj[objcount];
        //cout << "q_sort_front_obj middle 1.3" << endl;
        i = left-1;
        //cout << "q_sort_front_obj middle 2" << endl;
        for (j=left; j<right; j++)
        {
            if (pop->ind[obj_array[j]].obj[objcount] <= pivot)
            {
                i+=1;
                temp = obj_array[j];
                obj_array[j] = obj_array[i];
                obj_array[i] = temp;
            }
        }
        //cout << "q_sort_front_obj middle 3" << endl;
        index=i+1;
        temp = obj_array[index];
        obj_array[index] = obj_array[right];
        obj_array[right] = temp;
        //cout << "q_sort_front_obj middle 4" << endl;
        q_sort_front_obj (pop, objcount, obj_array, left, index-1);
        q_sort_front_obj (pop, objcount, obj_array, index+1, right);
        //cout << "q_sort_front_obj middle 5" << endl;
    }
    //cout << "q_sort_front_obj finished" << endl;
    return;
}

/* Randomized quick sort routine to sort a population based on crowding distance */
void quicksort_dist(population *pop, int *dist, int front_size)
{
    q_sort_dist (pop, dist, 0, front_size-1);
    return;
}

/* Actual implementation of the randomized quick sort used to sort a population based on crowding distance */
void q_sort_dist(population *pop, int *dist, int left, int right)
{
    int index;
    int temp;
    int i, j;
    double pivot;
    if (left<right)
    {
        index = rnd (left, right);
        temp = dist[right];
        dist[right] = dist[index];
        dist[index] = temp;
        pivot = pop->ind[dist[right]].crowd_dist;
        i = left-1;
        for (j=left; j<right; j++)
        {
            if (pop->ind[dist[j]].crowd_dist <= pivot)
            {
                i+=1;
                temp = dist[j];
                dist[j] = dist[i];
                dist[i] = temp;
            }
        }
        index=i+1;
        temp = dist[index];
        dist[index] = dist[right];
        dist[right] = temp;
        q_sort_dist (pop, dist, left, index-1);
        q_sort_dist (pop, dist, index+1, right);
    }
    return;
}
