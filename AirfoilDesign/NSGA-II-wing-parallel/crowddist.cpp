/* Crowding distance computation routines */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "global.h"
# include "rand.h"
#include <iostream>
using namespace std;
/* Routine to compute crowding distance based on ojbective function values when the population in in the form of a list */
void assign_crowding_distance_list (population *pop, list *lst, int front_size)
{
	//cout << "assign_rank_and_crowding_dist is called" << endl;
    int **obj_array;
    int *dist;
    int i, j;
    list *temp;
    temp = lst;
    if (front_size==1)
    {
        pop->ind[lst->index].crowd_dist = INF;
        return;
    }
    if (front_size==2)
    {
        pop->ind[lst->index].crowd_dist = INF;
        pop->ind[lst->child->index].crowd_dist = INF;
        return;
    }
    //cout << "assign_rank_and_crowding_dist in the middle 1" << endl;
    obj_array = (int **)malloc(nobj*sizeof(int*));
    dist = (int *)malloc(front_size*sizeof(int));
    for (i=0; i<nobj; i++)
    {
        obj_array[i] = (int *)malloc(front_size*sizeof(int));
    }
    //cout << "assign_rank_and_crowding_dist in the middle 2" << endl;
    for (j=0; j<front_size; j++)
    {
        dist[j] = temp->index;
        temp = temp->child;
    }
    //cout << "assign_rank_and_crowding_dist in the middle 3" << endl;
    assign_crowding_distance (pop, dist, obj_array, front_size);
    //cout << "assign_rank_and_crowding_dist in the middle 4" << endl;
    free (dist);
    for (i=0; i<nobj; i++)
    {
        free (obj_array[i]);
    }
    free (obj_array);
    //cout << "assign_rank_and_crowding_dist is finished" << endl;
    return;
}

/* Routine to compute crowding distance based on objective function values when the population in in the form of an array */
void assign_crowding_distance_indices (population *pop, int c1, int c2)
{
    int **obj_array;
    int *dist;
    int i, j;
    int front_size;
    front_size = c2-c1+1;
    if (front_size==1)
    {
        pop->ind[c1].crowd_dist = INF;
        return;
    }
    if (front_size==2)
    {
        pop->ind[c1].crowd_dist = INF;
        pop->ind[c2].crowd_dist = INF;
        return;
    }
    obj_array = (int **)malloc(nobj*sizeof(int*));
    dist = (int *)malloc(front_size*sizeof(int));
    for (i=0; i<nobj; i++)
    {
        obj_array[i] = (int *)malloc(front_size*sizeof(int));
    }
    for (j=0; j<front_size; j++)
    {
        dist[j] = c1++;
    }
    assign_crowding_distance (pop, dist, obj_array, front_size);
    free (dist);
    for (i=0; i<nobj; i++)
    {
        free (obj_array[i]);
    }
    free (obj_array);
    return;
}

/* Routine to compute crowding distances */
void assign_crowding_distance (population *pop, int *dist, int **obj_array, int front_size)
{
	//cout << "assign_crowding_distance started" << endl;
    int i, j;
    for (i=0; i<nobj; i++)
    {
        for (j=0; j<front_size; j++)
        {
            obj_array[i][j] = dist[j];
        }
        //cout << "assign_crowding_distance middle 0.1" << endl;
        quicksort_front_obj (pop, i, obj_array[i], front_size);
        //cout << "assign_crowding_distance middle 0.2" << endl;
    }
    //cout << "assign_crowding_distance middle 1" << endl;
    for (j=0; j<front_size; j++)
    {
        pop->ind[dist[j]].crowd_dist = 0.0;
    }
    for (i=0; i<nobj; i++)
    {
        pop->ind[obj_array[i][0]].crowd_dist = INF;
    }
    //cout << "assign_crowding_distance middle 2" << endl;
    for (i=0; i<nobj; i++)
    {
        for (j=1; j<front_size-1; j++)
        {
            if (pop->ind[obj_array[i][j]].crowd_dist != INF)
            {
                if (pop->ind[obj_array[i][front_size-1]].obj[i] == pop->ind[obj_array[i][0]].obj[i])
                {
                    pop->ind[obj_array[i][j]].crowd_dist += 0.0;
                }
                else
                {
                    pop->ind[obj_array[i][j]].crowd_dist += (pop->ind[obj_array[i][j+1]].obj[i] - pop->ind[obj_array[i][j-1]].obj[i])/(pop->ind[obj_array[i][front_size-1]].obj[i] - pop->ind[obj_array[i][0]].obj[i]);
                }
            }
        }
    }
    //cout << "assign_crowding_distance middle 3" << endl;
    for (j=0; j<front_size; j++)
    {
        if (pop->ind[dist[j]].crowd_dist != INF)
        {
            pop->ind[dist[j]].crowd_dist = (pop->ind[dist[j]].crowd_dist)/nobj;
        }
    }
    //cout << "assign_crowding_distance finished" << endl;
    return;
}
