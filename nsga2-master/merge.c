/* Routine for mergeing two populations */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Routine to merge two populations into one */
void merge(NSGA2Type *nsga2Params,  population *pop1, population *pop2, population *pop3)
{
    int i, k;

    for (i=0; i<nsga2Params->popsize; i++)
    {
        copy_ind (nsga2Params,  &(pop1->ind[i]), &(pop3->ind[i]));
    }

    for (i=0, k=nsga2Params->popsize; i<nsga2Params->popsize; i++, k++)
    {
        copy_ind (nsga2Params,  &(pop2->ind[i]), &(pop3->ind[k]));
    }
    return;
}

/* Routine to copy an individual 'ind1' into another individual 'ind2' */
void copy_ind (NSGA2Type *nsga2Params, individual *ind1, individual *ind2)
{

    int i,j;
    ind2->rank = ind1->rank;
    ind2->constr_violation = ind1->constr_violation;
    ind2->crowd_dist = ind1->crowd_dist;
    /*if (nsga2Params->nreal!=0)
    {
        for (i=0; i<nsga2Params->nreal; i++)
        {
            ind2->xreal[i] = ind1->xreal[i];
        }
    }*/
    /*if (nsga2Params->nbin!=0)
    {
        for (i=0; i<nsga2Params->nbin; i++)
        {
            ind2->xbin[i] = ind1->xbin[i];
            for (j=0; j<nsga2Params->nbits[i]; j++)
            {
                ind2->gene[i][j] = ind1->gene[i][j];
            }
        }
    }*/

    for (i=0; i<nsga2Params->nobj; i++)
    {
        ind2->obj[i] = ind1->obj[i];
    }
    if (nsga2Params->ncon!=0)
    {
        for (i=0; i<nsga2Params->ncon; i++)
        {
            ind2->constr[i] = ind1->constr[i];
        }
    }

    for(i=0; i<nsga2Params->encoding_len_ns; i++) 
        {
            
            ind2->absolute_encoding[i] = ind1->absolute_encoding[i];
        }
        


    for(int i=0; i<nsga2Params->sequence_len_ns; i++) {
           for(j=0; j<2; j++) {
            //printf("%f\t",ind1->coordinates[i][j]);
            ind2->coordinates[i][j] = ind1->coordinates[i][j];
            
        }
           
    }

    ind2->HHtc = ind1->HHtc;

    for(i=0; i<nsga2Params->nobj; i++) ind2->objectives[i] = ind1->objectives[i];  

    for(i=0; i<ind1->HHtc; i++){        
        ind2->contacts[i][0] = ind1->contacts[i][0];
        ind2->contacts[i][1] = ind1->contacts[i][1];
    }

    ind2->valid = ind1->valid;
   

    return;
}

