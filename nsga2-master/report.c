/* Routines for storing population data into files */

# include <stdio.h>
# include <stdlib.h>
# include <math.h>

# include "nsga2.h"
# include "rand.h"

/* Function to print the information of a population in a file */
void report_pop (NSGA2Type *nsga2Params,  population *pop, FILE *fpt)
{
    
    int i, j;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        for (j=0; j<nsga2Params->nobj; j++)
        {
            fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
        }
        /*if (nsga2Params->ncon!=0)
        {
            for (j=0; j<nsga2Params->ncon; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
            }
        }
        if (nsga2Params->nreal!=0)
        {
            for (j=0; j<nsga2Params->nreal; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
            }
        }
        if (nsga2Params->nbin!=0)
        {
            for (j=0; j<nsga2Params->nbin; j++)
            {
                for (k=0; k<nsga2Params->nbits[j]; k++)
                {
                    fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                }
            }
        }*/
        //fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
        fprintf(fpt,"%d\t",pop->ind[i].rank);
        fprintf(fpt,"%d\n",pop->ind[i].valid);
        //fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);

    }
    return;
}


void report_pop_enco(NSGA2Type *nsga2Params,  population *pop, FILE *fpt)
{
       int i, j;
    for (j=0; j<nsga2Params->popsize; j++)
    {
        for(i=0; i<nsga2Params->encoding_len_ns; i++) {
            if (pop->ind[j].rank==1)
            {
                fprintf(fpt,"%i\t", pop->ind[j].absolute_encoding[i]);
            }
            
        }


        fprintf(fpt,"\n");

    }
    return;
}

void report_pop_coor(NSGA2Type *nsga2Params,  population *pop, FILE *fpt)
{
    int i, j,k;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        for (j=0; j<nsga2Params->sequence_len_ns; j++)
        {
            for(k=0; k<nsga2Params->dimensions_ns; k++) {
                fprintf(fpt,"%e\t",pop->ind[i].coordinates[j][k]);
            }
                fprintf(fpt,"\n");
        }
    }
    return;
}

/* Function to print the information of feasible and non-dominated population in a file */
void report_feasible (NSGA2Type *nsga2Params,  population *pop, FILE *fpt)
{
    int i, j;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        if (pop->ind[i].constr_violation == 0.0 && pop->ind[i].rank==1)
        {
            for (j=0; j<nsga2Params->nobj; j++)
            {
                fprintf(fpt,"%e\t",pop->ind[i].obj[j]);
            }
            /*if (nsga2Params->ncon!=0)
            {
                for (j=0; j<nsga2Params->ncon; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].constr[j]);
                }
            }
            if (nsga2Params->nreal!=0)
            {
                for (j=0; j<nsga2Params->nreal; j++)
                {
                    fprintf(fpt,"%e\t",pop->ind[i].xreal[j]);
                }
            }
            if (nsga2Params->nbin!=0)
            {
                for (j=0; j<nsga2Params->nbin; j++)
                {
                    for (k=0; k<nsga2Params->nbits[j]; k++)
                    {
                        fprintf(fpt,"%d\t",pop->ind[i].gene[j][k]);
                    }
                }
            }*/
            //fprintf(fpt,"%e\t",pop->ind[i].constr_violation);
            fprintf(fpt,"%d\n",pop->ind[i].rank);
            //fprintf(fpt,"%e\n",pop->ind[i].crowd_dist);
        }
    }
    return;
}

void report_pop_here (NSGA2Type *nsga2Params,  population *pop)
{
    printf("\n");
    int i,j;
    for (i=0; i<nsga2Params->popsize; i++)
    {
        for (j=0; j<nsga2Params->nobj; j++)
        {
            printf("%e\t",pop->ind[i].obj[j]);
        }
        /*
        if (nsga2Params->ncon!=0)
        {
            for (j=0; j<nsga2Params->ncon; j++)
            {
                printf("%e\t",pop->ind[i].constr[j]);
            }
        }
        if (nsga2Params->nreal!=0)
        {
            for (j=0; j<nsga2Params->nreal; j++)
            {
                printf("%e\t",pop->ind[i].xreal[j]);
            }
        }
        if (nsga2Params->nbin!=0)
        {
            for (j=0; j<nsga2Params->nbin; j++)
            {
                for (k=0; k<nsga2Params->nbits[j]; k++)
                {
                    printf("%d\t",pop->ind[i].gene[j][k]);
                }
            }
        }*/
        /*for (j=0; j<nsga2Params->nobj; j++)
        {
            printf("%e\t",pop->ind[i].objectives[j]);
        }*/
       

        //printf("%e\t",pop->ind[i].constr_violation);
        printf("%d\t",pop->ind[i].rank);
        printf("%e\t",pop->ind[i].crowd_dist);
        printf("%d\n",pop->ind[i].valid);


    }
    return;
}
