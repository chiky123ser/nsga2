/* NSGA-II routine (implementation of the 'main' function) */

# include "nsga2.h"
# include "rand.h"





#include "../PSPlib/psp.h"

//FILE *fpt1;
FILE *fpt2;
FILE *fpt3;
//FILE *fpt4;
//FILE *fpt5;
//FILE *fpt6;
FILE *gp;
population *parent_pop;
population *child_pop;
population *mixed_pop;
char* No_archivo;
char* No_problema;
int type_Crossover;



NSGA2Type ReadParameters(int argc, char **argv){
    //int i;
    NSGA2Type nsga2Params;
    
    if (argc<2)
    {
        printf("\n Usage ./nsga2r random_seed \n");
        exit(1);
    }
    nsga2Params.seed = (double)atof(argv[1]);
    if (nsga2Params.seed<=0.0 || nsga2Params.seed>=1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(1);
    }

     
    //printf("\n\tUSAGE:    LATTICE_FILE   PROBLEM_INSTANCE   FORMULATION={SO, MO}   EVALUATION={D85, K99..., DEC1, DEC2...}   VARIANT = {DET, RND, DYN, DYN10..30}   MAX_EVAL   RUN\n\n"), exit(1);   
    
    // Load lattice configuration
    load_lattice_configuration(argv[2]);
    
    // Read problem instance
    read_problem_instance(argv[3]);
    
    // Problem formulation to be used
    strcpy(problem_formulation, argv[4]);
    
    // Evaluation function or multiobjectivization strategy
    strcpy(evaluation_strategy, argv[5]);   
    
    // Type of subsets to use: DET, RND, DYN
    strcpy(strategy_variant, argv[6]);

    // Allowed number of energy function evaluations
    //max_evaluations = atoi(argv[7]);    
    
    // Run and seed for random numbers
    /*srand(time(NULL));
    int r = 1+(rand()%9);*/

    run = atoi(argv[7]);
    seed = run*sequence_len*dimensions;
    initialize_rnd(seed);   
        
    // Parameters Adjustment
    formulation_adjustment();

    if(n_objectives>1){
        Hsubset = (int *)malloc(sizeof(int) * sequence_len);
        compute_subsets();  
    }

    nsga2Params.sequence_len_ns=sequence_len;               // Length of the sequence 
    nsga2Params.optimalHHtc_ns=optimalHHtc;
    nsga2Params.encoding_len_ns=encoding_len;
    nsga2Params.dimensions_ns=dimensions;

    //printf("encoding_len %i\n",encoding_len );
    //printf("sequence_len %i\n", sequence_len);
 
    /*printf("\n Enter the problem relevant and algorithm relevant parameters ... ");
    printf("\n Enter the population size (a multiple of 4) : ");*/
    nsga2Params.popsize = (int)atof(argv[8]);
    //scanf("%d",&nsga2Params.popsize);

    if (nsga2Params.popsize<4 || (nsga2Params.popsize%4)!= 0)
    {
        printf("\n population size read is : %d",nsga2Params.popsize);
        printf("\n Wrong population size entered, hence exiting \n");
        exit (1);
    }
    //printf("\n Enter the number of generations : ");
    nsga2Params.ngen = (int)atof(argv[9]);
    //scanf("%d",&nsga2Params.ngen);
    if (nsga2Params.ngen<1)
    {
        printf("\n number of generations read is : %d",nsga2Params.ngen);
        printf("\n Wrong nuber of generations entered, hence exiting \n");
        exit (1);
    }
    //printf("\n Enter the number of objectives : ");
    nsga2Params.nobj = (int)atof(argv[10]);
    //scanf("%d",&nsga2Params.nobj);
    if (nsga2Params.nobj<1)
    {
        printf("\n number of objectives entered is : %d",nsga2Params.nobj);
        printf("\n Wrong number of objectives entered, hence exiting \n");
        exit (1);
    }
   /* printf("\n Enter the number of constraints : ");
    scanf("%d",&nsga2Params.ncon);
    if (nsga2Params.ncon<0)
    {
        printf("\n number of constraints entered is : %d",nsga2Params.ncon);
        printf("\n Wrong number of constraints enetered, hence exiting \n");
        exit (1);
    }*/
    //printf("\n Enter the number of objective function HP : ");
    nsga2Params.nreal = (int)atof(argv[11]);
    //scanf("%d",&nsga2Params.nreal);
    if (nsga2Params.nreal<0)
    {
        printf("\n number of objective function HP entered is : %d",nsga2Params.nreal);
        printf("\n Wrong number of variables entered, hence exiting \n");
        exit (1);
    }
    if (nsga2Params.nreal != 0)
    {
        nsga2Params.min_realvar = (double *)malloc(nsga2Params.nreal*sizeof(double));
        nsga2Params.max_realvar = (double *)malloc(nsga2Params.nreal*sizeof(double));
        /*for (i=0; i<nsga2Params.nreal; i++)
        {
            printf ("\n Enter the lower limit of real variable %d : ",i+1);
            scanf ("%lf",&nsga2Params.min_realvar[i]);
            printf ("\n Enter the upper limit of real variable %d : ",i+1);
            scanf ("%lf",&nsga2Params.max_realvar[i]);
            if (nsga2Params.max_realvar[i] <= nsga2Params.min_realvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of real variable, hence exiting \n");
                exit(1);
            }
        }*/
        //printf ("\n Enter the probability of crossover  (0.6-1.0) : ");
        nsga2Params.pcross_real = (double)atof(argv[12]);
        //scanf ("%lf",&nsga2Params.pcross_real);
        if (nsga2Params.pcross_real<0.0 || nsga2Params.pcross_real>1.0)
        {
            printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_real);
            printf("\n Entered value of probability of crossover of objective is out of bounds, hence exiting \n");
            exit (1);
        }
        //printf ("\n Enter the probablity of mutation  (1/nsga2Params.nreal) : ");
        nsga2Params.pmut_real = (double)atof(argv[13]);
        //scanf ("%lf",&nsga2Params.pmut_real);
        if (nsga2Params.pmut_real<0.0 || nsga2Params.pmut_real>1.0)
        {
            printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_real);
            printf("\n Entered value of probability of mutation of real variables is out of bounds, hence exiting \n");
            exit (1);
        }
        /*printf ("\n Enter the value of distribution index for crossover (5-20): ");
        scanf ("%lf",&nsga2Params.eta_c);
        if (nsga2Params.eta_c<=0)
        {
            printf("\n The value entered is : %e",nsga2Params.eta_c);
            printf("\n Wrong value of distribution index for crossover entered, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the value of distribution index for mutation (5-50): ");
        scanf ("%lf",&nsga2Params.eta_m);
        if (nsga2Params.eta_m<=0)
        {
            printf("\n The value entered is : %e",nsga2Params.eta_m);
            printf("\n Wrong value of distribution index for mutation entered, hence exiting \n");
            exit (1);
        }*/
    }
    /*printf("\n Enter the number of binary variables : ");
    scanf("%d",&nsga2Params.nbin);
    if (nsga2Params.nbin<0)
    {
        printf ("\n number of binary variables entered is : %d",nsga2Params.nbin);
        printf ("\n Wrong number of binary variables entered, hence exiting \n");
        exit(1);
    }
    if (nsga2Params.nbin != 0)
    {
        nsga2Params.nbits = (int *)malloc(nsga2Params.nbin*sizeof(int));
        nsga2Params.min_binvar = (double *)malloc(nsga2Params.nbin*sizeof(double));
        nsga2Params.max_binvar = (double *)malloc(nsga2Params.nbin*sizeof(double));
        for (i=0; i<nsga2Params.nbin; i++)
        {
            printf ("\n Enter the number of bits for binary variable %d : ",i+1);
            scanf ("%d",&nsga2Params.nbits[i]);
            if (nsga2Params.nbits[i] < 1)
            {
                printf("\n Wrong number of bits for binary variable entered, hence exiting");
                exit(1);
            }
            printf ("\n Enter the lower limit of binary variable %d : ",i+1);
            scanf ("%lf",&nsga2Params.min_binvar[i]);
            printf ("\n Enter the upper limit of binary variable %d : ",i+1);
            scanf ("%lf",&nsga2Params.max_binvar[i]);
            if (nsga2Params.max_binvar[i] <= nsga2Params.min_binvar[i])
            {
                printf("\n Wrong limits entered for the min and max bounds of binary variable entered, hence exiting \n");
                exit(1);
            }
        }
        printf ("\n Enter the probability of crossover of binary variable (0.6-1.0): ");
        scanf ("%lf",&nsga2Params.pcross_bin);
        if (nsga2Params.pcross_bin<0.0 || nsga2Params.pcross_bin>1.0)
        {
            printf("\n Probability of crossover entered is : %e",nsga2Params.pcross_bin);
            printf("\n Entered value of probability of crossover of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
        printf ("\n Enter the probability of mutation of binary variables (1/nsga2Params.nbits): ");
        scanf ("%lf",&nsga2Params.pmut_bin);
        if (nsga2Params.pmut_bin<0.0 || nsga2Params.pmut_bin>1.0)
        {
            printf("\n Probability of mutation entered is : %e",nsga2Params.pmut_bin);
            printf("\n Entered value of probability  of mutation of binary variables is out of bounds, hence exiting \n");
            exit (1);
        }
    }*/
    if (nsga2Params.nreal==0 && nsga2Params.nbin==0)
    {
        printf("\n Number of real as well as binary variables, both are zero, hence exiting \n");
        exit(1);
    }
    nsga2Params.choice=0;
    /*printf("\n Do you want to use gnuplot to display the results realtime (0 for NO) (1 for yes) : ");
    scanf("%d",&nsga2Params.choice);
    if (nsga2Params.choice!=0 && nsga2Params.choice!=1)
    {
        printf("\n Entered the wrong nsga2Params.choice, hence exiting, nsga2Params.choice entered was %d\n",nsga2Params.choice);
        exit(1);
    }
    if (nsga2Params.choice==1)
    {
        gp = popen(GNUPLOT_COMMAND,"w");
        if (gp==NULL)
        {
            printf("\n Could not open a pipe to gnuplot, check the definition of GNUPLOT_COMMAND in file global.h\n");
            printf("\n Edit the string to suit your system configuration and rerun the program\n");
            exit(1);
        }
        if (nsga2Params.nobj==2)
        {
            printf("\n Enter the objective for X axis display : ");
            scanf("%d",&nsga2Params.obj1);
            if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
            {
                printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
                exit(1);
            }
            printf("\n Enter the objective for Y axis display : ");
            scanf("%d",&nsga2Params.obj2);
            if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
            {
                printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
                exit(1);
            }
            nsga2Params.obj3 = -1;
        }
        else
        {
            printf("\n #obj > 2, 2D display or a 3D display ?, enter 2 for 2D and 3 for 3D :");
            scanf("%d",&nsga2Params.choice);
            if (nsga2Params.choice!=2 && nsga2Params.choice!=3)
            {
                printf("\n Entered the wrong nsga2Params.choice, hence exiting, nsga2Params.choice entered was %d\n",nsga2Params.choice);
                exit(1);
            }
            if (nsga2Params.choice==2)
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d",&nsga2Params.obj1);
                if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                scanf("%d",&nsga2Params.obj2);
                if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
                    exit(1);
                }
                nsga2Params.obj3 = -1;
            }
            else
            {
                printf("\n Enter the objective for X axis display : ");
                scanf("%d",&nsga2Params.obj1);
                if (nsga2Params.obj1<1 || nsga2Params.obj1>nsga2Params.nobj)
                {
                    printf("\n Wrong value of X objective entered, value entered was %d\n",nsga2Params.obj1);
                    exit(1);
                }
                printf("\n Enter the objective for Y axis display : ");
                scanf("%d",&nsga2Params.obj2);
                if (nsga2Params.obj2<1 || nsga2Params.obj2>nsga2Params.nobj)
                {
                    printf("\n Wrong value of Y objective entered, value entered was %d\n",nsga2Params.obj2);
                    exit(1);
                }
                printf("\n Enter the objective for Z axis display : ");
                scanf("%d",&nsga2Params.obj3);
                if (nsga2Params.obj3<1 || nsga2Params.obj3>nsga2Params.nobj)
                {
                    printf("\n Wrong value of Z objective entered, value entered was %d\n",nsga2Params.obj3);
                    exit(1);
                }
                printf("\n You have chosen 3D display, hence location of eye required \n");
                printf("\n Enter the first angle (an integer in the range 0-180) (if not known, enter 60) :");
                scanf("%d",&nsga2Params.angle1);
                if (nsga2Params.angle1<0 || nsga2Params.angle1>180)
                {
                    printf("\n Wrong value for first angle entered, hence exiting \n");
                    exit(1);
                }
                printf("\n Enter the second angle (an integer in the range 0-360) (if not known, enter 30) :");
                scanf("%d",&nsga2Params.angle2);
                if (nsga2Params.angle2<0 || nsga2Params.angle2>360)
                {
                    printf("\n Wrong value for second angle entered, hence exiting \n");
                    exit(1);
                }
            }
        }
    }*/
    No_problema=argv[14];
    No_archivo=argv[15];
    type_Crossover = atoi(argv[16]);
    //printf("%i\n", type_Crossover);


    return nsga2Params;
}



void InitNSGA2(NSGA2Type *nsga2Params, void *inp, void *out)
{
    //int i;


    char *cruce_promblema=cross(type_Crossover);
     
    //cruce_promblema=cross(type_Crossover);
    //printf("%s\n", cruce_promblema);

    //printf("\n == InitNSGA2 ==");
    
    // Initialize the files...
    //char initial[]="/home/servando/Documents/Frankenstein/Resu/output/a1_initial_pop.out";
    
    char final[150]="/home/servando/Documents/Frankenstein/Resu/";
    char best[150]="/home/servando/Documents/Frankenstein/Resu/";
    //char final[150]="Dec/Cycle/final/a_final_pop.out";
    //char best[150]="Dec/Cycle/best/a_best_pop.out";

    //char all[100]="a_all_pop.out";
    //char params[]="/home/servando/Documents/Frankenstein/Resu/output/a1_params.out";
    //char coor[100]="/home/servando/Documents/Frankenstein/Resu/output/a_coor_pop.out";
    char t[]="_";
    strcat( final, evaluation_strategy );
    strcat( best, evaluation_strategy );

   strcat( final, cruce_promblema );
   strcat( best, cruce_promblema );

    strcat( final, "/final/a_final_pop.out" );
    strcat( best, "/best/a_best_pop.out" );
    
    //strcat( initial, No_problema );
    strcat( final, No_problema );
    strcat( best, No_problema );
    //strcat( all, No_problema );
    //strcat( params, No_problema );
    //strcat( coor, No_problema );

    strcat( final, t );
    strcat( best, t );
    //strcat( all, t );
    //strcat( coor, t );
     
    //strcat( initial, No_archivo );
    strcat( final, No_archivo );
    strcat( best, No_archivo );

    //printf("%s\n", best);
    //strcat( all, No_archivo );
    //strcat( params, No_archivo );
    //strcat( coor, No_archivo );

    //fpt1 = fopen(initial,"w");
    fpt2 = fopen(final,"w");
    fpt3 = fopen(best,"w"); 
    //fpt4 = fopen(all,"w");
    //fpt5 = fopen(params,"w");
    //fpt6 = fopen(coor,"w");
    //fprintf(fpt1,"#\tThis file contains the data of initial population\n");
    fprintf(fpt2,"#\tThis file contains the data of final population\n");
    fprintf(fpt3,"#\tThis file contains the data of final feasible population (if found)\n");
    /*fprintf(fpt4,"#\tThis file contains the data of all generations\n");
    fprintf(fpt5,"#\tThis file contains information about inputs as read by the program\n");*/
    //fprintf(fpt6,"#\tThis file contains the coordinates of all generations\n");/*
    /*fprintf(fpt5,"\n Population size = %d",nsga2Params->popsize);
    fprintf(fpt5,"\n Number of generations = %d",nsga2Params->ngen);
    fprintf(fpt5,"\n Number of objective functions = %d",nsga2Params->nobj);
    //fprintf(fpt5,"\n Number of constraints = %d",nsga2Params->ncon);
    fprintf(fpt5,"\n Number of real variables = %d",nsga2Params->nreal);*/
    /*if (nsga2Params->nreal!=0)
    {
        for (i=0; i<nsga2Params->nreal; i++)
        {
            fprintf(fpt5,"\n Lower limit of real variable %d = %e",i+1,nsga2Params->min_realvar[i]);
            fprintf(fpt5,"\n Upper limit of real variable %d = %e",i+1,nsga2Params->max_realvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of real variable = %e",nsga2Params->pcross_real);
        fprintf(fpt5,"\n Probability of mutation of real variable = %e",nsga2Params->pmut_real);
        fprintf(fpt5,"\n Distribution index for crossover = %e",nsga2Params->eta_c);
        fprintf(fpt5,"\n Distribution index for mutation = %e",nsga2Params->eta_m);
    }*/

    /*fprintf(fpt5,"\n Number of binary variables = %d",nsga2Params->nbin);
    if (nsga2Params->nbin!=0)
    {
        for (i=0; i<nsga2Params->nbin; i++)
        {
            fprintf(fpt5,"\n Number of bits for binary variable %d = %d",i+1,nsga2Params->nbits[i]);
            fprintf(fpt5,"\n Lower limit of binary variable %d = %e",i+1,nsga2Params->min_binvar[i]);
            fprintf(fpt5,"\n Upper limit of binary variable %d = %e",i+1,nsga2Params->max_binvar[i]);
        }
        fprintf(fpt5,"\n Probability of crossover of binary variable = %e",nsga2Params->pcross_bin);
        fprintf(fpt5,"\n Probability of mutation of binary variable = %e",nsga2Params->pmut_bin);
    }*/
    //fprintf(fpt5,"\n Seed for random number generator = %e",nsga2Params->seed);
    /*nsga2Params->bitlength = 0;
    if (nsga2Params->nbin!=0)
    {
        for (i=0; i<nsga2Params->nbin; i++)
        {
            nsga2Params->bitlength += nsga2Params->nbits[i];
        }
    }*/
    //fprintf(fpt1,"#\tof objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params->nobj,nsga2Params->ncon,nsga2Params->nreal,nsga2Params->bitlength);
    //fprintf(fpt2,"#\tof objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params->nobj,nsga2Params->ncon,nsga2Params->nreal,nsga2Params->bitlength);
    //fprintf(fpt3,"#\tof objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params->nobj,nsga2Params->ncon,nsga2Params->nreal,nsga2Params->bitlength);
    //fprintf(fpt4,"#\tof objectives = %d, # of constraints = %d, # of real_var = %d, # of bits of bin_var = %d, constr_violation, rank, crowding_distance\n",nsga2Params->nobj,nsga2Params->ncon,nsga2Params->nreal,nsga2Params->bitlength);
    
    nsga2Params->nbinmut = 0;
    nsga2Params->nrealmut = 0;
    nsga2Params->nbincross = 0;
    nsga2Params->nrealcross = 0;
    
    // Initializing the populations
    parent_pop = (population *)malloc(sizeof(population));
    child_pop = (population *)malloc(sizeof(population));
    mixed_pop = (population *)malloc(sizeof(population));

    allocate_memory_pop (nsga2Params, parent_pop, nsga2Params->popsize);

    allocate_memory_pop (nsga2Params, child_pop, nsga2Params->popsize);
    allocate_memory_pop (nsga2Params, mixed_pop, 2*nsga2Params->popsize);
    
    // Preparing first Population
    randomize(nsga2Params->seed);
   
initialize_pop (nsga2Params,  parent_pop);
    
    /*printf("\n Initialization done, now performing first generation");
    printf("\n");*/
    //decode_pop(nsga2Params, parent_pop);
    //evaluate_pop (nsga2Params, parent_pop, inp, out);
assign_rank_and_crowding_distance (nsga2Params, parent_pop);

    
    //printf("///////////////////////\n");
    //print_pop(nsga2Params,parent_pop,nsga2Params->popsize);  

    //report_pop (nsga2Params, parent_pop, fpt1);
    //fprintf(fpt4,"#\tgen = 1\n");
    //fprintf(fpt6,"#\tgen = 1\n");
    //report_pop_coor(nsga2Params, parent_pop,fpt6);
    //report_pop(nsga2Params, parent_pop,fpt4);
    //report_pop_here(nsga2Params, parent_pop);
    //printf("\n -- Generation 1 --");
    //printf("\n");
    
    fflush(stdout);
    if (nsga2Params->choice!=0)    onthefly_display (nsga2Params, parent_pop,gp,1);
    //fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    //fflush(fpt4);
    /*fflush(fpt5);*/
    //fflush(fpt6);
    
}




int NSGA2(NSGA2Type *nsga2Params, void *inp, void *out)
{   
    int i;
    
    for (i=2; i<=nsga2Params->ngen; i++)
    {
          //printf("in Selection\n");
       
        
        selection (nsga2Params,  parent_pop, child_pop);
        
        //printf("out Selection\n");
        //report_pop_here(nsga2Params, child_pop);
        //printf("in\n");
        mutation_pop (nsga2Params,  child_pop);
      // printf("out\n");
        //printf("//////////////////Mutacion///////////////\n");
        //report_pop_here(nsga2Params, child_pop);    
        //evaluate_pop(nsga2Params,  child_pop, inp, out);    
        //assign_rank_and_crowding_distance (nsga2Params, child_pop);
        //print_pop(nsga2Params,child_pop,nsga2Params->popsize);  
        //decode_pop(nsga2Params,  child_pop);
        //evaluate_pop(nsga2Params,  child_pop, inp, out);
        
        merge (nsga2Params,  parent_pop, child_pop, mixed_pop);
        //printf("///////////////merge///////////////////////\n");
        /*nsga2Params->popsize = nsga2Params->popsize*2;
        report_pop_here(nsga2Params, mixed_pop);*/
        //print_pop(nsga2Params,mixed_pop,nsga2Params->popsize*2);

        fill_nondominated_sort (nsga2Params,  mixed_pop, parent_pop);

//report_pop_here(nsga2Params, parent_pop);
        //printf("////////////new pop//////////////\n");
        //report_pop_here(nsga2Params, parent_pop);
        //print_pop(nsga2Params,parent_pop,nsga2Params->popsize); 
        //printf("\n");
        //report_pop_here(nsga2Params, parent_pop);
    //run .5 /home/servando/Documents/Multiobjetivizacion/Lattices/2D_Square.lat /home/servando/Documents/Multiobjetivizacion/Instances/2D_Square/2d4.hp MO DEC2 DYN30 2 200 3 2 2 .9 .5 4 2

        /* Comment following four lines if information for all
         generations is not desired, it will speed up the execution */
        //fprintf(fpt4,"#\tgen = %d\n",i);
        //report_pop_coor(nsga2Params, parent_pop,fpt6);
        
        
        //fflush(fpt4);
       // fflush(fpt6);
        if (nsga2Params->choice!=0)    onthefly_display (nsga2Params, parent_pop,gp,i);

        //printf("\n -- Generation %d --", i);
        //print_pop(nsga2Params,parent_pop,nsga2Params->popsize);
    }
printf("\n Generations finished, now reporting solutions\n");

    report_pop(nsga2Params,  parent_pop,fpt2);
    report_feasible(nsga2Params,  parent_pop,fpt3);
    


    /*if (nsga2Params->nreal!=0)
    {
        fprintf(fpt5,"\n Number of crossover of real variable = %d",nsga2Params->nrealcross);
        fprintf(fpt5,"\n Number of mutation of real variable = %d",nsga2Params->nrealmut);
    }
    if (nsga2Params->nbin!=0)
    {
        fprintf(fpt5,"\n Number of crossover of binary variable = %d",nsga2Params->nbincross);
        fprintf(fpt5,"\n Number of mutation of binary variable = %d",nsga2Params->nbinmut);
    }*/

    
    // Closing the files and freeing up memories...
    fflush(stdout);
    //fflush(fpt1);
    fflush(fpt2);
    fflush(fpt3);
    //fflush(fpt4);
    //fflush(fpt5);
    //fflush(fpt6);
    //fclose(fpt1);
    fclose(fpt2);
    fclose(fpt3);
    //fclose(fpt4);
    //fclose(fpt5);
    //fclose(fpt6);

    if (nsga2Params->choice!=0)
    {
        pclose(gp);
    }
    if (nsga2Params->nreal!=0)
    {
        free (nsga2Params->min_realvar);
        free (nsga2Params->max_realvar);
    }
    if (nsga2Params->nbin!=0)
    {
        free (nsga2Params->min_binvar);
        free (nsga2Params->max_binvar);
        free (nsga2Params->nbits);
    }

   
    deallocate_memory_pop (nsga2Params,  parent_pop, nsga2Params->popsize);
   
    deallocate_memory_pop (nsga2Params,  child_pop, nsga2Params->popsize);

    deallocate_memory_pop (nsga2Params,  mixed_pop, 2*nsga2Params->popsize);


    free (parent_pop);
    free (child_pop);
    free (mixed_pop);

    //printf("\n Routine successfully exited \n");
    return (0);
}


void print_nsga2Params(NSGA2Type *nsga2Params){
    int i;
    
    printf("NSGA2 Parameters:\n");
    printf("\n seed number is : %lf", nsga2Params->seed);
    printf("\n population size read is : %d",nsga2Params->popsize);
    printf("\n number of generations read is : %d",nsga2Params->ngen);
    printf("\n number of objectives entered is : %d",nsga2Params->nobj);
    printf("\n number of constraints entered is : %d",nsga2Params->ncon);
    printf("\n number of real variables entered is : %d",nsga2Params->nreal);
    printf("\n variables bounds: ");
    for (i=0; i<nsga2Params->nreal; i++){
        printf("[%lf", nsga2Params->min_realvar[i]);
        printf(" %lf], ", nsga2Params->max_realvar[i]);
    }
    printf("\n Probability of crossover entered is : %e",nsga2Params->pcross_real);
    printf("\n Probability of mutation entered is : %e",nsga2Params->pmut_real);
    printf("\n The eta_c entered is : %e",nsga2Params->eta_c);
    printf("\n The eta_m entered is : %e",nsga2Params->eta_m);
    if (nsga2Params->nbin != 0){
        printf ("\n number of binary variables entered is : %d",nsga2Params->nbin);
        printf("\n Probability of crossover entered is : %e",nsga2Params->pcross_bin);
        printf("\n Probability of mutation entered is : %e",nsga2Params->pmut_bin);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void New_Ind(individual *ind){

    struct conformation_def current_solution;
    // Allocating memory
    allocate_memory_conformation(&current_solution);    
     // Generate initial solution at random     
    generate_valid_random_conformation(&current_solution);  
    //Evaluate
    evaluate_conformation(&current_solution);   
    /*printf("current\n");
    print_conformation(&current_solution);
    printf("end current_solution\n");
    printf("%lf\t%lf \n", current_solution.objectives[0],current_solution.objectives[1]);
    xreal[0]=current_solution.objectives[0];
    xreal[1]=current_solution.objectives[1];
    printf("%lf\t%lf \n", xreal[0],xreal[1]);*/

    //print_conformation_compact(&current_solution);
    a_ind(ind,&current_solution);
    //frees the memory of the current solution
    free_memory_conformation(&current_solution);

    /*for (int i = 0; i < encoding_len; ++i)
    {
        printf("%i\n", xreal->absolute_encoding[i]);
    }
    printf("\n");*/
}


void a_psp(individual *ind,struct conformation_def *to_psp){
    int i;
    
    // Encoding
    copy_vector_int(    ind->absolute_encoding, 
                    to_psp->absolute_encoding, 
                    encoding_len);
        
     // Coordinates
    for(i=0; i<sequence_len; i++){
        copy_vector_double( ind->coordinates[i], 
                        to_psp->coordinates[i],
                        dimensions);
    }
    
    // Objectives and HHtc
    to_psp->HHtc = ind->HHtc;  

    for(i=0; i<n_objectives; i++) to_psp->objectives[i] = ind->objectives[i];  
    
    // Contacts
        for(i=0; i<ind->HHtc; i++){        
        to_psp->contacts[i][0] = ind->contacts[i][0];
        to_psp->contacts[i][1] = ind->contacts[i][1];
    }

    //~ to->energy = from->energy;
    to_psp->valid = ind->valid;
   
   return;
}

void a_ind(individual *ind,struct conformation_def *current_solution){

    int i;

    for(i=0; i<encoding_len; i++) 
        {
            
            ind->absolute_encoding[i] = current_solution->absolute_encoding[i];
        }
        

     // Coordinates
    for(i=0; i<sequence_len; i++){
        copy_vector_double( current_solution->coordinates[i],
                            ind->coordinates[i], 
                        dimensions);
    }

    ind->HHtc = current_solution->HHtc;

    for(i=0; i<n_objectives; i++) {
        ind->objectives[i] = current_solution->objectives[i];
        //ind->xreal[i] = current_solution->objectives[i];
        ind->obj[i]=current_solution->objectives[i];
    }  

    for(i=0; i<current_solution->HHtc; i++){        
        ind->contacts[i][0] = current_solution->contacts[i][0];
        ind->contacts[i][1] = current_solution->contacts[i][1];
    }

    ind->valid = current_solution->valid;
    return;

}


void Muta(individual *parent1){

   
   struct conformation_def current_solution,hijo1;

    allocate_memory_conformation(&current_solution);
    allocate_memory_conformation(&hijo1);
    //individual a conofrmacion
    //generate_valid_random_conformation(&current_solution);
    //generate_valid_random_conformation(&hijo1);
    /*for (int i = 0; i < encoding_len; ++i)
    {
        printf("%i\n", parent1->absolute_encoding[i]);
    }*/
    a_psp(parent1,&current_solution);

    //print_conformation(&current_solution);
    //printf("current_solution\n");
    evaluate_conformation(&current_solution);   
    
    //evaluate_conformation(&mutate);
    //printf("current_solution\n");
    //print_conformation(&current_solution);


    //printf("antes\n");
    mutation(&current_solution,&hijo1);

    evaluate_conformation(&hijo1);
    //printf("hijo\n");
    //print_conformation(&hijo1);

    /*  printf("hijo1\n");
    print_conformation(&hijo1);
    printf("hiko2\n");
    print_conformation(&hijo2);*/

    a_ind(parent1,&hijo1);

    free_memory_conformation(&current_solution);
    free_memory_conformation(&hijo1);


}   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Cruc(individual *parent1, individual *parent2, individual *child1, individual *child2){

    
    struct conformation_def current_solution,current_solution2,hijo1,hijo2;

    allocate_memory_conformation(&current_solution);
    allocate_memory_conformation(&current_solution2);

    allocate_memory_conformation(&hijo1);
    allocate_memory_conformation(&hijo2);
    //individual a conofrmacion

    a_psp(parent1,&current_solution);

    a_psp(parent2,&current_solution2);  
        

    switch(type_Crossover)
    {
        case 0:
                //printf("cycle\n");
        return cruce_Cycle(&current_solution,&current_solution2,&hijo1,&hijo2);
        break;
        case 1:
            //printf("cruce\n");
        return cruce(&current_solution,&current_solution2,&hijo1,&hijo2);
        break;

        case 2:
            //printf("Static\n");
        return cruce_static(&current_solution,&current_solution2,&hijo1,&hijo2);
        break;

        case 3:
            //printf("Non_wrap\n");
        return NON_WRAPPING_ORDER_CROSSOVER(&current_solution,&current_solution2,&hijo1,&hijo2);
        break;

        case 4:
               //printf("One\n");
        return  cruce_one_point(&current_solution,&current_solution2,&hijo1,&hijo2);
        break;

        default:
           // printf("two\n");
        return cruce_two_point(&current_solution,&current_solution2,&hijo1,&hijo2);
        break;
    }

    evaluate_conformation(&hijo1);
    evaluate_conformation(&hijo2);


    /*printf("hijo1\n");
    print_conformation(&hijo1);
    printf("hijo2\n");
    print_conformation(&hijo2);*/

    a_ind(child1,&hijo1);
    a_ind(child2,&hijo2);

    //print_conformation_compact(&hijo1);

    //print_conformation_compact(&hijo2);

    free_memory_conformation(&current_solution);
    free_memory_conformation(&hijo1);
    free_memory_conformation(&current_solution2);
    free_memory_conformation(&hijo2);

    
}   

void print_pop(NSGA2Type *nsga2Params,population *pop,int size){

    for (int i = 0; i < size; ++i)
    {
        print_ind(nsga2Params,&pop->ind[i]);
    }

}


void print_ind(NSGA2Type *nsga2Params,individual *ind){
    int i;
    
    printf("\nRank: %i\n",ind->rank );
    //printf("constr_violation: %f\n",ind->constr_violation);

    /*printf("Xreal:\t");
    for (i = 0; i < nsga2Params->nreal; ++i) printf("%f\t", ind->xreal[i]);*/

        printf("\nObjetivos: \t");
    for (i = 0; i < nsga2Params->nobj; ++i)
        printf("%f\t", ind->obj[i]);

        printf("\ncrowding_distance: %.2f\n", ind->crowd_dist);

    // Absolute encoding
    printf("\nAbsolute encoding: ");
    /*if (letters)*/
        for(i=0; i<encoding_len; i++) printf("%c ", letter[ind->absolute_encoding[i]]);
    /*else*/
    printf("\nAbsolute encoding: ");
        for(i=0; i<encoding_len; i++) printf("%d ", ind->absolute_encoding[i]);

    

    
    // Coordinates
    /*printf("\n\nCoordinates: \n");
    for(i=0; i<sequence_len; i++){ 
        //printf("[\t");
        for(int j=0; j<dimensions; j++) printf("%.2lf\t", ind->coordinates[i][j]);
        //printf("]");
        printf("\n");
    }*/
    
    printf("\nObjectives: [ ");
    for(i=0; i<n_objectives; i++) printf("%lf ", ind->objectives[i]);  
    printf("]\n");
    printf("\nHHtc: %d\n", ind->HHtc);
    printf("Valid: %d\n", ind->valid); 

    /*printf("\nConformation\n");
    for(i=0; i<ind->HHtc; i++){        
        printf("%i\n", ind->contacts[i][0]); 
        printf("%i\n", ind->contacts[i][1]);  
    }*/
}
void Prints(individual *child1, individual *child2){
    struct conformation_def current_solution,current_solution2;

    allocate_memory_conformation(&current_solution);
    allocate_memory_conformation(&current_solution2);

    a_psp(child1,&current_solution);
    a_psp(child2,&current_solution2);

    print_conformation_compact(&current_solution);
    print_conformation_compact(&current_solution2);

    free_memory_conformation(&current_solution);
    free_memory_conformation(&current_solution2);
}

char * cross(int type){

    char *Cycle = "/Cycle";
    char *Cruce = "/Cruce";
    char *Static = "/Static";
    char *Non_wrap = "/Non_wrap";
    char *One = "/One";
    char *Two = "/Two";
    switch(type)
    {
        case 0:
        return Cycle;
        break;
        case 1:
        return Cruce;
        break;

        case 2:
        return Static;
        break;

        case 3:
        return Non_wrap;
        break;

        case 4:
        return One;
        break;

        default:
        return Two;
        break;
    }

}