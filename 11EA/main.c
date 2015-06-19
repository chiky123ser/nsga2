////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	(1+1) EA for Protein Structure Prediction (PSP)  in the HP Lattice Model 
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	08/08/2011
	
	----------------------------------------------------------------------------------------------------------------
	
	Compilation: Use the "Makefile"
	
	----------------------------------------------------------------------------------------------------------------
		
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "11ea.c"
#include <time.h>
#include <stdlib.h>


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main function / User interface
int main(int argc, char *argv[]){	
	
	// Reading input parameters
	if(argc!=8) 
	printf("\n\tUSAGE:  ./11ea   LATTICE_FILE   PROBLEM_INSTANCE   FORMULATION={SO, MO}   EVALUATION={D85, K99..., DEC1, DEC2...}   VARIANT = {DET, RND, DYN, DYN10..30}   MAX_EVAL   RUN\n\n"), exit(1);	
	
	// Load lattice configuration
	load_lattice_configuration(argv[1]);
		
	// Read problem instance
	read_problem_instance(argv[2]);
	
	// Problem formulation to be used
	strcpy(problem_formulation, argv[3]);
	
	// Evaluation function or multiobjectivization strategy
	strcpy(evaluation_strategy, argv[4]);	
	
	// Type of subsets to use: DET, RND, DYN
	strcpy(strategy_variant, argv[5]);

	// Allowed number of energy function evaluations
	max_evaluations = atoi(argv[6]);	
	
	// Run and seed for random numbers
	srand(time(NULL)); 
	int variaable=rand() %5 ;
	run = atoi(argv[7]);
	seed = run*sequence_len*dimensions+variaable;
	
	initialize_rnd(seed);	
		
	// Parameters Adjustment
	formulation_adjustment();
	
	// Main procedure
	do_11ea();
	
	return 0;
}










