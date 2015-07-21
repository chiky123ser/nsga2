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
	
	Execution:
	
		./   ... 
		
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "../PSPlib/psp.h"

struct conformation_def current_solution,mutate,hijo1,hijo2;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Implements the main cycle of (1+1) EA */
void do_11ea(){
	
	// Allocating memory
	allocate_memory_conformation(&current_solution);	
	allocate_memory_conformation(&mutate);
	allocate_memory_conformation(&hijo1);	
	allocate_memory_conformation(&hijo2);	
	
	 // Generate initial solution at random		
	/*for (int i = 0; i < 40; ++i)
	{*/
	generate_valid_random_conformation(&current_solution);	
	
	
	evaluate_conformation(&current_solution);	
	print_conformation(&current_solution);

	//cruce_Cycle();
	mutation(&current_solution, &mutate);
	print_conformation(&mutate);
	//printf("%i\n", i);


	/*
	generate_valid_random_conformation(&mutate);
	evaluate_conformation(&mutate);
	
	print_conformation(&current_solution);

	print_conformation(&mutate);
	//cruce(&current_solution,&mutate,&hijo1,&hijo2);
	//cruce_static(&current_solution,&mutate,&hijo1,&hijo2);
	//cruce_Cycle(&current_solution,&mutate,&hijo1,&hijo2);
	//NON_WRAPPING_ORDER_CROSSOVER(&current_solution,&mutate,&hijo1,&hijo2);
	//cruce_one_point(&current_solution,&mutate,&hijo1,&hijo2);
	//cruce_two_point(&current_solution,&mutate,&hijo1,&hijo2);

	evaluate_conformation(&hijo1);
	
	print_conformation(&hijo1);

	evaluate_conformation(&hijo2);
	
	print_conformation(&hijo2);*/
	
}	











