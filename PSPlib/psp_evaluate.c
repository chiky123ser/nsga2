////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Protein Structure Prediction (PSP) in the HP Lattice Model
	
	(psp_evaluate.c): Energy (evaluation) functions and related ones
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	08/08/2011
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Invokes the evaluation of the given conformation
void evaluate_conformation(struct conformation_def *conformation){
	
	// Compute the number of HH topological contacts (the original objective)
	getHHTopologicalContacts(conformation);
	
	// Check validity of the conformation
	conformation->valid = isValid(conformation);	
	
	// Evaluate the conformation according to the given alternative strategy
	if(n_objectives == 1)	 
		evaluate_conformation_singleobjective(conformation);
				
	else 
		evaluate_conformation_multiobjective(conformation);			
	
	
	// Increase evaluations counter
	/*evaluation++;
	conformation->evaluation = evaluation;
	
	// Save the best evaluated solution
	if( conformation->valid == 1  &&  conformation->HHtc > evaluated_best.HHtc ){
		copy_conformation(conformation, &evaluated_best);
	}*/
	
	// If evaluation number is a multiple of "report_each", then report progress
	//~ if( report_each!=0 && evaluation%100 == 0 ){			
		//~ printf("%ld\t%d\t%d\t%ld\n", evaluation, evaluated_best.HHtc, evaluated_best.valid, evaluated_best.evaluation);		
	//~ }
	
	//~ int current_energy =  evaluated_best.HHtc*-1;
	//~ if( current_energy != online_energy ){
		
		//~ if( evaluation > 1 )
			//~ printf("%ld\t%d\n", evaluation-1, online_energy);		
		//~ online_energy = current_energy;
		//~ printf("%ld\t%d\n", evaluation, online_energy);		
	//~ }
	
	
}
