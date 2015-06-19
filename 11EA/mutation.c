////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	(1+1) EA for Protein Structure Prediction (PSP)  in the HP Lattice Model 
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	08/08/2011
	
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Uniform mutation: change each position according to a probability of p_mut
// If no change in this position leads to a valid conformation, then restore the original value
void mutation(struct conformation_def *current_solution, struct conformation_def *mutated_solution){
	int i, j, rnddir[ndirections], original, valid, changed = 0;	
	double p_mut = (1.0 / encoding_len);
	
	// Make a copy of the conformation
	copy_conformation(current_solution, mutated_solution);
	
	// Initial permutation of directions
	for(i=0; i<ndirections; i++) rnddir[i] = i;
	
	do{
	
		// Mutate each position with prob. p_mut
		for(i=0; i<encoding_len; i++){
			if( rndpsp.RandU01() <= p_mut ){
				
				// Backup original value
				original = mutated_solution->absolute_encoding[i];
				
				// Get a permutation of directions
				shuffle(rnddir, ndirections);
				
				// Randomly explore the possible mutations for the position
				// until a valid configuration is obtained
				for(j=0, valid=0; j<ndirections && valid==0; j++){
					if(rnddir[j] != original){
						mutated_solution->absolute_encoding[i] = rnddir[j];
						update_coordinates(mutated_solution, i, encoding_len-1);
						valid = isValid(mutated_solution); // Verify					
					}				
				}
				
				// If no mutations in this position leads to a valid configuration
				// return to the original value
				if(!valid) {
					mutated_solution->absolute_encoding[i] = original;
					update_coordinates(mutated_solution, i, encoding_len-1);
				}else{
					changed = 1;
				}
			}
		}
		
	}while(changed == 0);
}

