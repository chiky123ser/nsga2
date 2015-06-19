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
// Pareto dominance. 
// Returns 1 if vec1 dominates vec2, 2 if vec2 dominates vec1, 0 otherwise (incomparable)
int pareto_dominance(double *vec1, double *vec2, int size){
	int i, vec1cont=0, vec2cont=0;
	
	for(i=0; i<size; i++){			
		if(vec1[i]  < vec2[i]) vec1cont++;
		else if(vec1[i]  > vec2[i]) vec2cont++;		
	}
	
	if(vec2cont==0 && vec1cont>0) return 1;
	if(vec1cont==0 && vec2cont>0) return 2;
	
	return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns 1 if c1 is noninferior (not worse, nondominated) with regard c2
int noninferior(struct conformation_def *c1, struct conformation_def *c2){
	
	// Multiobjective (general) case
	return ( pareto_dominance(c1->objectives, c2->objectives, n_objectives) != 2 ) ? 1 : 0;
}

