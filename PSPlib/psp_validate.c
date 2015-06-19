////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Protein Structure Prediction (PSP) in the HP Lattice Model
	
	(psp_generate_validate.c): Validation of conformations
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	15/02/2011
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Verifies if the given solution is valid (collision-free)
// Returns 1 if valid, 0 otherwise
 int isValid(struct conformation_def *conformation){
	int i, j;
	
	// A solution is valid if it is free of collisions; that is, if there is no location
	// occupied by more than one amino acid. This verification implicitly
	// considers reverse moves.
	
	// Compare each pair of coordinates
	for(i=0; i<sequence_len; i++){
		for(j=i+1; j<sequence_len; j++){
			if( isEqual_double(conformation->coordinates[i], conformation->coordinates[j], dimensions) ) 
				return 0;
		}
	}
	
	return 1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Verifies if the given partial solution is valid (collision-free)
// Returns 1 if valid, 0 otherwise
int isValid_partial(struct conformation_def *conformation, int ncoord){
	int i, j;
	
	// A solution is valid if it is free of collisions; that is, if there is no location
	// occupied by more than one amino acid. This verification implicitly
	// considers reverse moves.
	
	// Compare each pair of coordinates
	for(i=0; i<ncoord; i++){
		for(j=i+1; j<ncoord; j++){
			if( isEqual_double(conformation->coordinates[i], conformation->coordinates[j], dimensions) ) 
				return 0;
		}
	}
	
	return 1;
}
 
