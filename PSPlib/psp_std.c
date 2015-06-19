/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Protein Structure Prediction (PSP) in the HP Lattice Model
	
	(psp_functions.c): PSP standard functions (load settings, printing, etc.)
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	15/02/2011
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int continue_searching(){
	
	if(	(max_evaluations == 0 || evaluation < max_evaluations)  // If the maximum number of evaluations has not been reached (or no maximum defined)
	&& 	(evaluated_best.HHtc < optimalHHtc || evaluated_best.valid==0) // and the optimum (feasible) solution has not been found
	) 
		return 1;
	
	// Otherwise
	return 0;	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initialize_rnd(unsigned long seed){
	unsigned long seeds[6] = {seed, seed, seed, seed, seed, seed};
	rndpsp.SetSeed(seeds);
	rndpsp.IncreasedPrecis(true);
	rndpsp.ResetStartStream();	
	rndpsp.RandInt(0,3); rndpsp.RandU01();		
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Reads the sequence to be solved
void read_problem_instance(char *filename){
	int i;
	char c;
	
	FILE *fp=fopen(filename,"r");
	
	if(fp){
		// Length of the sequence
		fscanf(fp, "%d", &sequence_len); 
	
		// Absolute encoding is (sequence_len - 1) long
		// Since the position of the first amino acid is fixed, the encoding only
		// includes the reimaining n-1 decisions 
		encoding_len = sequence_len - 1;
		
		// Relative encoding is (sequence_len - 1) long
		// Since the first move is assumed to be 'forward', the encoding only
		// includes the reimaining n-1 decisions 
		//~ relative_encoding_len = sequence_len - 2;
	
		// Reading sequence
		sequence = (char *)malloc(sizeof(char) * sequence_len);
		for(i=0, NH=0, NP=0; i<sequence_len; ) {
			fscanf(fp, "%c", &c);
			if(c=='H' || c=='P') {
				sequence[i] = c;
				if(c=='H') NH++;
				else NP++;
				i++;
			}
		}
	
		// Best known energy (for reference)
		fscanf(fp, "%d", &optimalHHtc); 
	}else
		printf("\n\n\tERROR: Problem instance file not foud\n\n"), exit(1);
	
	fclose(fp);	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Loads lattice configuration from the input file
void load_lattice_configuration(char *filename){
	int i,j;
	char c;
	
	FILE *fp=fopen(filename,"r");
	
	if(fp){		
		// Dimensions of the lattice
		fscanf(fp, "%d", &dimensions); 
	
		// Number of direction vectors
		fscanf(fp, "%d", &ndirections); 
	
		// Reading  direction vectors
		directions = (double **)malloc(sizeof(double *) * ndirections);
		for(i=0; i<ndirections; i++) {
			directions[i] = (double *)malloc(sizeof(double) * dimensions);		
			for(j=0; j<dimensions; j++) fscanf(fp, "%lf", &directions[i][j]);		
		}
	
		// Letter equivalences
		fscanf(fp, "%d", &letters); 
		if(letters){
			letter=(char *)malloc(sizeof(char)*ndirections);
			for(i=0; i<ndirections; i++){
				fscanf(fp, "%c", &c);
				if(c=='\n') i--;
				else letter[i]=c;
			}
		}
	}else 
		printf("\n\n\tERROR: Lattice configuration file not foud\n\n"), exit(1);		
	
	fclose(fp);		
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Allocates memory for the given conformation and initializes its properties
void allocate_memory_conformation(struct conformation_def *conformation){
	int i;
	
	conformation->absolute_encoding = (int *)malloc(sizeof(int) * encoding_len);
	conformation->HHtc=0;		
	//~ conformation->energy=0;	
	conformation->objectives = (double *)malloc(sizeof(double) * n_objectives);	
	for(i=0; i<n_objectives; i++) conformation->objectives[i]=0;	
	
	// We can store as many contacts as optimalHHtc
	conformation->contacts = (int **)malloc(sizeof(int *) * 2*optimalHHtc);
	for(i=0; i<2*optimalHHtc; i++) conformation->contacts[i] = (int *)malloc(sizeof(int) * 2);

	conformation->contactsP = (int **)malloc(sizeof(int *) * 2*optimalHHtc);
	for(i=0; i<2*optimalHHtc; i++) conformation->contactsP[i] = (int *)malloc(sizeof(int) * 2);
	
	conformation->rankpsp=0;		
	conformation->cwd=0.0;		
	conformation->valid=0;		
	conformation->evaluation=0;	
		
	// Coordinates for each amino acid in the sequence
	conformation->coordinates = (double **)malloc(sizeof(double *) * sequence_len);
	for(i=0; i<sequence_len; i++) 
		conformation->coordinates[i] = (double *)malloc(sizeof(double) * dimensions);
	
	// The first residue is fixed at the origin
	for(i=0; i<dimensions; i++) conformation->coordinates[0][i] = 0;		
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Liberates the memory of the given conformation
void free_memory_conformation(struct conformation_def *conformation){
	int i;
	
	free(conformation->absolute_encoding);
	free(conformation->objectives);
	for(i=0; i<sequence_len; i++)  free( conformation->coordinates[i] );
	free( conformation->coordinates );
	for(i=0; i<2*optimalHHtc; i++) free( conformation->contacts[i] );
	free( conformation->contacts );

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Prints all the information of the given solution
void print_conformation(struct conformation_def *conformation){
	int i;
	
	// Absolute encoding
	printf("\nAbsolute encoding: ");
	/*if (letters)*/
		for(i=0; i<encoding_len; i++) printf("%c ", letter[conformation->absolute_encoding[i]]);
	/*else*/
		printf("\nAbsolute encoding: ");
		for(i=0; i<encoding_len; i++) printf("%d ", conformation->absolute_encoding[i]);

	/*printf("\nConformation\n");
	for(i=0; i<conformation->HHtc; i++){		
		printf("%i\n", conformation->contacts[i][0]); 
		printf("%i\n", conformation->contacts[i][1]);  
	}*/
	
	// Coordinates
	printf("\n\nCoordinates: \n");
	for(i=0; i<sequence_len; i++){ 
		//printf("[\t");
		for(int j=0; j<dimensions; j++) printf("%.2lf\t", conformation->coordinates[i][j]);
		//printf("]");
		printf("\n");
	}
	
	printf("\nObjectives: [ ");
	for(i=0; i<n_objectives; i++) printf("%lf ", conformation->objectives[i]);	
	printf("]\n");
	printf("\nHHtc: %d\n", conformation->HHtc);
	printf("Valid: %d\n", conformation->valid);	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Prints the information of the given solution using one line
void print_conformation_compact(struct conformation_def *conformation){
	int i;
	
	for(i=0; i<sequence_len; i++){ 
		//printf("[\t");
		for(int j=0; j<dimensions; j++) printf("%.2lf\t", conformation->coordinates[i][j]);
		//printf("]");
		printf("\n");
	}
	
	// Encoding
	/*if (letters)
		for(i=0; i<encoding_len; i++) printf("%c ", letter[conformation->absolute_encoding[i]]);
	else
		for(i=0; i<encoding_len; i++) printf("%d ", conformation->absolute_encoding[i]);*/
		
	//~ printf("\tE: %lf", conformation->energy);
	
	/*for(i=0; i<n_objectives; i++) {
		
		printf("%lf ", conformation->objectives[i]);	
		printf("\t");}
		printf("\n");
	//printf("\tHH: %d", conformation->HHtc);	
	printf("\t%s\n", (conformation->valid)? "VALID" : "INVALID");*/	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Makes a copy of a conformation
 void copy_conformation(struct conformation_def *from, struct conformation_def *to){
	int i;
	
	// Encoding
	copy_vector_int(	from->absolute_encoding, 
					to->absolute_encoding, 
					encoding_len);
		
	 // Coordinates
	for(i=0; i<sequence_len; i++){
		copy_vector_double(	from->coordinates[i], 
						to->coordinates[i],
						dimensions);
	}
	
	// Objectives and HHtc
	to->HHtc = from->HHtc;	
	to->PPtc = from->PPtc;
	for(i=0; i<n_objectives; i++) to->objectives[i] = from->objectives[i];	
	
	// Contacts
	for(i=0; i<from->HHtc; i++){		
		to->contacts[i][0] = from->contacts[i][0];
		to->contacts[i][1] = from->contacts[i][1];
		to->contactsP[i][0] = from->contactsP[i][0];
		to->contactsP[i][1] = from->contactsP[i][1];
	}
	
	//~ to->energy = from->energy;
	to->valid = from->valid;
	to->rankpsp = from->rankpsp;	
	to->cwd = from->cwd;	
	to->evaluation = from->evaluation;	
}





