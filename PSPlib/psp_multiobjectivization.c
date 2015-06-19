////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	NSGA-II for Protein Structure Prediction (PSP)  in the HP Lattice Model 
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	22/08/2011
		
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_DYNRES_OBJ(){
	int i;
	
	for(i=0; i<DYNRES_MAXHH; i++){
		DYNRES_OBJ[i] = rndpsp.RandInt(0, n_objectives-1);
	}
	
}	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void compute_subsets(){
	int i;
	// DET - Deterministic: split Hs into 'n_subsets' equally-sized subsets
	// RND - Random: randomly organize Hs into 'n_subsets' groups
	// DYNk - Dynamic: like RND, but the subsets are re-computed after k generations without move
	
	if(strcmp(strategy_variant, "DET") == 0){
		
		if(n_subsets == 2){
			
				int n = ( rndpsp.RandU01() <= 0.5 ) ? (int)ceil(NH/2.0) : (int)floor(NH/2.0);
			
				for(i=0; n>0 && i<sequence_len; i++) 
					if (sequence[i]=='H'){
						Hsubset[i] = 1;
						n--;
					}else Hsubset[i] = 0;
					
				for(; i<sequence_len; i++) 
					Hsubset[i] = (sequence[i]=='H') ? 2 : 0;
					
		}else if(n_subsets == 3){
			
				int n1 = ( rndpsp.RandU01() <= 0.5 ) ? (int)ceil(NH/3.0) : (int)floor(NH/3.0);
				int n2 = ( rndpsp.RandU01() <= 0.5 ) ? (int)ceil((NH-n1)/2.0) : (int)floor((NH-n1)/2.0);
			
				for(i=0; n1>0 && i<sequence_len; i++) 
					if (sequence[i]=='H'){
						Hsubset[i] = 1;
						n1--;
					}else Hsubset[i] = 0;				
				
				for(; n2>0 && i<sequence_len; i++) 
					if (sequence[i]=='H'){
						Hsubset[i] = 2;
						n2--;
					}else Hsubset[i] = 0;
					
				for(; i<sequence_len; i++) 
					Hsubset[i] = (sequence[i]=='H') ? 3 : 0;

			

		}	
		
	}else{
		
		for(i=0; i<sequence_len; i++) 
			Hsubset[i] = (sequence[i]=='H') ? rndpsp.RandInt(1, n_subsets) : 0;
		
	}
	
	//~ printf("\n");
	//~ for(i=0; i<sequence_len; i++) printf("%d ", Hsubset[i]);
	//~ printf("\n");	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Makes a copy of a conformation
void formulation_adjustment(){

	if(strcmp(problem_formulation, "SO") == 0){
		
		n_objectives=1;
		
	}else if(strcmp(problem_formulation, "MO") == 0){
		
		// Decomposition Strategies
		if(strcmp(evaluation_strategy, "DEC1") == 0){
			 
			n_subsets=2;
			n_objectives=2;
			 
		}else if(strcmp(evaluation_strategy, "DEC2") == 0){
			
			n_subsets=2;
			n_objectives=2;
			
		}else if(strcmp(evaluation_strategy, "DEC3") == 0){
			
			n_subsets=2;
			n_objectives=3;
			
		}else if(strcmp(evaluation_strategy, "DEC4") == 0){
			
			n_subsets=3;
			n_objectives=2;
			
		}else if(strcmp(evaluation_strategy, "DEC5") == 0){
			
			n_subsets=3;
			n_objectives=4;
			
		}else if(strcmp(evaluation_strategy, "PARITY") == 0){
			
			n_objectives=2;
			
		}else if(strcmp(evaluation_strategy, "LOCALITY") == 0){
			
			n_objectives=2;
			
			localitysize =7;
			
		}else if(strcmp(evaluation_strategy, "MOK99_LOCALITY") == 0){
			
			n_objectives=2;
			
			localitysize = (int)atoi(strategy_variant);
			
		}else if(strcmp(evaluation_strategy, "MOK99_PARITY") == 0){
			
			n_objectives=2;
					
		}else if(strcmp(evaluation_strategy, "RESIDUE") == 0){			
			
			if(dimensions == 2) n_objectives=3;
			else if(dimensions == 3) n_objectives=5;
			
		}else if(strcmp(evaluation_strategy, "RESIDUE2") == 0){			
			
			if(dimensions == 2) n_objectives=2;
			else if(dimensions == 3) n_objectives=3;
			
		}else if(strcmp(evaluation_strategy, "K99M") == 0){			
			
			n_objectives=2;
			
		}else if(strcmp(evaluation_strategy, "I09_M") == 0){			
			
			n_objectives=2;
			
		}else if(strcmp(evaluation_strategy, "K99M2") == 0){			
			
			n_objectives=2;
			
		}else if(strcmp(evaluation_strategy, "B08M") == 0){			
			
			n_objectives=2;
			n_subsets=3;
			
		}else if(strcmp(evaluation_strategy, "DYNRES2") == 0){
			
			n_objectives=2; // Always two objectives,  but some can be not considered (=0.0)
			
			if(dimensions == 2) DYNRES_MAXHH=3;
			else if(dimensions == 3) DYNRES_MAXHH=5; 
			
		}else if(strcmp(evaluation_strategy, "DYNRES3") == 0){
			
			n_objectives=3; // Always three objectives,  but some can be not considered (=0.0)
			
			if(dimensions == 2) DYNRES_MAXHH=3;
			else if(dimensions == 3) DYNRES_MAXHH=5; 
			
		}			
		
		// Subsets type for Decomposition strategies
		if(strncmp(evaluation_strategy, "DEC", 3) ==0){
			
			if(strcmp(strategy_variant, "DYN") == 0) max_iwm=0;
			else if(strcmp(strategy_variant, "DYN10") == 0) max_iwm=10;
			else if(strcmp(strategy_variant, "DYN20") == 0) max_iwm=20;
			else if(strcmp(strategy_variant, "DYN30") == 0) max_iwm=30;
			else if(strcmp(strategy_variant, "DYN25") == 0) max_iwm=25;
			else if(strcmp(strategy_variant, "DYN50") == 0) max_iwm=50;	
			else max_iwm=-1;
		}
		
		// Multiobjectivization of alternative evaluation functions	
		if(strcmp(evaluation_strategy, "I09a") == 0){
			 
			n_objectives=3;
			 
		}else if(strcmp(evaluation_strategy, "I09b") == 0){
			 
			n_objectives=2;
			 
		}else if(strcmp(evaluation_strategy, "I09c") == 0){
			 
			n_objectives=2;
			 
		}else if(strcmp(evaluation_strategy, "I09d") == 0){
			 
			n_objectives=2;
			 
		}
		
	}
	

	
}

