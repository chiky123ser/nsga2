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
/* *Dynamic Residue-level decomposition */
void DYNRES(struct conformation_def *conformation){		
	int i;	
	
	// Initialize objective values	
	for(i=0; i<n_objectives; i++)  conformation->objectives[i] = 0.0;
	
	// En cuantos contactos esta involucrado cada aminoacido H
	int cuantos[sequence_len];	
	for(i=0; i<sequence_len; i++) cuantos[i]=0;	
	for(i=0; i<conformation->HHtc; i++){			 
		cuantos[conformation->contacts[i][0]]++;
		cuantos[conformation->contacts[i][1]]++;			 
	}
	
	int obi[6]={0,0,0,0,0,0};			
	for(i=0; i<sequence_len; i++) obi[cuantos[i]]++;	
	
	if(dimensions == 2){	
		
		conformation->objectives[DYNRES_OBJ[0]] -= 1.0*(obi[1]/2.0);
		conformation->objectives[DYNRES_OBJ[1]] -= 1.0*obi[2];
		conformation->objectives[DYNRES_OBJ[2]] -= 3.0*obi[3]/2.0;
		
	}else{
	
		conformation->objectives[DYNRES_OBJ[0]] -= 1.0*(obi[1]/2.0);
		conformation->objectives[DYNRES_OBJ[1]] -= 1.0*obi[2];
		conformation->objectives[DYNRES_OBJ[2]] -= 3.0*(obi[3]/2.0);
		conformation->objectives[DYNRES_OBJ[3]] -= 4.0*(obi[4]/2.0);		
		conformation->objectives[DYNRES_OBJ[4]] -= 5.0*(obi[5]/2.0);		
		
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* *Residue-level decomposition */
void RESIDUE2(struct conformation_def *conformation){		
	int i;	
	
	// En cuantos contactos esta involucrado cada aminoacido H
	int cuantos[sequence_len];	
	for(i=0; i<sequence_len; i++) cuantos[i]=0;	
	for(i=0; i<conformation->HHtc; i++){			 
		cuantos[conformation->contacts[i][0]]++;
		cuantos[conformation->contacts[i][1]]++;			 
	}
	
	int obi[6]={0,0,0,0,0,0};			
	for(i=0; i<sequence_len; i++) obi[cuantos[i]]++;	
	
	if(dimensions == 2){
	
		conformation->objectives[0] = -1.0*(obi[1]/2.0) - 3.0*obi[3]/2.0;
		conformation->objectives[1] = -1.0*obi[2];
		
	}else{
	
		conformation->objectives[0] = -1.0*(obi[1]/2.0)  - 5.0*(obi[5]/2.0);
		conformation->objectives[1] = -1.0*obi[2] - 4.0*(obi[4]/2.0);
		conformation->objectives[2] = -3.0*(obi[3]/2.0);
		
	}

}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* *Residue-level decomposition */
void RESIDUE(struct conformation_def *conformation){		
	int i;	
	
	// En cuantos contactos esta involucrado cada aminoacido H
	int cuantos[sequence_len];	
	for(i=0; i<sequence_len; i++) cuantos[i]=0;	
	for(i=0; i<conformation->HHtc; i++){			 
		cuantos[conformation->contacts[i][0]]++;
		cuantos[conformation->contacts[i][1]]++;			 
	}
	
	int obi[6]={0,0,0,0,0,0};			
	for(i=0; i<sequence_len; i++) obi[cuantos[i]]++;	
	
	if(dimensions == 2){
	
		// Compute objective values
		conformation->objectives[0] = -1.0*(obi[1]/2.0);
		conformation->objectives[1] = -1.0*obi[2];
		conformation->objectives[2] = -3.0*obi[3]/2.0;
	}else{
	
		conformation->objectives[0] = -1.0*(obi[1]/2.0);
		conformation->objectives[1] = -1.0*obi[2];
		conformation->objectives[2] = -3.0*(obi[3]/2.0);
		conformation->objectives[3] = -4.0*(obi[4]/2.0);		
		conformation->objectives[4] = -5.0*(obi[5]/2.0);		
		
	}
		
	// Compute objective values
	//~ conformation->objectives[0] = -1.0*(obi[1]/2.0) - 3.0*obi[3]/2.0;
	//~ conformation->objectives[1] = -1.0*obi[2];	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* LOCALITY Decomposition applied to function K99 */
void MOK99_LOCALITY(struct conformation_def *conformation){		
	int i, j, k;
	double d, energy_contribution;
	
	// Initialize objective values
	conformation->objectives[0] = 0.0;
	conformation->objectives[1] = 0.0;
	
	// k=4 for the square lattice and k=5 for the cubic lattice
	k = (dimensions == 2) ? 4 : 5;	
	
	// Compute energy/objective values
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
			for(j=i+2; j<sequence_len; j++){
				if(sequence[j]=='H'){
					
					// Compute energy/objective values contriburion of this (H,H) pair
					d = euclidean_distance(conformation->coordinates[i], conformation->coordinates[j], dimensions);					
					energy_contribution = (ceil(d)>1.0) ? (1.0/(pow(d, k)*NH)) : 1;
					
					// Multiobjectivization
					if( (j - i) <= localitysize ) 
						conformation->objectives[0] -= energy_contribution;
					else 
						conformation->objectives[1] -= energy_contribution;
					
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* PARITY Decomposition applied to function K99 */
void MOK99_PARITY(struct conformation_def *conformation){		
	int i, j, k;
	double d, energy_contribution;
	
	// Initialize objective values
	conformation->objectives[0] = 0.0;
	conformation->objectives[1] = 0.0;
	
	// k=4 for the square lattice and k=5 for the cubic lattice
	k = (dimensions == 2) ? 4 : 5;	
	
	// Compute energy/objective values
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
			for(j=i+2; j<sequence_len; j++){
				if(sequence[j]=='H'){
					
					// Compute energy/objective values contriburion of this (H,H) pair
					d = euclidean_distance(conformation->coordinates[i], conformation->coordinates[j], dimensions);					
					energy_contribution = (ceil(d)>1.0) ? (1.0/(pow(d, k)*NH)) : 1;
					
					// Multiobjectivization
					if( (i%2) == 0 )
						conformation->objectives[0] -= energy_contribution;
					else 
						conformation->objectives[1] -= energy_contribution;
					
				}
			}
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* LOCALITY Decomposition */
void LOCALITY(struct conformation_def *conformation){		
	int i, r1, r2;
	
	// Initialize objective values
	conformation->objectives[0] = 0.0;
	conformation->objectives[1] = 0.0;
	
	for(i=0; i<conformation->HHtc; i++){
		
		r1 = conformation->contacts[i][0];
		r2 = conformation->contacts[i][1];
		
		if( (r2 - r1) <= localitysize ) 
			conformation->objectives[0]--;
		else 
			conformation->objectives[1]--;
	}	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Decomposition 1: 2 Hsubsets, Strategy 1, 2 Objectives
f1 = | HH contacts the first residue index is even |
f2 = | all other HH contacts |
*/
void PARITY(struct conformation_def *conformation){		
	int i, r1;
	
	// Initialize objective values
	conformation->objectives[0] = 0.0;
	conformation->objectives[1] = 0.0;
	
	for(i=0; i<conformation->HHtc; i++){
		
		r1 = conformation->contacts[i][0];
		//printf("%i\n",conformation->contacts[i][0] );
				
		if( (r1%2) == 0 ) 
			conformation->objectives[0]--;
		else 
			conformation->objectives[1]--;
	}		
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Decomposition 1: 2 Hsubsets, Strategy 1, 2 Objectives
f1 = | HH contacts where both Hs are in the first subset |
f2 = | all other HH contacts |
*/
void DEC1(struct conformation_def *conformation){		
	int i, r1, r2;
	
	// Initialize objective values
	conformation->objectives[0] = 0.0;
	conformation->objectives[1] = 0.0;
	
	for(i=0; i<conformation->HHtc; i++){
		
		r1 = conformation->contacts[i][0];
		r2 = conformation->contacts[i][1];
		
		if(Hsubset[r1]==1 and Hsubset[r2]==1) 
			conformation->objectives[0]--;
		else 
			conformation->objectives[1]--;
	}		
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Decomposition 2: 2 Hsubsets, Strategy 2, 2 Objectives
Increase f1 if both H are in the same subset, increase f2 otherwise
f1 = | HH contacts where both Hs are in the same subset (local interaction) |
f2 = | all other HH contacts (nonlocal interaction) |
*/
void DEC2(struct conformation_def *conformation){		
	int i, r1, r2;
	
	// Initialize objective values
	conformation->objectives[0] = 0.0;
	conformation->objectives[1] = 0.0;
	
	for(i=0; i<conformation->HHtc; i++){
		
		r1 = conformation->contacts[i][0];
		r2 = conformation->contacts[i][1];
		
		if( Hsubset[r1] == Hsubset[r2] ) 
			conformation->objectives[0]--;
		else 
			conformation->objectives[1]--;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Decomposition 3: 2 Hsubsets, Strategy 3, 3 Objectives
f1 = | HH contacts where both Hs are in the first subset (local interaction) |
f2 = | HH contacts where both Hs are in the second subset  (local interaction) |
f3 = | all other HH contacts (nonlocal interaction) |
*/
void DEC3(struct conformation_def *conformation){
	int i, r1, r2;
	
	// Initialize objective values
	conformation->objectives[0] = 0.0;
	conformation->objectives[1] = 0.0;
	conformation->objectives[2] = 0.0;
	
	for(i=0; i<conformation->HHtc; i++){
		
		r1 = conformation->contacts[i][0];
		r2 = conformation->contacts[i][1];
		
		if(Hsubset[r1]==1 and Hsubset[r2]==1) 			
			conformation->objectives[0]--;		
		else if(Hsubset[r1]==2 and Hsubset[r2]==2) 			
			conformation->objectives[1]--;		
		else 			
			conformation->objectives[2]--;
	}	
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Decomposition 5: 3 Hsubsets, Strategy 2, 4 Objectives
f1 = | HH contacts where both Hs are in the first subset (local interaction) |
f2 = | HH contacts where both Hs are in the second subset  (local interaction) |
f3 = | HH contacts where both Hs are in the third subset  (local interaction) |
f4 = | all other HH contacts (nonlocal interaction) |
*/
void DEC5(struct conformation_def *conformation){
	int i, r1, r2;
	
	// Initialize objective values
	conformation->objectives[0] = 0.0;
	conformation->objectives[1] = 0.0;
	conformation->objectives[2] = 0.0;
	conformation->objectives[3] = 0.0;
	
	for(i=0; i<conformation->HHtc; i++){
		
		r1 = conformation->contacts[i][0];
		r2 = conformation->contacts[i][1];
		
		if(Hsubset[r1]==1 and Hsubset[r2]==1) 			
			conformation->objectives[0]--;		
		else if(Hsubset[r1]==2 and Hsubset[r2]==2) 			
			conformation->objectives[1]--;	
		else if(Hsubset[r1]==3 and Hsubset[r2]==3) 			
			conformation->objectives[2]--;		
		else 			
			conformation->objectives[3]--;
	}			
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
K. Islam and M. Chetty. Novel Memetic Algorithm for Protein Structure Prediction. In Ann
Nicholson and Xiaodong Li, editors, AI 2009: Advances in Artificial Intelligence, volume 5866
of Lecture Notes in Computer Science, pages 412-421. Springer Berlin / Heidelberg, 2009.
*/
void Islam2009_multiobjective(struct conformation_def *conformation, int type){
	int i, j;
	double HP=0.0, H_min_max[dimensions][2], P_min_max[dimensions][2];
	double Hc[dimensions], Hcompliance=0.0, Pcompliance=0.0, list[dimensions*2];
	
	// Initialice min_max	
	for(j=0; j<dimensions; j++){
		H_min_max[j][0] = INF; H_min_max[j][1] = -INF;
		P_min_max[j][0] = INF; P_min_max[j][1] = -INF;
	}	
		
	// Compute hypotetical rectangles enclosing H and P residues
	for(i=0; i<sequence_len; i++){
		if(sequence[i]=='H'){				
			for(j=0; j<dimensions; j++){ 				
				if(conformation->coordinates[i][j] < H_min_max[j][0]) 
					H_min_max[j][0] = conformation->coordinates[i][j];				
				if(conformation->coordinates[i][j] > H_min_max[j][1]) 
					H_min_max[j][1] = conformation->coordinates[i][j];				
			}			
		}else{			
			for(j=0; j<dimensions; j++){ 				
				if(conformation->coordinates[i][j] < P_min_max[j][0]) 
					P_min_max[j][0] = conformation->coordinates[i][j];				
				if(conformation->coordinates[i][j] > P_min_max[j][1]) 
					P_min_max[j][1] = conformation->coordinates[i][j];				
			}
		}
	}
	
	// Compute H-compliance
	for(j=0; j<dimensions; j++){
		//~ Hc[j] = ( (H_min_max[j][1] - H_min_max[j][0]) / 2.0); // The original definition is wrong
		Hc[j] = ( (H_min_max[j][1] + H_min_max[j][0]) / 2.0); // This is the rigth expression
	}
	
	for(i=0; i<sequence_len; i++){
		if(sequence[i]=='H'){				
			for(j=0; j<dimensions; j++){ 				
				Hcompliance += pow(conformation->coordinates[i][j] - Hc[j], 2);
			}			
		}else{			
			for(j=0; j<dimensions; j++){ 
				list[2*j] = fabs( P_min_max[j][0] - conformation->coordinates[i][j] ); 
				list[2*j+1] = fabs( P_min_max[j][1] - conformation->coordinates[i][j] );
			}
			Pcompliance += min_double(list, dimensions*2);
		}
	}
	
	Hcompliance /= NH;
	Pcompliance /= NP;
	
	// Compute HP model original energy	
	HP = (-1.0*conformation->HHtc);
		
	// Final energy
	if(type==0){
		conformation->objectives[0] = HP;
		conformation->objectives[1] = Hcompliance;
		conformation->objectives[2] = Pcompliance;
	}else if(type==1){
		conformation->objectives[0] = HP;
		conformation->objectives[1] = Hcompliance;
	}else if(type==2){
		conformation->objectives[0] = HP;
		conformation->objectives[1] = Pcompliance;
	}else if(type==3){
		conformation->objectives[0] = Hcompliance;
		conformation->objectives[1] = Pcompliance;
	}
}

void Islam2009_energy_1(struct conformation_def *conformation){
	int i, j;
	double HP=0.0, H_min_max[dimensions][2], P_min_max[dimensions][2];
	double Hc[dimensions], Hcompliance=0.0, Pcompliance=0.0, list[dimensions*2];
	double alpha=10000.0;
	
	// Initialice min_max	
	for(j=0; j<dimensions; j++){
		H_min_max[j][0] = INF; H_min_max[j][1] = -INF;
		P_min_max[j][0] = INF; P_min_max[j][1] = -INF;
	}	
		
	// Compute hypotetical rectangles enclosing H and P residues
	for(i=0; i<sequence_len; i++){
		if(sequence[i]=='H'){				
			for(j=0; j<dimensions; j++){ 				
				if(conformation->coordinates[i][j] < H_min_max[j][0]) 
					H_min_max[j][0] = conformation->coordinates[i][j];				
				if(conformation->coordinates[i][j] > H_min_max[j][1]) 
					H_min_max[j][1] = conformation->coordinates[i][j];				
			}			
		}else{			
			for(j=0; j<dimensions; j++){ 				
				if(conformation->coordinates[i][j] < P_min_max[j][0]) 
					P_min_max[j][0] = conformation->coordinates[i][j];				
				if(conformation->coordinates[i][j] > P_min_max[j][1]) 
					P_min_max[j][1] = conformation->coordinates[i][j];				
			}
		}
	}
	
	// Compute H-compliance
	for(j=0; j<dimensions; j++){
		//~ Hc[j] = ( (H_min_max[j][1] - H_min_max[j][0]) / 2.0); // The original definition is wrong
		Hc[j] = ( (H_min_max[j][1] + H_min_max[j][0]) / 2.0); // This is the rigth expression
	}
	
	for(i=0; i<sequence_len; i++){
		if(sequence[i]=='H'){				
			for(j=0; j<dimensions; j++){ 				
				Hcompliance += pow(conformation->coordinates[i][j] - Hc[j], 2);
			}			
		}else{			
			for(j=0; j<dimensions; j++){ 
				list[2*j] = fabs( P_min_max[j][0] - conformation->coordinates[i][j] ); 
				list[2*j+1] = fabs( P_min_max[j][1] - conformation->coordinates[i][j] );
			}
			Pcompliance += min_double(list, dimensions*2);
		}
	}
	
	Hcompliance /= NH;
	Pcompliance /= NP;
	
	// Compute HP model original energy	
	HP = (-1.0*conformation->HHtc);
		
	// Final energy
	conformation->objectives[0] = Hcompliance + Pcompliance;
	conformation->objectives[1] = alpha*HP ;
}


void Krasnogor1999_energy_Multy(struct conformation_def *conformation){
	int i, j, k;
	double energy=0.0,energy1=0.0, d;
	// k=4 for the square lattice and k=5 for the cubic lattice
	k = (dimensions == 2) ? 4 : 5;	
	
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
			for(j=i+2; j<sequence_len; j++){
				if(sequence[j]=='H'){
					//~ d = manhattan_distance(conformation->coordinates[i], conformation->coordinates[j], dimensions);
					d = euclidean_distance(conformation->coordinates[i], conformation->coordinates[j], dimensions);
					//printf("%f\t%f\t%f\t%f\t%f\n",d ,conformation->coordinates[i][0], conformation->coordinates[i][1],conformation->coordinates[j][0], conformation->coordinates[j][1]);
					energy -= (ceil(d)>1.0) ? (1.0/(pow(d, k)*NH)) : 1;
				}
			}
		}
	}
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='P'){			
			for(j=i+2; j<sequence_len; j++){
				if(sequence[j]=='P'){
					//~ d = manhattan_distance(conformation->coordinates[i], conformation->coordinates[j], dimensions);
					d = euclidean_distance(conformation->coordinates[i], conformation->coordinates[j], dimensions);
					//printf("%f\t%f\t%f\t%f\t%f\n",d ,conformation->coordinates[i][0], conformation->coordinates[i][1],conformation->coordinates[j][0], conformation->coordinates[j][1]);
					energy1 -= (ceil(d)>1.0) ? (1.0/(pow(d, k)*NH)) : 1;
				}
			}
		}
	}
	
	conformation->objectives[0] = energy;
	//printf("%f\n",energy1 );
	conformation->objectives[1] = energy1;
	printf("%f\t%f\n",conformation->objectives[0],conformation->objectives[1]);

}
void Krasnogor1999_energy_Multy2(struct conformation_def *conformation){

	int i, j, k,px,py;
	double x,y;
	double energy=0.0,energy1=0.0, d;
	// k=4 for the square lattice and k=5 for the cubic lattice
	k = (dimensions == 2) ? 4 : 5;	
	
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
			for(j=i+2; j<sequence_len; j++){
				if(sequence[j]=='H'){
					//~ d = manhattan_distance(conformation->coordinates[i], conformation->coordinates[j], dimensions);
					d = euclidean_distance(conformation->coordinates[i], conformation->coordinates[j], dimensions);
					//printf("%f\t%f\t%f\t%f\t%f\n",d ,conformation->coordinates[i][0], conformation->coordinates[i][1],conformation->coordinates[j][0], conformation->coordinates[j][1]);
					energy -= (ceil(d)>1.0) ? (1.0/(pow(d, k)*NH)) : 1;
				}
			}
		}
	}
	d=0.0;
	for (int i = 0; i < sequence_len; ++i)
		{
		/* code */
		x+=conformation->coordinates[i][0];
		
		y+=conformation->coordinates[i][1];
		}
		px=x/sequence_len;
		py=y/sequence_len;

		int filas = 1;
		int columnas = 2;	
		double **C;	

	
		C = (double **)malloc(filas*sizeof(double*)); 
	
		for (i=0;i<filas;i++) 
			C[i] = (double*)malloc(columnas*sizeof(double)); 


	C[0][0] = px; 
	C[0][1] = py; 
	
for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
				
					d = euclidean_distance(conformation->coordinates[i], C[0], dimensions);
					//printf("%f\t%f\t%f\t%f\t%f\n",d ,conformation->coordinates[i][0], conformation->coordinates[i][1],conformation->coordinates[j][0], conformation->coordinates[j][1]);
					energy1 -= (ceil(d)>1.0) ? (1.0/(pow(d, k)*NH)) : 1;
			
		}
	}

	x=px+py;
	conformation->objectives[0] = energy;
	
	conformation->objectives[1] = energy1;



}

void Berenboym2008_energyM(struct conformation_def *conformation){

	int i, k,j;
	double energy=0.0,energy1=0.0, e,e2;
	
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
			for(j=i+2; j<sequence_len; j++){
				if(sequence[j]=='H'){

					if(Hsubset[i] == Hsubset[j]) {
						for(k=0, e=0.0; k<dimensions; k++)
							e+= pow(conformation->coordinates[i][k] - conformation->coordinates[j][k], 2);

						energy-=(1.0/e);
						//printf("%f\n",e );
						e=0;

					}else{

						for(k=0, e=0.0; k<dimensions; k++)
							e2+= pow(conformation->coordinates[i][k] - conformation->coordinates[j][k], 2);

						energy1-=(1.0/e2);
						//printf("%f\n",e2 );
						e2=0;
					}

										
				}
			}
		}
	}

	/*for(i=0; i<conformation->HHtc; i++){
		
		r1 = conformation->contacts[i][0];
		r2 = conformation->contacts[i][1];

		if(Hsubset[r1] == Hsubset[r2]) {
				printf("hola\n");

			for(k=0, e=0.0; k<dimensions; k++)
					e+= pow(conformation->coordinates[r1][k] - conformation->coordinates[r2][k], 2);

			printf("%f\n",e );
			energy-=(1.0/e);
		}
		else {
				printf("ho22la\n");

			for(k=0, e=0.0; k<dimensions; k++)
					e2+= pow(conformation->coordinates[r1][k] - conformation->coordinates[r2][k], 2);

			energy1-=(1.0/e2);
			printf("%f\n",e2 );
		}
	
	}*/

	conformation->objectives[0] = energy;
	conformation->objectives[1] = energy1;

		/*printf("%i\n", conformation->contacts[i][0]);
		printf("%f\t%f\n",conformation->coordinates[r1][0], conformation->coordinates[r1][1]);
		printf("%i\n", conformation->contacts[i][1]);
		printf("%f\t%f\n",conformation->coordinates[r2][0], conformation->coordinates[r2][1]);
		*/
		
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Invokes the evaluation of the given solution with the specified function
void evaluate_conformation_multiobjective(struct conformation_def *conformation){	
	
	if(strcmp(evaluation_strategy, "DEC1") == 0)	 		DEC1(conformation);		
	else if(strcmp(evaluation_strategy, "DEC2") == 0)		DEC2(conformation);	
	else if(strcmp(evaluation_strategy, "DEC3") == 0)		DEC3(conformation);	
	else if(strcmp(evaluation_strategy, "DEC4") == 0)		DEC2(conformation); // DEC4=DEC2, but n_subsets=3
	else if(strcmp(evaluation_strategy, "DEC5") == 0)		DEC5(conformation);	
	else if(strcmp(evaluation_strategy, "PARITY") == 0)	 PARITY(conformation);	
	else if(strcmp(evaluation_strategy, "LOCALITY") == 0)	LOCALITY(conformation);	
	else if(strcmp(evaluation_strategy, "RESIDUE") == 0)	RESIDUE(conformation);	
	else if(strcmp(evaluation_strategy, "RESIDUE2") == 0)	RESIDUE2(conformation);	
	else if(strcmp(evaluation_strategy, "DYNRES2") == 0)	DYNRES(conformation);	
	else if(strcmp(evaluation_strategy, "DYNRES3") == 0)	DYNRES(conformation);	

	else if(strcmp(evaluation_strategy, "K99M") == 0)	Krasnogor1999_energy_Multy(conformation);
	else if(strcmp(evaluation_strategy, "K99M2") == 0)	Krasnogor1999_energy_Multy2(conformation);
	else if(strcmp(evaluation_strategy, "B08M") == 0)	Berenboym2008_energyM(conformation);	

	else if(strcmp(evaluation_strategy, "MOK99_LOCALITY") == 0)	MOK99_LOCALITY(conformation);	
	else if(strcmp(evaluation_strategy, "MOK99_PARITY") == 0)		MOK99_PARITY(conformation);	
	
	else if(strcmp(evaluation_strategy, "I09_M") == 0)	Islam2009_energy_1(conformation);	
	
	else if(strcmp(evaluation_strategy, "I09a") == 0)	Islam2009_multiobjective(conformation, 0);	
	else if(strcmp(evaluation_strategy, "I09b") == 0)	Islam2009_multiobjective(conformation, 1);	
	else if(strcmp(evaluation_strategy, "I09c") == 0)	Islam2009_multiobjective(conformation, 2);	
	else if(strcmp(evaluation_strategy, "I09d") == 0)	Islam2009_multiobjective(conformation, 3);

		

}
