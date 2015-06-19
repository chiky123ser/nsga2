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
// Returns 1 if the given points are adjacent in the lattice, 0 otherwise
int isAdjacent(double *c1, double *c2){
	
	// Option 1 
	// Two points c1 and c2 are adjacent in the lattice if the distance 
	// between them is 1. This is correct for the square and cubic lattices, 
	// but what happends in the case of, for example, the FCC lattice?
	
	if(ceil(manhattan_distance(c1, c2, dimensions)) == 1) return 1;	
	return 0;

	// Option 2 - More expensive
	// Two points c1 and c2 are adjacent if the vector v = (c1 - c2) 
	// is one of the direction vectors in the lattice definition
	
	//~ int i;
	//~ double v[dimensions];	
	//~ // Compute vector v = (c1 - c2) 
	//~ for(i=0; i<dimensions; i++) v[i] = (c1[i] - c2[i]);	
	//~ // Determine if v is one of the direction vectors of the lattice
	//~ for(i=0; i<ndirections; i++) if(isEqual_double(v, directions[i], dimensions)) return 1;
	//~ return 0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the number of H-H topological contacts in a conformation
void getHHTopologicalContacts(struct conformation_def *conformation){
	int i, j, HHtc=0;
	
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
			for(j=i+2; j<sequence_len; j++){
				//~ if(sequence[j]=='H') HHtc += isAdjacent(conformation->coordinates[i], conformation->coordinates[j]);
				
				if(sequence[j]=='H'){
					
					if( isAdjacent(conformation->coordinates[i], conformation->coordinates[j]) ){
						
						// Store contacts
						conformation->contacts[HHtc][0] = i;
						conformation->contacts[HHtc][1] = j;
						/*printf("%i\n", conformation->contacts[HHtc][0]);
						printf("%f\t%f\n",conformation->coordinates[conformation->contacts[HHtc][0]][0], conformation->coordinates[conformation->contacts[HHtc][0]][1]);
						printf("%i\n", conformation->contacts[HHtc][1]);
						printf("%f\t%f\n",conformation->coordinates[conformation->contacts[HHtc][1]][0], conformation->coordinates[conformation->contacts[HHtc][1]][1]);
						*/
						// Increase contacts counter
						HHtc += 1;						
					}					
				}
			}
		}
	}	
	

	conformation->HHtc = HHtc;
	conformation->PPtc = 0;
	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Conventional HP model energy function proposed by K. Dill in 1985
// The energy of a structure is given by the negative of the number of HH contacts
// HH contact occurs when non-consecutive H amino acids are adjacent in the lattice.
void Dill1985_energy(struct conformation_def *conformation){

	conformation->objectives[0] = (-1.0*conformation->HHtc);
	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
H. Li, R. Helling, C. Tang, and N. Wingreen. Emergence of Preferred Structures in a Simple
Model of Protein Folding. Science, 273(5275):666-669, August 1996.

eHH=-2.3, eHP=-1, ePP=0
*/
void Li1996_energy(struct conformation_def *conformation){
	int i, j, a;
	double energy=0.0;
	
	for(i=0; i<sequence_len-2; i++){					
		for(j=i+2; j<sequence_len; j++){
			if(sequence[i]=='H'){
				a = isAdjacent(conformation->coordinates[i], conformation->coordinates[j]);				
				energy -= (sequence[j]=='H')? 2.3*a : a;				
			}else{
				if(sequence[j]=='H')
					energy -= isAdjacent(conformation->coordinates[i], conformation->coordinates[j]);
			}
		}		
	}
	
	conformation->objectives[0] = energy;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Distance-based potential proposed in:

N. Krasnogor, W.E. Hart, J. Smith, and D.A. Pelta. Protein Structure Prediction With Evo-
lutionary Algorithms. In Wolfgang Banzhaf, Jason M. Daida, A. E. Eiben, Max H. Garzon,
Vasant Honavar, Mark J. Jakiela, and Robert E. Smith, editors, Proceedings of the Genetic
and Evolutionary Computation Conference (GECCO 1999), Orlando, Florida, USA, July 1999.
Morgan Kaufman.
*/
void Krasnogor1999_energy(struct conformation_def *conformation){
	int i, j, k;
	double energy=0.0, d;
	
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
	
	conformation->objectives[0] = energy;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Given a lattice location, returns 1 if the provided neighboring location
// is an H residue, 2 if it is a P residue, and 3 if is an empty location
int locationStatus(struct conformation_def *conformation, double *location){
	int i;
	
	// Verify if the provided location is one of the coordinates used by the solution
	// In positive case, verify if this location is occupied by an H or P residue.
	for(i=0; i<sequence_len; i++){		
		if(isEqual_double(location, conformation->coordinates[i], dimensions)) {
			return (sequence[i]=='H') ? 1 : 2;
		}		
	}	
	return 3; // Empty location, i.e. HS contact
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
F.L. Custodio, H.J.C. Barbosa, and L.E. Dardenne. Investigation of the three-dimensional
lattice HP protein folding model using a genetic algorithm. Genetics and Molecular Biology,
27:611-615, 2004.
*/
void Custodio2004_energy(struct conformation_def *conformation){
	int i,j,k, HHc=0, HPc=0, HSc=0, s;
	double location[dimensions], w1=0.0, w2=10.0, w3=40.0;
	
	// Explore the neighboring locations of each H residue
	for(i=0; i<sequence_len; i++){
		if(sequence[i]=='H'){		
			
			// Verify all directions
			for(j=0; j<ndirections; j++){
				
				// Get location
				for(k=0; k<dimensions; k++) location[k] = conformation->coordinates[i][k] + directions[j][k];
				
				// Verify location
				s = locationStatus(conformation, location);
				
				if(s==1) HHc++;
				else if(s==2) HPc++;
				else HSc++;
			}			
		}
	}	
	
	conformation->objectives[0] = (w1*HHc + w2*HPc + w3*HSc);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Energy function based on the concept of Radius of Gyration. Proposed in:

H. Lopes and M. Scapin. An Enhanced Genetic Algorithm for Protein Structure Prediction
Using the 2D Hydrophobic-Polar Model. In El-Ghazali Talbi, Pierre Liardet, Pierre Collet,
Evelyne Lutton, and Marc Schoenauer, editors, Artificial Evolution, volume 3871 of Lecture
Notes in Computer Science, pages 238-246. Springer Berlin / Heidelberg, 2006.
*/
void Lopes2006_energy(struct conformation_def *conformation){
	int i, j, HnLB=0;//, NC=0;
	double meanHcoord[dimensions], meanPcoord[dimensions], meanHmax[dimensions]; //PW
	double RadiusH=0.0, RadiusP=0.0, RgH=0.0, RgP=0.0, MaxRgH=0.0, diffRG;
	
	// Compute the number of HH topological contacts (non-local bonds)
	HnLB = conformation->HHtc;
	
	// Penalization according to the number of collisions
	//~ NC=0; // Collisions are not taken into account in this study
	//~ PW = (0.033*sequence_len)+1.33;
	//~ HnLB -= (NC*PW);
		
	// Compute mean coordinates of H and P residues	
	for(j=0; j<dimensions; j++)  meanHmax[j]=meanHcoord[j]=meanPcoord[j]=0.0;	

	for(i=0; i<sequence_len; i++){
		if(sequence[i]=='H'){
			for(j=0; j<dimensions; j++){ 
				meanHcoord[j] += conformation->coordinates[i][j];
				
				// The totally unfolded conformation is assumed to be
				// composed only by moves in the first direction of the
				// lattice specification
				meanHmax[j] += directions[0][j]*i;
			}
		}else{
			for(j=0; j<dimensions; j++) meanPcoord[j] += conformation->coordinates[i][j];
		}
	}
	
	for(j=0; j<dimensions; j++) {
		meanHcoord[j] /= NH;
		meanPcoord[j] /= NP;
		meanHmax[j] /= NH;
	}
	
	// Compute RgH and RgP, and MaxRgH
	for(i=0; i<sequence_len; i++){
		if(sequence[i]=='H'){			
			for(j=0; j<dimensions; j++){ 
				RgH += pow((conformation->coordinates[i][j] - meanHcoord[j]), 2);	
				
				// The totally unfolded conformation is assumed to be
				// composed only by moves in the first direction of the
				// lattice specification				
				MaxRgH += pow((directions[0][j]*i - meanHmax[j]), 2);			
			}
		}else{
			for(j=0; j<dimensions; j++) 
				RgP += pow((conformation->coordinates[i][j] - meanPcoord[j]), 2);
		}
	}
	
	RgH = sqrt(RgH/(double)NH);
	RgP = sqrt(RgP/(double)NP);
	MaxRgH = sqrt(MaxRgH/(double)NH);		
	
	// Compute RadiusH term
	RadiusH = MaxRgH - RgH;
	
	// Compute RadiusP term
	diffRG = (RgP - RgH);	
	RadiusP = (diffRG >= 0) ? 1.0 : (1.0 / (1.0-diffRG));
	
	// Final energy
	conformation->objectives[0] = -(HnLB * RadiusH * RadiusP);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
Global Energy proposed in:

I. Berenboym and M. Avigal. Genetic algorithms with local search optimization for protein
structure prediction problem. In GECCO '08: Proceedings of the 10th annual conference on
Genetic and evolutionary computation, pages 1097-1098, New York, NY, USA, 2008. ACM.
*/
void Berenboym2008_energy(struct conformation_def *conformation){
	int i, j, k;
	double energy=0.0, e;
	
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
			for(j=i+2; j<sequence_len; j++){
				if(sequence[j]=='H'){
					
					for(k=0, e=0.0; k<dimensions; k++)
						e+= pow(conformation->coordinates[i][k] - conformation->coordinates[j][k], 2);
					
					energy -= (1.0/e);						
				}
			}
		}
	}
	
	conformation->objectives[0] = energy;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
M. Cebrian, I. Dotu, P. Van Hentenryck, and P. Clote. Protein structure prediction on the
face centered cubic lattice by local search. In Proceedings of the 23rd national conference on
Artificial intelligence - Volume 1, pages 241-246. AAAI Press, 2008.
*/
void Cebrian2008_energy(struct conformation_def *conformation){
	int i, j, k, l;
	double energy=0.0, distance;
	
	// The authors tested with k={1,2,3}, k=2 performed the best
	k = 2;
	
	for(i=0; i<sequence_len-2; i++){
		if(sequence[i]=='H'){			
			for(j=i+2; j<sequence_len; j++){
				if(sequence[j]=='H'){
					
					for(l=0, distance=0.0; l<dimensions; l++){
						distance += pow(conformation->coordinates[i][l] - conformation->coordinates[j][l], 2);					
					}				
					energy += pow((distance - 1.0), k);	
				}
			}
		}
	}
	
	conformation->objectives[0] = energy;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
K. Islam and M. Chetty. Novel Memetic Algorithm for Protein Structure Prediction. In Ann
Nicholson and Xiaodong Li, editors, AI 2009: Advances in Artificial Intelligence, volume 5866
of Lecture Notes in Computer Science, pages 412-421. Springer Berlin / Heidelberg, 2009.
*/
void Islam2009_energy(struct conformation_def *conformation){
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
	conformation->objectives[0] = alpha*HP + Hcompliance + Pcompliance;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Invokes the evaluation of the given solution with the specified function
void evaluate_conformation_singleobjective(struct conformation_def *conformation){
	
	// Evaluate the conformation by using the selected energy function
	//functions = ["D85", "L96", "K99", "C04", "L06", "B08", "C08", "I09"]
	
	if(strcmp(evaluation_strategy, "D85") == 0)	 	Dill1985_energy(conformation);		
	else if(strcmp(evaluation_strategy, "L96") == 0)	Li1996_energy(conformation);	
	else if(strcmp(evaluation_strategy, "K99") == 0)	Krasnogor1999_energy(conformation);	
	else if(strcmp(evaluation_strategy, "C04") == 0)	Custodio2004_energy(conformation);	
	else if(strcmp(evaluation_strategy, "L06") == 0)	Lopes2006_energy(conformation);	
	else if(strcmp(evaluation_strategy, "B08") == 0)	Berenboym2008_energy(conformation);	
	else if(strcmp(evaluation_strategy, "C08") == 0)	Cebrian2008_energy(conformation);	
	else if(strcmp(evaluation_strategy, "I09") == 0)	Islam2009_energy(conformation);	

}

