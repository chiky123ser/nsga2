////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Protein Structure Prediction (PSP) in the HP Lattice Model
	
	(psp_generate_validate.c): Random generation and validation of conformations
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	15/02/2011
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the absolute coordinates of the amino acids with respect to 
// the the given range of internal coordinates [from, to]
void update_coordinates(struct conformation_def *c, int from, int to){
	int i, j;
	
	for(i=from; i<=to; i++){
		for(j=0; j<dimensions; j++){
			c->coordinates[i+1][j] = c->coordinates[i][j] + directions[c->absolute_encoding[i]][j];
			//~ c->coordinates[i][j] = c->coordinates[i-1][j] + directions[c->absolute_encoding[i-1]][j];
		}
	}	
}
////////////////////////////////////////////////////////////////////////// //////////////////////////////////////////////
// Generates a random solution accepting only valid solutions
void generate_valid_random_conformation(struct conformation_def *conformation){
	int i, d, rnddir[ndirections], valid, ncoord;
	
	// Initial permutation of directions
	for(i=0; i<ndirections; i++) rnddir[i] = i;
	
	do{

		// First move without restrictions
		conformation->absolute_encoding[0] = rndpsp.RandInt(0, ndirections-1); // Select a random direction		
		//~ conformation->absolute_encoding[0] = 0; // Force the same initial direction
		update_coordinates(conformation, 0, 0); // Compute coordinates
		//~ update_coordinates(conformation, 1, 1); // Compute coordinates
		ncoord=2; // Increase coordinates computed
		valid=1;		

		for(i=1; i<encoding_len && valid==1; i++){		
			
			// Get a permutation of directions
			shuffle(rnddir, ndirections);
			
			// Explore the directions in the obtained order until a valid moved is found	
			valid=0; // Asume is not valid
			for(d=0; d<ndirections && valid==0; d++){
				
				conformation->absolute_encoding[i] = rnddir[d]; // Select a random direction
				update_coordinates(conformation, i, i);	// Compute coordinates
				//~ update_coordinates(conformation, i+1, i+1);	// Compute coordinates
				if(d==0) ncoord++; // Increase coordinates computed
				
				valid = isValid_partial(conformation, ncoord); // Verify
				
				// NOTE: It is computationally expensive to verify the whole solution. 
				// Alternatively, it can be verified only the new candidate move.
			}				
		}
		
	}while(!valid);	
	
	conformation->valid = 1;
} 

//funcion que evalua si el random devuelto no se encuentra ya en la lista de indices randoms
int isValidTemp(int intTemp,int *cr){
		for (int i = 0; i < encoding_len; ++i)
		{
			if (cr[i]==intTemp)
			{
				
				return 1;
			}
		}
		
	return 0;
}

int isValidTemp2(int intTemp,int *cr,int lrn){
	for (int i = 0; i < lrn; ++i)
		{
			if (cr[i]==intTemp)
			{
				return 1;
			}
		}
	return 0;
}
 
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
				//printf("valid: %i\n",valid);
				if(!valid) {
					mutated_solution->absolute_encoding[i] = original;
					update_coordinates(mutated_solution, i, encoding_len-1);
				}else{
					changed = 1;
				}
			}
		}
		
	}while(changed == 0 );
}


void cruce(struct conformation_def *current_solution, struct conformation_def *current_solution2,struct conformation_def *hijo1, struct conformation_def *hijo2){
	
	int *cr,temp=0,valid1,valid2,intTemp=0;
	// Make a copy of the conformation
	copy_conformation(current_solution, hijo1);
	copy_conformation(current_solution2, hijo2);
	

	cr = (int *)malloc(sizeof(int) * 8);
	
	for (int i = 0; i < 8; ++i)
	{	
		intTemp=rndpsp.RandInt(0, encoding_len-1);
		
		if (isValidTemp2(intTemp,cr,8)==1)
		{
			--i;
		}
		else{
			cr[i]=intTemp;
			
		}	
	}
	
	

	
	for (int i = 0; i < 8; ++i)
	{
		temp=hijo1->absolute_encoding[cr[i]];
		hijo1->absolute_encoding[cr[i]]=hijo2->absolute_encoding[cr[i]];
		hijo2->absolute_encoding[cr[i]]=temp;	
	}
		
	

		
	update_coordinates(hijo1,0,encoding_len-1);
	update_coordinates(hijo2,0,encoding_len-1);
		

		
	valid1=isValid(hijo1);
	valid2=isValid(hijo2);
		
	if (!valid1)
	{
		copy_conformation(current_solution, hijo1);
	}

		if (!valid2)
	{
		copy_conformation(current_solution2, hijo2);
	}
	

	free(cr);

}


void cruce_Cycle(struct conformation_def *current_solution, struct conformation_def *current_solution2,struct conformation_def *hijo1, struct conformation_def *hijo2){
	//void cruce_Cycle(){
	//int i, j, rnddir[ndirections], original, valid, changed = 0;	
	//double p_mut = (1.0 / encoding_len);
 

	int *cr,*cr1,intTemp=0,j=0,temp=0,valid1=0,valid2=0;

	// Make a copy of the conformation
	
	copy_conformation(current_solution, hijo1);
	copy_conformation(current_solution2, hijo2);
	
	cr = (int *)malloc(sizeof(int) * encoding_len);
	cr1 = (int *)malloc(sizeof(int) * encoding_len);

	for (int i = 0; i < encoding_len; ++i)
	{
		cr[i]=i;
	
	}


	for (int i = 0; i < encoding_len; ++i)
	{
		cr1[i]=encoding_len+1;
	}
	
	for (int i = 0; i < encoding_len; ++i)
	{
		intTemp=rndpsp.RandInt(0, encoding_len-1);
		
		if (isValidTemp(intTemp,cr1)==1)
		{
			--i;
			
		}
		else{
			cr1[i]=intTemp;	
		}	
	}

	while(cr1[j]!=cr[0]){
		temp=hijo1->absolute_encoding[cr[j]];
		
		hijo1->absolute_encoding[cr[j]]=hijo2->absolute_encoding[cr[j]];
		
		hijo2->absolute_encoding[cr[j]]=temp;	
		
		j=cr1[cr[j]];

	}

	update_coordinates(hijo1,0,encoding_len-1);
	update_coordinates(hijo2,0,encoding_len-1);

	   /* printf("current_solution\n");
	    print_conformation(current_solution);
	    printf("hijo1\n");
	    print_conformation(hijo1);
	    printf("current_solution2\n");
	    print_conformation(current_solution2);
	    printf("hijo2\n");
	    print_conformation(hijo2);*/



	valid1=isValid(hijo1);
	valid2=isValid(hijo2);

	
	

	if (!valid1)
	{
		copy_conformation(current_solution, hijo1);
		
	}

		if (!valid2)
	{
		copy_conformation(current_solution2, hijo2);

	}

	free(cr);
	free(cr1);
	cr1=NULL;
	cr=NULL;
	
}



void NON_WRAPPING_ORDER_CROSSOVER(struct conformation_def *current_solution, struct conformation_def *current_solution2,struct conformation_def *hijo1, struct conformation_def *hijo2){

	

	int techo=0,piso=0,posicion=0,valid1,valid2,intTemp;

	// Make a copy of the conformation 
	copy_conformation(current_solution, hijo1);
	copy_conformation(current_solution2, hijo2);
	
	
	techo=ceil(encoding_len*.20);
	
	piso=floor(encoding_len/2);
	posicion=ceil(piso-(techo/2));


	int *crm,*crm2,*cri,*cri2,*crf,*crf2,*ran,*ran2;

	crm = (int *)malloc(sizeof(int) * techo);
	crm2 = (int *)malloc(sizeof(int) * techo);
	ran = (int *)malloc(sizeof(int) * techo);
	ran2 = (int *)malloc(sizeof(int) * techo);
	cri = (int *)malloc(sizeof(int) * posicion);
	cri2 = (int *)malloc(sizeof(int) * posicion);
	crf= (int *)malloc(sizeof(int) * (encoding_len-(posicion+techo)));
	crf2= (int *)malloc(sizeof(int) * (encoding_len-(posicion+techo)));

	
	
	for (int i = 0; i < techo; ++i)
	{	

			intTemp=rndpsp.RandInt(0, encoding_len-1);
		if (intTemp>=posicion && intTemp<=posicion+techo)
		{
			--i;
		}
		else
		{
			/* code */
			if (isValidTemp2(intTemp,ran,techo)==1)
			{
				--i;
			}
			else{
				ran[i]=intTemp;
				
			}
		}

	}

	for (int i = 0; i < techo; ++i)
	{	

			intTemp=rndpsp.RandInt(0, encoding_len-1);
		if (intTemp>=posicion && intTemp<=posicion+techo)
		{
			--i;
		}
		else
		{
			/* code */
			if (isValidTemp2(intTemp,ran2,techo)==1)
			{
				--i;
			}
			else{
				ran2[i]=intTemp;
				
			}
		}

	}


	for (int i = 0; i < techo; ++i)
	{
		hijo1->absolute_encoding[ran[i]]=100;
		hijo2->absolute_encoding[ran2[i]]=100;

	}
	
	

	for (int i = 0; i < posicion; ++i)
	{	
			cri[i]=hijo1->absolute_encoding[i];
			cri2[i]=hijo2->absolute_encoding[i];
			//printf("%c\t",letter[cri[j]] );
			//printf("%c\n", letter[cri2[j]]);
	}


	for (int i = posicion,j=0; i < posicion+techo; ++i,++j)
	{	

			crm[j]=hijo1->absolute_encoding[i];
			crm2[j]=hijo2->absolute_encoding[i];
	
	}

	for (int i = posicion+techo,j=0; i < encoding_len; ++i,++j)
	{	

			crf[j]=hijo1->absolute_encoding[i];
			crf2[j]=hijo2->absolute_encoding[i];
			//printf("%c\t",letter[crf[j]] );
			//printf("%c\n", letter[crf2[j]]);
	}

	int temp,temp2,flag=1;
	while(flag==1){
		flag=0;
		for (int i = 0; i < posicion; ++i)
		{
			if (cri[i]==100 && i+1!=posicion && cri[i+1]!=100 )
			{
				temp=cri[i];
				cri[i]=cri[i+1];
				cri[i+1]=temp;
				flag=1;
			}
			if (cri2[i]==100 && i+1!=posicion && cri2[i+1]!=100 )
			{
				temp2=cri2[i];
				cri2[i]=cri2[i+1];
				cri2[i+1]=temp2;
				flag=1;
			}
		}

	}




	flag=1;
	while(flag==1){
		flag=0;
		for (int i = (encoding_len-(posicion+techo)); i > 0; --i)
		{
			if (crf[i]==100 && i-1!=posicion && crf[i-1]!=100 )
			{
				temp=crf[i];
				crf[i]=crf[i-1];
				crf[i-1]=temp;
				flag=1;
			}
			if (crf2[i]==100 && i-1!=posicion && crf2[i-1]!=100 )
			{
				temp2=crf2[i];
				crf2[i]=crf2[i-1];
				crf2[i-1]=temp2;
				flag=1;
			}
		}

	}


	int k=0,l=0;

	for (int i = 0; i < posicion; ++i)
	{
		if (cri[i]==100)
		{
			cri[i]=crm[k];
			k++;
		}
		if (cri2[i]==100)
		{
			cri2[i]=crm2[l];
			l++;
		}
	}




	for (int i = 0; i <(encoding_len-(posicion+techo)); ++i)
	{
		if (crf[i]==100)
		{
			crf[i]=crm[k];
			k++;
		}
		if (crf2[i]==100)
		{
			crf2[i]=crm2[l];
			l++;
		}
	}

	k=0;l=0;
	

	for (int i = 0; i < encoding_len; ++i)
	{
		if (i<posicion)
		{
			hijo1->absolute_encoding[i]=cri[i];
			hijo2->absolute_encoding[i]=cri2[i];
		}
		else if(i>=posicion && i<posicion+techo){
			hijo1->absolute_encoding[i]=crm2[k];
			
			hijo2->absolute_encoding[i]=crm[k];
			k++;
			
		}
		else{
			
			hijo1->absolute_encoding[i]=crf[l];
		
			hijo2->absolute_encoding[i]=crf2[l];
			l++;
		}

	}




	update_coordinates(hijo1,0,encoding_len-1);
	update_coordinates(hijo2,0,encoding_len-1);

	valid1=isValid(hijo1);
	valid2=isValid(hijo2);

	if (!valid1)
	{
		
		copy_conformation(current_solution, hijo1);
		
	}

		if (!valid2)
	{
		
		copy_conformation(current_solution2, hijo2);

	}
		
	free(crm);
	free(crm2);
	free(ran);
	free(ran2);
	free(cri);
	free(cri2);
	free(crf);
	free(crf2);

}

void cruce_static(struct conformation_def *current_solution, struct conformation_def *current_solution2,struct conformation_def *hijo1, struct conformation_def *hijo2){

	int *cr,temp=0,valid1,valid2,punto;
	// Make a copy of the conformation

	punto=floor(encoding_len/8);
	copy_conformation(current_solution, hijo1);

	copy_conformation(current_solution2, hijo2);
	
	cr = (int *)malloc(sizeof(int) * 8);

	int i=0,j=0;
	while(i!=8){
		
		j=j+punto;
		cr[i]=j;
		
		++i;
	}

	
	


	for (int i = 0; i < 8; ++i)
	{
		temp=hijo1->absolute_encoding[cr[i]];
		hijo1->absolute_encoding[cr[i]]=hijo2->absolute_encoding[cr[i]];
		hijo2->absolute_encoding[cr[i]]=temp;	
	}
		
	

	update_coordinates(hijo1,0,encoding_len-1);
	update_coordinates(hijo2,0,encoding_len-1);

	valid1=isValid(hijo1);
	valid2=isValid(hijo2);

	if (!valid1)
	{
		copy_conformation(current_solution, hijo1);
	}

		if (!valid2)
	{
		copy_conformation(current_solution2, hijo2);
	}

	free(cr);

}

void cruce_one_point(struct conformation_def *current_solution, struct conformation_def *current_solution2,struct conformation_def *hijo1, struct conformation_def *hijo2){
	
	int temp=0,valid1,valid2,intTemp=0;
	// Make a copy of the conformation

	copy_conformation(current_solution, hijo1);
	copy_conformation(current_solution2, hijo2);
	
	intTemp=rndpsp.RandInt(0, encoding_len-1);
	

	
	for (int i = intTemp; i < encoding_len; ++i)
	{
		temp=hijo1->absolute_encoding[i];
		hijo1->absolute_encoding[i]=hijo2->absolute_encoding[i];
		hijo2->absolute_encoding[i]=temp;	
	}
		
	

		
	update_coordinates(hijo1,0,encoding_len-1);
	update_coordinates(hijo2,0,encoding_len-1);
		

		
	valid1=isValid(hijo1);
	valid2=isValid(hijo2);
		
	if (!valid1)
	{
		copy_conformation(current_solution, hijo1);
		
	}

		if (!valid2)
	{
		copy_conformation(current_solution2, hijo2);
		
	}

}

void cruce_two_point(struct conformation_def *current_solution, struct conformation_def *current_solution2,struct conformation_def *hijo1, struct conformation_def *hijo2){
	
	int temp=0,valid1,valid2,intTemp1=0,intTemp2=0;
	// Make a copy of the conformation

	copy_conformation(current_solution, hijo1);
	copy_conformation(current_solution2, hijo2);
	
	intTemp1=rndpsp.RandInt(0, encoding_len-1);
	intTemp2=rndpsp.RandInt(0, encoding_len-1);
	
	while(intTemp1>=intTemp2){
		
		intTemp1=rndpsp.RandInt(0, encoding_len-1);
		intTemp2=rndpsp.RandInt(0, encoding_len-1);
	}

	
	for (int i = intTemp1; i <= intTemp2; ++i)
	{
		temp=hijo1->absolute_encoding[i];
		hijo1->absolute_encoding[i]=hijo2->absolute_encoding[i];
		hijo2->absolute_encoding[i]=temp;	
	}
		
	

		
	update_coordinates(hijo1,intTemp1,intTemp2);
	update_coordinates(hijo2,intTemp1,intTemp2);
		

		
	valid1=isValid(hijo1);
	valid2=isValid(hijo2);
		
	if (!valid1)
	{
		copy_conformation(current_solution, hijo1);
		
	}

		if (!valid2)
	{
		copy_conformation(current_solution2, hijo2);
		
	}

	
		
}

