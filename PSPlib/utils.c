////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Protein Structure Prediction (PSP) in the HP Lattice Model
	
	(utils.c)  	Generic utility functions
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	15/02/2011
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Receives the vector of values and the elements' ids
// order=0: Ascending    order=1: Descending
void quicksort(double *value, int *ids, long int inf, long int sup, int order){
	double elem_div = value[sup], temp;
	int tmp, i = inf - 1, j = sup, cont = 1;

	if (inf >= sup)   return;

	while (cont){
		if(order==0){
			while (i< j && value[++i] < elem_div);
			while (j>i && value[--j] > elem_div);
		}else{		
			while (i< j && value[++i] > elem_div);
			while (j>i && value[--j] < elem_div);
		}
		
		if (i < j){
			temp = value[i];
			value[i] = value[j];
			value[j] = temp;
			tmp=ids[i];
			ids[i]=ids[j];
			ids[j]=tmp;
		}else
			cont = 0;
	}	

	temp = value[i];
	value[i] = value[sup];
	value[sup] = temp;	
	tmp=ids[i];
	ids[i]=ids[sup];
	ids[sup]=tmp;
	quicksort (value, ids, inf, i - 1, order);
	quicksort (value, ids, i + 1, sup, order);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Does not receive the ids vector
// order=0: Ascending    order=1: Descending
void quicksortWID(double *value, long int inf, long int sup, int order){ 
	double elem_div = value[sup], temp;
	int i = inf - 1, j = sup, cont = 1;

	if (inf >= sup)   return;

	while (cont){
		if(order==0){
			while (i< j && value[++i] < elem_div);
			while (j>i && value[--j] > elem_div);
		}else{		
			while (i< j && value[++i] > elem_div);
			while (j>i && value[--j] < elem_div);
		}
		
		if (i < j){
			temp = value[i];
			value[i] = value[j];
			value[j] = temp;
		}else
			cont = 0;
	}	

	temp = value[i];
	value[i] = value[sup];
	value[sup] = temp;	
	quicksortWID (value, inf, i - 1, order);
	quicksortWID (value, i + 1, sup, order);
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Shuffles (generates a permutation of) the elements of a given integer vector
void shuffle(int *vector, int size){
	int i, temp, r;	
	for(i=0; i<size; i++){
		r = rndpsp.RandInt(i, size-1);
		temp = vector[i];
		vector[i] = vector[r];
		vector[r] = temp;
	}
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares two double vectors. Returns 1 if equal, 0 otherwise
int isEqual_double(double *v1, double *v2, int size){
	int i;	
	for(i=0; i<size; i++){
		if(v1[i] != v2[i]) return 0;
	}
	return 1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares two integer vectors. Returns 1 if equal, 0 otherwise
int isEqual_int(int *v1, int *v2, int size){
	int i;	
	for(i=0; i<size; i++){
		if(v1[i] != v2[i]) return 0;
	}
	return 1;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the minimum of  the given double vector
double min_double(double *list, int size){
	int i, min=0;
	
	for(i=1; i<size; i++)
		if(list[i] < list[min]) min = i;	
	
	return list[min];
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Returns the minimum of  the given integer vector
int min_int(int *list, int size){
	int i, min=0;
	
	for(i=1; i<size; i++)
		if(list[i] < list[min]) min = i;	
	
	return list[min];
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the Manhattan distance (1-norm) of the given pair of vectors
double manhattan_distance(double *vector1, double *vector2, int size){
	int i;
	double dist=0.0;
	for(i=0; i<size; i++) dist+= fabs(vector1[i]-vector2[i]);
	return dist;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Computes the euclidean distance (2-norm) of the given pair of vectors
double euclidean_distance(double *vector1, double *vector2, int size){
	int i;
	double dist=0.0;
	
	for(i=0; i<size; i++) dist+= (vector1[i]-vector2[i])*(vector1[i]-vector2[i]); 
	return sqrt(dist);	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hamming distance: the number of different positions
int hamming_distance(int *v1, int *v2, int size){
	int distance=0, i;
	for(i=0; i<size; i++)
		if( v1[i] != v2[i]) distance++;
	return distance;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Makes a copy of an array of doubles
 void copy_vector_double(double *from, double *to, int size){
	int i;
	for(i=0; i<size; i++) to[i] = from[i];	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Makes a copy of an array of integers
 void copy_vector_int(int *from, int *to, int size){
	int i;
	for(i=0; i<size; i++) to[i] = from[i];	
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int *getVector(int n){
  int *m;
  m=(int *)malloc(n*sizeof(int));
  if (!m){fprintf(stderr,"Vector Memory error!\n");exit(-1);}
  return m;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int **getMatrix(int a,int b){
  int i,**m;
  m=(int **)malloc(a*sizeof(int *));
  if (!m){fprintf(stderr,"Matriz Memory error in step 1!\n");exit(-1);}
  for (i=0;i<a;i++){
    m[i]=(int *)malloc(b*sizeof(int));
    if (!m[i]){fprintf(stderr,"Matrix Memory error in step 2!\n");exit(-1);}
  }
  return m;
}
