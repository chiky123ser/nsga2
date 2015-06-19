////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Protein Structure Prediction (PSP) in the HP Lattice Model
	
	(psp_global.h):  Global declarations
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	15/02/2011
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CONSTANTS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define PI 3.14159265358979
#define EPS 1.0e-14
#define INF 1.0e14

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// LATTICE CONFIGURATION 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int dimensions;				// Dimensions of the lattice (2D or 3D)
int ndirections;				// Number of direction vectors
double **directions;			// Direction vectors
int letters;					// Indicate if  letter equivalences are specified for the lattice vectors
char *letter;					// Letter equivalences

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PROBLEM SETTINGS 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

char *sequence;				// HP Sequence to be solved
int sequence_len;				// Length of the sequence 
int optimalHHtc;				// Best known number of HH top. contacts
int NH, NP;					// Number of H and P residues in the sequence

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CONFORMATIONS ENCODING
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct conformation_def{		// Definition of the CONFORMATION data type
	int *absolute_encoding;		// Internal coordinates encoding - absolute moves
	//~ int *relative_encoding;		// Internal coordinates encoding - relative moves
	double **coordinates;		// Absolute coordinates of each amino acid
	//~ int ncoord;				// Number of  coordinates computed (useful in the case of partial solutions)
	int HHtc;					// Number of H-H topological contacts
	int PPtc;
	//~ double energy;			// Energy (objective) value 
	double *objectives;			// Objective values (to be minimized). Single-objective implementations will refer the energy value to as objectives[0]
	int **contacts;				// Stores the pairs of residues forming an HH contact
	int **contactsP;				// Stores the pairs of residues forming an HH contact
	int valid;					// 1 if valid, 0 otherwise
	int rankpsp;					// Rank of the individual = {1:best, N:worst} (used by some algorithms)
	double cwd;				// Crowding distance
	unsigned long evaluation;	// Total number of evaluations until this solution was reached
};
int encoding_len;				// Length of the encoded solution

//~ int absolute_encoding_len;				// Length of the encoded solution
//~ int relative_encoding_len;				// Length of the encoded solution

struct conformation_def evaluated_best;

/* 
	A conformation or a group of conformations is to be defined as:
		struct conformation_def *conformation;	
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// OTHER VARIABLES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int run;							// Number of execution (only for reports and for the seed of the randoms number generator)

unsigned long seed; 				// Seed for random numbers generator
RngStream rndpsp("global stream"); 	// Stream for randoms number generator

unsigned long evaluation=0;			// Function evaluations counter
unsigned long max_evaluations=0;	// Allowed number of objective function evaluations
unsigned long report_each=50000;	// Each N evaluations the results are reported

//~ unsigned long iteration=0;			// Function iterations counter
//~ unsigned long max_iterations=0;		// Allowed number of iterations

char problem_formulation[20]; 		// SO: single objective / MO: Multiobjective
char evaluation_strategy[20]; 			// Evaluation strategy to be used (energy function or multiobjectivization approach)
char strategy_variant[20]; 				// Type of subsets to use: DET, RND, DYN

char acceptance_criterion[20]; 		

int n_objectives = 1;				// Number of objective functions to be used

int n_subsets=1;	// Number of Hsubsets to organize Hs
int *Hsubset;	// For each H in the sequence, store the corresponding subset

int localitysize;	// Max distance between amino acids to be considered neighbors (local interaction)

int DYNRES_OBJ[5] = {0, 1, 2, 1, 0}; // This initialization corresponds to RESIDUE2 decomposition
int DYNRES_MAXHH = 0;	// Each residue could be involved in at most 5 contacts in 3D (3 contacts in 2D)

int iterations_without_move=0;
int max_iwm = -1;

int online_energy=1;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FUNCTION PROTOTYPES
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

