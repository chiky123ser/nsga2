////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Protein Structure Prediction (PSP) in the HP Lattice Model
	
	(psp.h): 	This the main file, to be included as a library in the 
			implementation of search algorithms
	
	Author: 	Mario Garza Fabre / CINVESTAV-Tamaulipas
	Contact: 	mgarza@tamps.cinvestav.mx
			garzafabre@gmail.com
	
	Updated:	15/02/2011
*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Standard libraries
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>

// Random numbers generator
#include "RngStream.cpp"

// PSP framework files
#include "psp_global.h"
#include "utils.c"
#include "psp_std.c"
#include "psp_validate.c"
#include "psp_generate.c"
#include "psp_evaluate_singleobjective.c"
#include "psp_evaluate_multiobjective.c"
#include "psp_evaluate.c"
#include "psp_multiobjectivization.c"



