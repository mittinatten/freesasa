#ifndef INTEGRATE_H
#define INTEGRATE_H

#include <stdio.h>
#include "atomclassifier.h"
#include "protein.h"

/**  */
void integrate_sasa_per_atomclass(FILE*, atomclassifier*, 
				  protein*, double *sasa);


#endif
