/**
 * 
 * @author Petrucci Enrico
*/


#ifndef preprocessing_H_  
#define preprocessing_H_


#include "scp.h"


int preprocessinglegacycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                         const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);
void performPreprocessing(CPXENVptr env, CPXLPptr lp, instance* inst);

int scpPreprocessing(instance *inst);


#endif   /* preprocessing_H_ */
