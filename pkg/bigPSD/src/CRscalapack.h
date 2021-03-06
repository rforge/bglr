/*======================================================================
 *           R-ScaLAPACK version 0.4.x:  ScaLAPACK interface to R
 *              Oak Ridge National Laboratory, Oak Ridge TN.
 *        Authors: David Bauer, Guruprasad Kora, Nagiza. F. Samatova, 
 *                            Srikanth Yoginath.
 *     Contact: Nagiza F. Samatova; (865) 241-4351; samatovan@ornl.gov
 *     Contact: Guruprasad Kora; (865) 576-6210; koragh@ornl.gov
 *                 Computer Science and Mathematics Division
 *             Oak Ridge National Laboratory, Oak Ridge TN 37831 
 *                   (C) 2004 All Rights Reserved
 *
 *                              NOTICE
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted
 * provided that the above copyright notice appear in all copies and
 * that both the copyright notice and this permission notice appear in
 * supporting documentation.
 *
 * Neither the Oak Ridge National Laboratory nor the Authors make any
 * representations about the suitability of this software for any
 * purpose.  This software is provided ``as is'' without express or
 * implied warranty.
 *
 * RScaLAPACK (http://www.aspect-sdm.org/Parallel-R) was funded
 * as part of the Scientific Data Management Center
 * (http://sdm.lbl.gov/sdmcenter) under the Department of Energy's 
 * Scientific Discovery through Advanced Computing (DOE SciDAC) program
 * (http://www.scidac.org ). 
=========================================================================*/
#ifndef _RSCALAPACK_H_
#define _RSCALAPACK_H_
#define USE_RINTERNALS
#include  "mpi.h"
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <R_ext/PrtUtil.h>
#include <R_ext/RS.h>

#include <dlfcn.h>

#define DONT_SPAWN_R

#if 0
#define DEBUG_RSCALAPACK
#endif

#define PC_SEND_DATA_DIM 501

#define CR_RecvVectorFromPA crrecvvectorfrompa_ 
#define CR_SendVectorToPA crsendvectortopa_

#define CR_SendIntToPA crsendinttopa_
#define CR_SendDoubleToPA crsenddoubletopa_
#define CR_CheckFailFlag crcheckfailflag_

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))

void F77_NAME(blacs_pinfo)(int *pid, int *nprocs);
void F77_NAME(blacs_get)(int *,int *, int *);
void F77_NAME(blacs_exit)(int *ipExitFlag);
int F77_NAME(numroc)(int *, int *, int *, int *, int *);
int F77_NAME(indxg2p)(int *, int *, int *, int *, int *);

void F77_NAME(callpdgesv) (int *, double *, int *);
void F77_NAME(callpdgesvd)(int *, double *, int *);
void F77_NAME(callpdpotrf)(int *, double *, int *);
void F77_NAME(callpdgeqrf)(int *, double *, int *);
void F77_NAME(callpdsyevd)(int *, double *, int *);
void F77_NAME(callpdpotri)(int *, double *, int *);
void F77_NAME(callpdgemm) (int *, double *, int *);
void F77_NAME(callpdgetcp)(int *, double *, int *);

SEXP CR_Exec();
int CR_InitializeEnv(int *myrank, int *nprocs);
int CR_CallScalapackFn(int *, int);
int CR_GetInputParams(MPI_Comm parent, int *dim);

void CR_RecvVectorFromPA (int *ib, int *ia, double *A, int *mb);
void CR_SendVectorToPA (int *ib, int *ia, double *work, int *mb);

void CR_SendIntToPA(int *, int *, int *);
void CR_SendDoubleToPA(double *, int *, int *);

int CRSF_qr (int dim[], int myRank);
int CRSF_chol (int dim[], int myrank);
int CRSF_chol2inv (int dim[], int myrank);
int CRSF_eigen (int dim[], int myrank);
int CRSF_solve (int dim[], int myrank);
int CRSF_svd (int dim[], int myrank);
int CRSF_multiply (int dim[], int myrank);
int CRSF_tcrossprod(int dim[], int myrank);

#endif
