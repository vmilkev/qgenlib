/*******************************************************************************
!   Copyright(C) 1999-2012 Intel Corporation. All Rights Reserved.
!   
!   The source code, information  and  material ("Material") contained herein is
!   owned  by Intel Corporation or its suppliers or licensors, and title to such
!   Material remains  with Intel Corporation  or its suppliers or licensors. The
!   Material  contains proprietary information  of  Intel or  its  suppliers and
!   licensors. The  Material is protected by worldwide copyright laws and treaty
!   provisions. No  part  of  the  Material  may  be  used,  copied, reproduced,
!   modified, published, uploaded, posted, transmitted, distributed or disclosed
!   in any way  without Intel's  prior  express written  permission. No  license
!   under  any patent, copyright  or  other intellectual property rights  in the
!   Material  is  granted  to  or  conferred  upon  you,  either  expressly,  by
!   implication, inducement,  estoppel or  otherwise.  Any  license  under  such
!   intellectual  property  rights must  be express  and  approved  by  Intel in
!   writing.
!   
!   *Third Party trademarks are the property of their respective owners.
!   
!   Unless otherwise  agreed  by Intel  in writing, you may not remove  or alter
!   this  notice or  any other notice embedded  in Materials by Intel or Intel's
!   suppliers or licensors in any way.
!
!*******************************************************************************
!  Content:
!      D H A M M  Example Program Text ( C Interface )
!******************************************************************************/

#include "mkl.h"

#ifndef min
#define min(a,b) ((a) <= (b) ? (a) : (b) )
#endif


// The threshold problem size for the base case
#ifndef NB
#define NB 512
#endif


// The interface of the DHAMM routine
void dhamm( char *uplo, char *transa, char *transb, MKL_INT *m, 
            MKL_INT *n, MKL_INT *k, double *alpha, double *A, 
            MKL_INT *lda, double *B, MKL_INT *ldb, double *beta, 
            double *C, MKL_INT *ldc )
{
   dhamm_recur_(uplo, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
}

// The recursive routine called by DHAMM. It divides an upper or lower
// triangular matrix into two smaller triangular parts and one rectangular part.
// It computes the rectangular part using DGEMM, and then recursively computes
// the two smaller triangular parts, respectively.
void dhamm_recur_( char *uplo, char *transa, char *transb, MKL_INT *m, 
                   MKL_INT *n, MKL_INT *k, double *alpha, double *A, 
                   MKL_INT *lda, double *B, MKL_INT *ldb, double *beta, 
                   double *C, MKL_INT *ldc )
{
   double *Aptr, *Bptr, *Cptr; 
   MKL_INT i, j, trA, trB, m1, n1, m2, n2;
   extern void dhamm_naive_();

   // The base case, solved using the naive algorithm
   if ( min(*m,*n) <= NB ) 
   {
      dhamm_naive_ ( uplo, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);
      return ; 
   }

   trA = ( *transa=='T' || *transa=='t' );
   trB = ( *transb=='T' || *transb=='t' );
   if ( *uplo=='L' || *uplo=='l' ) // If the lower triangular part of C is needed
   {
      m1 = *m / 2;
      n1 = min(m1,*n);
      Aptr = A; 
      Bptr = B;
      Cptr = C;
      
      // Start recursion on one smaller triangular part
      dhamm_recur_( uplo, transa, transb, &m1, &n1, k, alpha, Aptr, lda, Bptr, ldb, beta, Cptr, ldc);

      m2 = *m - m1;
      if ( trA ) Aptr = &A[m1*(*lda)]; else Aptr = &A[m1];
      Bptr = B;
      Cptr = &C[m1];

      // Compute the rectangular part
      dgemm ( transa, transb, &m2, &n1, k, alpha, Aptr, lda, Bptr, ldb, beta, Cptr, ldc);

      n2 = *n - n1;
      if ( trA ) Aptr = &A[m1*(*lda)]; else Aptr = &A[m1];
      if ( trB ) Bptr = &B[n1]; else Bptr = &B[n1*(*ldb)];
      Cptr = &C[m1+n1*(*ldc)];

      // Start recursion on the other smaller triangular part
      dhamm_recur_( uplo, transa, transb, &m2, &n2, k, alpha, Aptr, lda, Bptr, ldb, beta, Cptr, ldc);
   } else {                         // If the upper triangular part of C is needed
      m1 = *m / 2;
      n1 = min(m1,*n);
      Aptr = A; 
      Bptr = B;
      Cptr = C;
      
      // Start recursion on one smaller triangular part
      dhamm_recur_( uplo, transa, transb, &m1, &n1, k, alpha, Aptr, lda, Bptr, ldb, beta, Cptr, ldc);

      m2 = *m - m1;
      n2 = *n - n1;
      Aptr = A;
      if ( trB ) Bptr = &B[n1]; else Bptr = &B[n1*(*ldb)];
      Cptr = &C[n1*(*ldc)];

      // Compute the rectangular part
      dgemm ( transa, transb, &m1, &n2, k, alpha, Aptr, lda, Bptr, ldb, beta, Cptr, ldc);

      if ( trA ) Aptr = &A[m1*(*lda)]; else Aptr = &A[m1];
      Cptr = &C[m1+m1*(*ldc)];

      // Start recursion on the other smaller triangular part
      dhamm_recur_( uplo, transa, transb, &m2, &n2, k, alpha, Aptr, lda, Bptr, ldb, beta, Cptr, ldc);
   }
}

// The naive algorithm for the base case: Call DGEMM and then discard the upper
// triangular elements if only the lower triangular part is needed. Likewise,
// discard the lower triangular elements if only the upper triangular part is
// needed.
void dhamm_naive_ ( char *uplo, char *transa, char *transb, MKL_INT *m, 
                    MKL_INT *n, MKL_INT *k, double *alpha, double *A, 
                    MKL_INT *lda, double *B, MKL_INT *ldb, double *beta, 
                    double *C, MKL_INT *ldc )
{
   double *D, *Cp=&C[0], *Dp;
   MKL_INT i, j;
   const double d_zero = 0.0, d_one = 1.0;
  
   if ( *m <= 0 || *n <= 0 ) return ;
   D = (double *) mkl_malloc ( (*m)*(*n)*sizeof(double), 64 );
   Dp = &D[0];

   // Call DGEMM
   dgemm ( transa, transb, m, n, k, alpha, A, lda, B, ldb, &d_zero, D, m);   

   // Discard half of the result, depending on if lower or upper triangular part
   // is needed.
   if ( *uplo=='L' || *uplo=='l' ) 
   {
        for ( i = j ; i < *m ; i++ ) Cp[i] = *beta*Cp[i] + Dp[i];
   } else {
        for ( i = 0 ; i <= min(j,*m-1) ; i++ ) Cp[i] = *beta*Cp[i] + Dp[i];
   }

   mkl_free(D);
}


