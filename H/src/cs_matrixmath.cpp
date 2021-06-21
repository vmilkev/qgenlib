/*
 * cs_matrixmath.cpp
 *
 *  Created on: Sep 7, 2017
 *      Author: vimi
 */

#include "cs_matrixmath.hpp"


//===============================================================================================================

MATRIXMATH::MATRIXMATH() {
	// TODO Auto-generated constructor stub

}

//===============================================================================================================

int MATRIXMATH::make_pivot_array (lapack_int *&ar, lapack_int row) {

	try
	{
		/*
		 * Prepare the initialized vector which is used during matrix inversion.
		 */

		// allocate memory
		ar = (lapack_int *)mkl_malloc( row*sizeof( lapack_int ), 64 );
		if (ar == NULL) {
			mkl_free(ar);
			throw 40;
		}
		//initialize by '1'
		for (lapack_int i = 0; i < (row); i++) {
			ar[i] = 1;
		}
	}
	catch (int ex) {
		return ex;
	}
	catch (std::exception const& e) {
		return 1;
	}

	return 0;
}

//===============================================================================================================

int MATRIXMATH::matr_prod (double *A, double *B, double *C, MKL_INT rowA, MKL_INT rowB, MKL_INT colA, MKL_INT colB) {

	/*
	 * Matrix product:
	 * C = A * B;
	 * 
	 * void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,...
	 * 					const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,...
	 * 					const MKL_INT k, const double alpha, const double *a,...
	 * 					const MKL_INT lda, const double *b, const MKL_INT ldb,...
	 * 					const double beta, double *c, const MKL_INT ldc);
	 */

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, colB, colA, 1.0, A, colA, B, colB, 0.0, C, colB);

	return 0;
}

//===============================================================================================================

int MATRIXMATH::matr_prodTr (double *A, double *B, double *C, MKL_INT rowA, MKL_INT rowB, MKL_INT colA, MKL_INT colB) {

	/*
	 * Matrix product:
	 * C = A' * B;
	 * 
	 * void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,...
	 * 					const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,...
	 * 					const MKL_INT k, const double alpha, const double *a,...
	 * 					const MKL_INT lda, const double *b, const MKL_INT ldb,...
	 * 					const double beta, double *c, const MKL_INT ldc);
	 */

	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, rowA, colB, colA, 1.0, A, rowA, B, colB, 0.0, C, colB);

	return 0;
}

//===============================================================================================================

int MATRIXMATH::matr_add (double *A, double *B, double *C, bool transpA, bool transpB, size_t row, size_t col) {

	/*
	 * matrix addition :
	 * C = A + B;
	 */

	char trA, trB;
	size_t lda, ldb, ldc;

	if (transpA) {
		trA = 't';
		lda = row;
	}
	else {
		trA = 'n';
		lda = col;
	}
	if (transpB) {
		trB = 't';
		ldb = row;
	}
	else {
		trB = 'n';
		ldb = col;
	}
	ldc = col;

	mkl_domatadd ('r', trA, trB, row, col, 1.0, A, lda, 1.0, B, ldb, C, ldc);

	return 0;
}

//===============================================================================================================

int MATRIXMATH::matr_sub (double *A, double *B, double *C, bool transpA, bool transpB, size_t row, size_t col) {

	/*
	 * matrix subtraction :
	 * C = A + B;
	 */

	char trA, trB;
	size_t lda, ldb, ldc;

	if (transpA) {
		trA = 't';
		lda = row;
	}
	else {
		trA = 'n';
		lda = col;
	}
	if (transpB) {
		trB = 't';
		ldb = row;
	}
	else {
		trB = 'n';
		ldb = col;
	}
	ldc = col;

	mkl_domatadd ('r', trA, trB, row, col, 1.0, A, lda, -1.0, B, ldb, C, ldc);

	return 0;
}

//===============================================================================================================

int MATRIXMATH::matr_inv (double *A, lapack_int *ipiv, lapack_int row, lapack_int col) {

	lapack_int info = 0;

	try
	{
		/*
		 * Matrix inversion.
		 * Here:
		 * 		A - matrix to invert;
		 * 		row - number of rows in A;
		 * 		col - number of columns in A;
		 * 		ipiv - empty vector of size lda.
		 *
		 * 1. First, make LU factorizatioon of A: A_fact = P*L*U; ipiv = P.
		 * 2. Second, make inversion of factorized A: A = inv(A_fact).
		 * The result of inversion is in matrix A.
		 *
		 */

		int matrix_order = LAPACK_ROW_MAJOR;

		//auto start = std::chrono::high_resolution_clock::now();

		info = LAPACKE_dgetrf (matrix_order, row, col, A, col, ipiv);
		if (info != 0) throw info;

/*
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Elapsed time of LU decomposition: " << elapsed.count() << " s\n";
*/

		//start = std::chrono::high_resolution_clock::now();

		info = LAPACKE_dgetri (matrix_order, row, A, row, ipiv);
		if (info != 0) throw info;

/*
		finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;
		std::cout << "Elapsed time of LU-based inverse: " << elapsed.count() << " s\n";
*/
	}
	catch (lapack_int ex) {
		return ex;
	}
	catch (std::exception const& e) {
		return 1;
	}

	return 0;
}

//===============================================================================================================

int MATRIXMATH::matr_inv_h (double *A, lapack_int lda) {

	lapack_int info = 0;;

	try
	{
		int matrix_order = LAPACK_ROW_MAJOR;

		//auto start = std::chrono::high_resolution_clock::now();

		info = LAPACKE_dpptrf (matrix_order, 'L', lda, A);

/*
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << "Elapsed time of Cholesky decomposition: " << elapsed.count() << " s\n";
*/

		if (info != 0) {
			std::cout<<"error in Cholesky decomposition: "<<info<<std::endl;
			throw info;
		}

		//start = std::chrono::high_resolution_clock::now();

		info = LAPACKE_dpptri (matrix_order, 'L', lda, A);

/*
		finish = std::chrono::high_resolution_clock::now();
		elapsed = finish - start;
		std::cout << "Elapsed time of symmetric inverse: " << elapsed.count() << " s\n";
*/

		if (info != 0) {
			throw info;
			std::cout<<"error in inverse: "<<info<<std::endl;
		}

	}
	catch (lapack_int ex) {
		return ex;
	}
	catch (std::exception const& e) {
		return 1;
	}

	return 0;
}

//===============================================================================================================

MATRIXMATH::~MATRIXMATH() {
	// TODO Auto-generated destructor stub
}

//===============================================================================================================


