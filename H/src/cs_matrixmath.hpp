/*
 * cs_matrixmath.hpp
 *
 *  Created on: Sep 7, 2017
 *      Author: vimi
 */
#pragma once

#include "structures.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>
#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <sys/resource.h>

#define MKL_INT size_t
#include "mkl.h"

#include <chrono>
#include <ctime>

//#ifndef min
//	#define min(x,y) (((x) < (y)) ? (x) : (y))
//#endif

class MATRIXMATH {

public:

	MATRIXMATH();
	virtual ~MATRIXMATH();

	int matr_prod (double *A, double *B, double *C, MKL_INT rowA, MKL_INT rowB, MKL_INT colA, MKL_INT colB);
	int matr_prodTr (double *A, double *B, double *C, MKL_INT rowA, MKL_INT rowB, MKL_INT colA, MKL_INT colB);
	int matr_add (double *A, double *B, double *C, bool transpA, bool transpB, size_t row, size_t col);
	int matr_sub (double *A, double *B, double *C, bool transpA, bool transpB, size_t row, size_t col);
	int matr_inv (double *A, lapack_int *ipiv, lapack_int row, lapack_int col);
	int matr_inv_h (double *A, lapack_int lda);

	template <typename T>
	int make_res_array (T *&ar, size_t row, size_t col);

	template <typename T>
	int make_res_array_h (T *&ar, size_t lda);

	int make_pivot_array (lapack_int *&ar, lapack_int row);

private:

};

//===============================================================================================================

template <typename T>
int MATRIXMATH::make_res_array_h (T *&ar, size_t lda) {

	try
	{
		/*
		 * Prepare the zero initialized matrix 'ar' which will hold results of other two matrices operations:
		 * ar = a + b; ar = a - b; ar = a * b; ar <- a (insert a into bigger ar).
		 */

		// check the data type to allocate memory: 32 or  64
		int dType = sizeof(T)*8;

		size_t sz = static_cast<size_t>((lda*lda + lda)/2);

		// allocate memory
		ar = (T *)mkl_malloc( sz*sizeof( T ), dType );

		if (ar == NULL) {
			mkl_free(ar);
			throw 40;
		}

		//initialize by '0.0'
		for (size_t i = 0; i < sz; i++) {
			ar[i] = 0.0;
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

template <typename T>
int MATRIXMATH::make_res_array (T *&ar, size_t row, size_t col) {

	try
	{
		/*
		 * Prepare the zero initialized matrix 'ar' which will hold results of other two matrices operations:
		 * ar = a + b; ar = a - b; ar = a * b; ar <- a (insert a into bigger ar).
		 */

		// check the data type to allocate memory: 32 or  64
		auto dType = sizeof(T)*8;

		// allocate memory
		ar = (T *)mkl_malloc( row*col*sizeof( T ), dType );
		if (ar == NULL) {
			mkl_free(ar);
			throw 40;
		}
		//initialize by '0.0'
		for (size_t i = 0; i < (row*col); i++) {
			ar[i] = 0.0;
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

