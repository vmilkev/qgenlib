/*
    cs_matrix.hpp

    General matrix class intended to work with Intel® Math Kernel Library
    in terms of matrix storage formats, memory functions and inversion & multiplication routines.

*/

#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <algorithm>
#include <utility>
#include <iterator>
#include <cstring>

#define MKL_INT size_t
#include "mkl.h"

#ifndef _min
    #define _min(x,y) (((x) < (y)) ? (x) : (y))
#endif

#ifndef _max
    #define _max(x,y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef worksize
    #define worksize 10000
#endif

namespace qgen
{

    template <typename T>
    class matrix
    {

        public:

            /* CONSTRUCTORS & DESTRUCTOR */

            matrix(size_t row, size_t col);                /* Constructor for rectangular matrix. */
            matrix(size_t lda);                            /* Constructor for square symmetric half-store (compact) matrix. */
            matrix(const char *type);                      /* Default constructor: rectangular (type = 'r'|| 'R')/symmetric (type = 's'|| 'S') matrix, no memmory allocated. */
            matrix();                                      /* Default constructor: rectangular matrix, no memmory allocated. */
            matrix( const matrix& obj);                    /* Copy constructor. */
            ~matrix();                                     /* Destructor. */


            /* METHODS */

            T & at(size_t atRow, size_t atCol);            /* Get/put element from/in a matrix. */
            void resize(size_t row, size_t col);           /* Resizes/allocates memmory for A. */
            void resize(size_t lda);                       /* Resizes/allocates memmory for A; overloaded method for symmetrical matrix. */
            void print(std::string whiichMatrix);         /* Prints part of a matrix into a LOG file. */
            void scale(T val);                             /* Scaling matrix by scalar: A = A*val. */
            size_t size();                                 /* Gives total number of elements in a matrix. */
            size_t capacity();                             /* Gives total number of allocated elements in a matrix. */
            void clear();                                  /* Frees an allocated memmory of A. */
            bool empty();                                  /* Checks if memmory is allocated to A. */
            void symtorec();                               /* Transform symmetric matrix in compact form to rectangular form. */
            void rectosym();                               /* Transform square matrix to triangular compact form (only for symmetric matrices). */
            void transpose();                              /* Transpose matrix. */
            void fwrite();                                 /* Move matrix to the disk and clear memory. */
            void invert();                                 /* Matrix inversion. */
            void fread();                                  /* Restore matrix from the disk into the memory. */
            bool eq(const matrix& rhs);                    /* Compare dimensins and shapes of two matrix objects. */


            /* OPERATORS */

            matrix operator+(const matrix& rhs);           /* Overloaded '+' operator to add two matrix objects. */
            matrix operator-(const matrix& rhs);           /* Overloaded '-' operator to substract two matrix objects. */
            matrix operator-(const T val);                 /* Overloaded '-' operator to substitute a scalar from a matrix object. */
            matrix operator^(const int val);               /* Overloaded '^' operator to multiply matrix by itself and find inversion. */
            matrix operator^(const char *val);             /* Overloaded '^' operator to transpose matrix. */
            matrix operator*(const matrix& rhs);           /* Overloaded '*' operator to multiply two matrix objects. */
            bool operator==(const matrix& rhs);            /* Compare complete equality of two matrix objects. */
            matrix & operator=(const matrix& rhs);         /* Overloaded assignment '=' operator. */
            T & operator()(size_t atRow, size_t atCol);    /* Access element of a matrix, memory reallocation is allowed. */
            T operator()(size_t atRow, size_t atCol) const;/* Access element of a matrix. */
            T & operator[](size_t i);                      /* Access element of a matrix. */
            T operator[](size_t i) const;                  /* Access element of a matrix. */


            /* VARIABLES & CONSTANTS */

            bool failbit;                                  /* Error state flag. */
            int failinfo;


        private:

            /* METHODS */

            int allocate (size_t row, size_t col);         /* Allocate memory for a rectangular matrix. */	    
            int allocate (size_t lda);                     /* Allocate memory for a half-store (compact) symmetric matrix. */
            void resize();                                 /* Resizes memmory allocated for A. */
            unsigned long long rdtsc();                    /* Seed for random number generator. */
            
            /* Interfaces to MKL routines */

            void dotprod (double *A, double *B, double *C, MKL_INT rowA, MKL_INT rowB, MKL_INT colA, MKL_INT colB);
            void dotprod (float *A, float *B, float *C, MKL_INT rowA, MKL_INT rowB, MKL_INT colA, MKL_INT colB);
            void inv_rec (double *A, MKL_INT rowA, MKL_INT colA);
            void inv_rec (float *A, MKL_INT rowA, MKL_INT colA);
            void inv_sym (double *A, MKL_INT colA);
            void inv_sym (float *A, MKL_INT colA);
            void gemmt_intrf (double *A, double *B, MKL_INT rowA, MKL_INT colA, MKL_INT colB);
            void gemmt_intrf (float *A, float *B, MKL_INT rowA, MKL_INT colA, MKL_INT colB);
            

            /* VARIABLES & CONSTANTS */

            T *A;                                          /* The main matrix container. */
            size_t numRow;                                 /* Number of rows in A. */
            size_t numCol;                                 /* Number of columns in A. */
            size_t resizedElements;                        /* Number of allocated/resized elements in A. */
            bool allocated;                                /* Flag indicated that the memory for A have been allocated. */
            bool rectangular;                              /* TRUE if the matrix is rectangular. */
            bool symetric;                                 /* TRUE if the matrix is symmetric. */
            bool compact;                                  /* TRUE if for the symmetric matrix only a lower triangular part is stored. */		
            std::string debug_file = "MATRIX.log";
            std::string binFilename;                       /* Name of binary file to store A on disck. */
            std::fstream fA;


        public:

        /* FRIEND OPERATORS */

        /*
         Such operator's definitions are inlined in the class definition because operators don't exist outside the class definition
         until a template instantiation of the class is generated (which happens during the compilation process).
         Here we perform implicit type conversion on non-templated non-member functions.
         We create the non-template friend functions which are automatically created for each template instantiation of the class.
        */
        friend matrix operator*(const T val, const matrix& rhs){
            /* 
                Overloaded '*' operator to multiply matrix by scalar.

                Return value: matrix object.
                
                Example:

                    qgen::matrix <double> obj;    // empty matrix
                    qgen::matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                    for (auto i = 0; i < M.size(); i++)
                        M[i] = 1.0;
                    
                    obj = -2.0 * M;  // obj become (n,m) matrix where all elements are -2.0;
                                     // M remains unchanged (all elements of the matrix are 1.0)
            */

            return rhs * val;
        }

        friend matrix operator*(const matrix& lhs, const T val){
            /* 
                Overloaded '*' operator to multiply matrix by scalar.

                Return value: matrix object.
                
                Example:

                    qgen::matrix <double> obj;    // empty matrix
                    qgen::matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.

                    for (auto i = 0; i < M.size(); i++)
                        M[i] = 1.0;
                    
                    obj = M * (-2.0); // obj become (n,m) matrix where all elements are -2.0;
                                      // M remains unchanged (all elements of the matrix are 1.0)
            */

            matrix <T> C;
            if( !lhs.compact ){
                int status = C.allocate(lhs.numRow, lhs.numCol);
                if(status != 0) {
                    C.failbit = true;
                    throw std::string("Memory allocation error");
                }
                
                C.allocated = true;
                C.resizedElements = lhs.numRow*lhs.numCol;
                C.numCol = lhs.numCol;
                C.numRow = lhs.numRow;
            }
            else{
                int status = C.allocate(lhs.numRow);
                if(status != 0) {
                    C.failbit = true;
                    throw std::string("Memory allocation error");
                }
                
                C.allocated = true;
                C.resizedElements = (lhs.numRow*lhs.numRow+lhs.numRow)/2;
                C.numCol = lhs.numRow;
                C.numRow = lhs.numRow;
            }

            auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int> (C.size()/(n_threads));
            
            if(block_size < worksize){
                block_size = C.size();
                n_threads = 1;
            }

        #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < C.size(); i++)
                C.A[i] = lhs.A[i] * val;

            return matrix(C);
        }

        friend matrix operator+(const T lhs, const matrix& rhs){
            /* 
                Overloaded '+' operator to add a scalar to a matrix.

                Return value: matrix object.
                
                Example:

                    qgen::matrix <double> obj;    // empty matrix
                    qgen::matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.
                    
                    obj = -2.0 + M; // obj become (n,m) matrix where all elements are -2.0; M remains unchanged
            */

            return rhs + lhs;
        }

        friend matrix operator+(const matrix& lhs, const T rhs){
            /* 
                Overloaded '+' operator to add a scalar to a matrix.

                Return value: matrix object.
                
                Example:

                    qgen::matrix <double> obj;    // empty matrix
                    qgen::matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.
                    
                    obj = M + (-2.0); // obj become (n,m) matrix where all elements are -2.0; M remains unchanged
            */

            matrix <T> C;
            if( !lhs.compact ){
                int status = C.allocate(lhs.numRow, lhs.numCol);
                if(status != 0) {
                    C.failbit = true;
                    throw std::string("Memory allocation error");
                }
                
                C.allocated = true;
                C.resizedElements = lhs.numRow*lhs.numCol;
                C.numCol = lhs.numCol;
                C.numRow = lhs.numRow;
            }
            else{
                int status = C.allocate(lhs.numRow);
                if(status != 0) {
                    C.failbit = true;
                    throw std::string("Memory allocation error");
                }
                
                C.allocated = true;
                C.resizedElements = (lhs.numRow*lhs.numRow+lhs.numRow)/2;
                C.numCol = lhs.numRow;
                C.numRow = lhs.numRow;
            }

            auto n_threads = std::thread::hardware_concurrency();
            auto block_size = static_cast<unsigned int> (C.size()/(n_threads));
            
            if(block_size < worksize){
                block_size = C.size();
                n_threads = 1;
            }

        #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for (size_t i = 0; i < C.size(); i++)
                C.A[i] = lhs.A[i] + rhs;

            return matrix(C);
        }

    };

}

//===============================================================================================================

template <typename T>
unsigned long long qgen::matrix<T>::rdtsc(){
    /* Seed for random number generator. */
    
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    
    return ((unsigned long long)hi << 32) | lo;
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::transpose(){
    /* 
       Handles only rectangular matrices (including full symmetric).
       No transpose of triangular matrix since we do not use 'U' format.
       
       Return value: none; modifies the container A of the calling object.
       
       Example:

            qgen::matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.
            qgen::matrix <double> res;    // empty matrix
            
            for (auto i = 0; i < M.size(); i++)
                M[i] = i;
            
            M.transpose(); // now M is (m,n) matrix
            res = M;       // now res is (m,n) matrix

    */
    if( !compact ){

        matrix<T> C(numCol, numRow);
                
        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int> (size()/(n_threads));

        if(block_size < worksize){
            block_size = size();
            n_threads = 1;
        }
        
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < numRow; i++) {
            for (size_t j = 0; j < numCol; j++) {
                C.A[j*numRow+i] = A[i*numCol+j];
            }
        }
        resize(numCol, numRow);
        *this = C;
        C.clear();
    }
    else{
        /*
            We also have to accept symmetric matrices in compact form.
            Due to the reasons described in the upper comment, we return the same matrix
            but in restored (not compact) format.
        */
        symtorec();
    }
}

//===============================================================================================================

template <typename T>
T & qgen::matrix<T>::operator()(size_t atRow, size_t atCol){
    /*
        Fast access operator (array-like).

        Return value: the element at the specified position in the container.
        
        Example:
                 qgen::matrix <double> M(2,2);
                 M(1,1) = 0.8;
                 double val = M(1,1);   // val = 0.8.
    */

    if (!compact ) return A[atRow*numCol+atCol];
    else return A[atRow*(atRow+1)/2+atCol];
}

//===============================================================================================================

template <typename T>
T qgen::matrix<T>::operator()(size_t atRow, size_t atCol) const {
    /*
        Fast access operator (array-like).

        Return value: the element at the specified position in the container.
        
        Example:
                 qgen::matrix <double> M(2,2);
                 M(1,1) = 0.8;
                 double val = M(1,1);   // val = 0.8.
    */

    if (!compact ) return A[atRow*numCol+atCol];
    else return A[atRow*(atRow+1)/2+atCol];
}


//===============================================================================================================

template <typename T>
T & qgen::matrix<T>::operator[](size_t i){
    /*
        Fast access operator (direct).

        Return value: the element at the specified position in the container.
        
        Example:

            qgen::matrix <double> M(n,m);
            double val = M[ M.size()-1 ]; // val = 0.0
            
            for (auto i = 0; i < M.size(); i++)
                M[ i ] = 1.0;
            
            val = M[ M.size()-1 ]; // val = 1.0

    */

    return A[i];
}

//===============================================================================================================

template <typename T>
T qgen::matrix<T>::operator[](size_t i) const {
    /*
        Fast access operator (direct).

        Return value: the element at the specified position in the container.
        
        Example:

            qgen::matrix <double> M(n,m);
            double val = M[ M.size()-1 ]; // val = 0.0

            for (auto i = 0; i < M.size(); i++)
                M[ i ] = 1.0;

            val = M[ M.size()-1 ]; // val = 1.0

    */

    return A[i];
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::resize(){
    /*
        Private member. Reallocates memmorey for the container A.
        Note: row and column of a matrix have to be modified before calling the resize() method.
        
        Return value: none.
    */

    size_t sz;
    if ( !compact )
        sz = numRow*numCol;
    else{
        size_t lda = _max(numRow, numCol);
        sz = (lda*lda + lda)/2;
    }
    
    A = (T*)mkl_realloc(A,sz*sizeof(T));
    if (A == NULL) {
        mkl_free(A);
        allocated = false;
        failbit = true;
        throw std::string("Memory allocation error");
    }
    resizedElements = sz;
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::resize(size_t row, size_t col){
    /*
        Reallocates memmorey for the member array A.
        
        Return value: none.
        Example:
        
            qgen::matrix <double> M(3,5);
            size_t val;
            val = M.size(); // val = 15.
            M.resize(2,4);
            val = M.size(); // val = 8.
    */

    size_t sz;
    size_t lda;
    if ( !compact ){
        sz = row*col;
    }
    else{
        lda = _max(row, col);
        sz = (lda*lda + lda)/2;
    }

    if (allocated){
        A = (T*)mkl_realloc(A,sz*sizeof(T));
        if (A == NULL) {
            mkl_free(A);
            allocated = false;
            failbit = true;
            throw std::string("Memory reallocation error");
        }
    }
    else{
        if ( !compact ){
            int status = allocate (row, col);
            if (status != 0){
                failbit = true;
                throw std::string("Memory allocation error");
            }
        }
        else{
            int status = allocate (lda);
            if (status != 0){
                failbit = true;
                throw std::string("Memory allocation error");
            }
        }
        allocated = true;
    }
    resizedElements = sz;
    if( !compact ){
        numRow = row;
        numCol = col;
    }
    else{
        numRow = numCol = _max(row, col);
    }
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::resize(size_t lda){
    /*
        Reallocates memmorey for the member array A.
        Overloaded method for symmetrical matrix in compact form.
        
        Return value: none.
        Example:
        
            qgen::matrix <double> M(3);
            size_t val;
            val = M.size(); // val = 6.
            M.resize(4);
            val = M.size(); // val = 10.
    */

    size_t sz = (lda*lda + lda)/2;

    if (allocated){
        A = (T*)mkl_realloc(A,sz*sizeof(T));
        if (A == NULL) {
            mkl_free(A);
            allocated = false;
            failbit = true;
            throw std::string("Memory reallocation error");
        }
    }
    else{
        compact = true;
        int status = allocate (lda);
        if (status != 0){
            failbit = true;
            throw std::string("Memory allocation error");
        }
    allocated = true;
    compact = true;
    }

    resizedElements = sz;
    numRow = numCol = lda;
}

//===============================================================================================================

template <typename T>
int qgen::matrix<T>::allocate (size_t lda) {
    /*
        Private member.
        Allocates memmory for a symmetric matrix in compact form.

        Return value: integer value; if 0 - allocation is successfull, otherwise - error.
    */

    size_t sz = static_cast<size_t>((lda*lda + lda)/2);
    A = (T *)mkl_malloc( sz*sizeof( T ), sizeof(T)*8 );
    if (A == NULL) {
        mkl_free(A);
        allocated = false;
        failbit = true;
        throw std::string("Memory allocation error");
    }
    for (size_t i = 0; i < sz; i++) {
        A[i] = 0.0;
    }

	return 0;
}

//===============================================================================================================

template <typename T>
int qgen::matrix<T>::allocate (size_t row, size_t col) {
    /*
        Private member.
        Allocates memmory for an arbitrary rectangular matrix.

        Return value: integer value; if 0 - allocation is successfull, otherwise - error.
    */

    A = (T *)mkl_malloc( row*col*sizeof( T ), sizeof(T)*8 );
    if (A == NULL) {
        mkl_free(A);
        allocated = false;
        failbit = true;
        throw std::string("Memory allocation error");
    }
    for (size_t i = 0; i < (row*col); i++) {
        A[i] = 0.0;
    }

	return 0;
}

//===============================================================================================================

template <typename T>
qgen::matrix<T>::matrix(size_t row, size_t col){
    /* 
       Constructor for regular rectangular matrix
    */

    failbit = false;
    failinfo = 0;

    srand ( rdtsc() );
    int iNum = rand() % 100000;
    binFilename = "matrix_"+std::to_string(iNum);

    rectangular = true;
    compact = false;
    symetric = false;
    allocated = true;

    int status = allocate(row, col);
    if (status != 0){
        failbit = true;
        allocated = false;
        failinfo = -1;
        throw std::string("Memory allocation error");
    }
    
    numRow = row;
    numCol = col;
    resizedElements = row*col;
}

//===============================================================================================================

template <typename T>
qgen::matrix<T>::matrix(size_t lda){
    /* 
       Constructor for symmetrical matrix in a compact form
    */

    failbit = false;
    failinfo = 0;

    srand ( rdtsc() );
    int iNum = rand() % 100000;
    binFilename = "matrix_"+std::to_string(iNum);

    rectangular = false;
    symetric = true;
    allocated = true;
    compact = true;

    int status = allocate (lda);
    if (status != 0){
        failbit = true;
        allocated = false;
        failinfo = -1;
        throw std::string("Memory allocation error");
    }
    
    numRow = numCol = lda;
    resizedElements = (lda*lda+lda)/2;
}

//===============================================================================================================

template <typename T>
qgen::matrix<T>::matrix(){
    /* 
       Default constructor for rectangular matrix;
       no memmory allocation.
    */

    failbit = false;
    failinfo = 0;

    srand ( rdtsc() );
    int iNum = rand() % 100000;
    binFilename = "matrix_"+std::to_string(iNum);

    allocated = false;
    rectangular = true;
    symetric = false;
    compact = false;
    numRow = numCol = 0;
    resizedElements = 0;
}

//===============================================================================================================

template <typename T>
qgen::matrix<T>::matrix(const char *type){
    /* 
       Constructor: explicite determination of a matrix type;
       no memmory allocation.
    */

    failbit = false;
    failinfo = 0;

    srand ( rdtsc() );
    int iNum = rand() % 100000;
    binFilename = "matrix_"+std::to_string(iNum);

    allocated = false;
    if (*type == 's' || 'S') {
        rectangular = false;
        symetric = true;
        compact = true;
    }
    else if (*type == 'r' || 'R'){
        rectangular = true;
        symetric = false;
        compact = false;
    }
    numRow = numCol = 0;
    resizedElements = 0;
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::scale(T val){
	/*
        Matrix scaling by scalar.

        Return value: none.

        Example:

            qgen::matrix <double> M(3);
            double val;
            for(auto i = 0; i < M.size(); i++)
                M[i] = 2.0;

            M.scale(0.5);
            val = M[ 0 ]; // val = 1.0
            val = M[ 1 ]; // val = 1.0
            ...
            val = M[ M.size()-1 ]; // val = 1.0

	*/

    if (allocated){

        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int> (size()/(n_threads));

        if(block_size < worksize){
            block_size = size();
            n_threads = 1;
        }

    #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < size(); i++)
            A[i] = A[i] * val;
    }
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::symtorec(){
    /*
        Trnsforms symmetric matrix in compact form to regular rectangular matrix (not compact).

        Return value: none.

        Example:
                qgen::matrix <mytype> M(3); // symmetric matrix in compact form
                qgen::matrix <mytype> res; // empty matrix
                M.symtorec(); // M now is regular square (symmetric) matrix
                res = M; // res now is regular square (symmetric) matrix
    */

    if ( allocated && compact ){
        
        matrix <T> C(numRow, numRow);

        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int> (size()/(n_threads));

        if(block_size < worksize){
            block_size = size();
            n_threads = 1;
        }

    #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < numRow; i++) {
            for (size_t j = 0; j <= i; j++) {
                C.A[i*numRow+j] = C.A[j*numRow+i] = A[i*(i + 1)/2 + j];
            }
        }

        compact = false;
        resize(numRow, numRow);
        *this = C;
    }
    else{
        failbit = true;
        if (!allocated)
            throw std::string("Matrix is empty");
        else
            throw std::string("Matrix is not in a compact form");
            
    }        
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::rectosym(){
    /*
        Trnsforms regular rectangular (symmetric) matrix into compact form.

        Return value: none.

        Example:
                qgen::matrix <mytype> M(3); // symmetric matrix in compact form
                qgen::matrix <mytype> res; // empty matrix
                M.symtorec();
                res = M; // res now is regular square (symmetric) matrix
                res.rectosym(); // res now is symmetric matrix in a compact form
    */

    if ( allocated && !compact ){
        
        matrix <T> C(numCol);

        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int> (size()/(n_threads));

        if(block_size < worksize){
            block_size = size();
            n_threads = 1;
        }

    #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for (size_t i = 0; i < numCol; i++) {
            for (size_t j = 0; j <= i; j++) {
                C.A[i*(i + 1)/2 + j] = A[i*numCol+j];
            }
        }

        compact = true;
        resize(numCol);
        *this = C;
    }
    else{
        failbit = true;
        if (!allocated)
            throw std::string("Matrix is empty");
        else
            throw std::string("Matrix is in a compact form");
    }        
}

//===============================================================================================================

template <typename T>
T & qgen::matrix<T>::at(size_t atRow, size_t atCol){
    /*
        Safe access to the elements of a matrix.
        Indexes atRow and atCol caunted starting from 0.

        Return value: the element at the specified position in the container.
        Example:
                 qgen::matrix <double> M; // M is empty matrix.
                 M.at(10,10) = 0.8;          // M now is (11,11) matrix.
                 double val = M.at(10,10);   // val = 0.8.
    */

    size_t index = 0;

    if( !compact ){

        if ((atRow+1 > numRow) || (atCol+1 > numCol)){
            numRow = _max(numRow, atRow+1);
            numCol = _max(numCol, atCol+1);
            if (allocated){
                resize();
            }
            else{
                int status = allocate (numRow,numCol);
                if (status != 0){
                    failbit = true;
                    throw std::string("Memory allocation error");
                }
                allocated = true;
                resizedElements = numRow*numCol;
            }
        }            
        index = atRow*numCol+atCol;

    }
    else{
        if (atRow < atCol){
            failbit = true;
            throw std::string("For a matrix in a compact form 'columns > rows' is not allowed");
        }
        if( atRow+1 > numRow ){
            numRow = numCol = atRow+1;
            if (allocated){
                resize();
            }
            else{
                int status = allocate (numRow);
                if (status != 0){
                    failbit = true;
                    throw std::string("Memory allocation error");
                }
                allocated = true;
                resizedElements = (numRow*numRow + numRow)/2;
            }
        }           
        index = atRow*(atRow+1)/2+atCol;
    }

    return A[index];
}

//===============================================================================================================

template <typename T>
size_t qgen::matrix<T>::size(){
    /*  
        Returns number of elements allocated for A.

        Return value: size_t value, number of elements in a matrix.
        Example:
            qgen::matrix <double> M(3); // symmetric matrix in a compact form.
            double val = M.size(); // val = 6.
            
    */

    if ( !compact )
        return numRow*numCol;
    else
        return (numCol*numCol+numCol)/2;
}

//===============================================================================================================

template <typename T>
qgen::matrix<T> qgen::matrix<T>::operator+(const matrix<T>& rhs) {
    /*
        Matrix addition.
        Accepts a general rectangular matrix as well as symmetric matrix in a compact form.

        Return value: matrix.
        Example:
                 qgen::matrix <double> M1(2,3);
                 qgen::matrix <double> M2(2,3);
                 qgen::matrix <double> res;
                 res = M1 + M2;
    */

    matrix <T> C;

    if ( !eq(rhs) ){
        C.failbit = true;
        throw std::string("Matrices dimensions are not equal");
    }
        
    if( !rhs.compact ){
        int status = C.allocate(rhs.numRow, rhs.numCol);
        if(status != 0) {
            C.failbit = true;
            throw std::string("Memory allocation error");
        }
        
        C.allocated = true;
        resizedElements = rhs.numRow*rhs.numCol;
        C.numCol = rhs.numCol;
        C.numRow = rhs.numRow;
    }
    else{
        int status = C.allocate(rhs.numRow);
        if(status != 0) {
            C.failbit = true;
            throw std::string("Memory allocation error");
        }
            
        C.allocated = true;
        resizedElements = (rhs.numRow*rhs.numRow+rhs.numRow)/2;
        C.numCol = rhs.numRow;
        C.numRow = rhs.numRow;
    }

    auto n_threads = std::thread::hardware_concurrency();
    auto block_size = static_cast<unsigned int> (C.size()/(n_threads));
    
    if(block_size < worksize){
        block_size = C.size();
        n_threads = 1;
    }

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
    for (size_t i = 0; i < C.size(); i++)
        C.A[i] = A[i] + rhs.A[i];

    return matrix(C);
}

//===============================================================================================================

template <typename T>
qgen::matrix<T> qgen::matrix<T>::operator-(const matrix<T>& rhs) {
    /*
        Matrix substitution.
        Accepts a general rectangular matrix as well as symmetric matrix in a compact form.

        Return value: matrix.
        Example:
                 qgen::matrix <double> M1(2,3);
                 qgen::matrix <double> M2(2,3);
                 qgen::matrix <double> res;
                 res = M1 - M2;
    */

    matrix <T> C;

    if ( !eq(rhs) ){
        C.failbit = true;
        throw std::string("Matrices dimensions are not equal");
    }
            
    if( !rhs.compact ){
        int status = C.allocate(rhs.numRow, rhs.numCol);
        if(status != 0) {
            C.failbit = true;
            throw std::string("Memory allocation error");
        }
        
        C.allocated = true;
        resizedElements = rhs.numRow*rhs.numCol;
        C.numCol = rhs.numCol;
        C.numRow = rhs.numRow;
    }
    else{
        int status = C.allocate(rhs.numRow);
        if(status != 0) {
            C.failbit = true;
            throw std::string("Memory allocation error");
        }
        
        C.allocated = true;
        resizedElements = (rhs.numRow*rhs.numRow+rhs.numRow)/2;
        C.numCol = rhs.numRow;
        C.numRow = rhs.numRow;
    }

    auto n_threads = std::thread::hardware_concurrency();
    auto block_size = static_cast<unsigned int> (C.size()/(n_threads));
    
    if(block_size < worksize){
        block_size = C.size();
        n_threads = 1;
    }

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
    for (size_t i = 0; i < C.size(); i++)
        C.A[i] = A[i] - rhs.A[i];

    return matrix(C);
}

//===============================================================================================================

template <typename T>
qgen::matrix<T> qgen::matrix<T>::operator-(const T val) {
    /* 
        Overloaded '-' operator to substract a scalar from a matrix.

        Return value: matrix object.
        Example:

            qgen::matrix <double> res;    // empty matrix
            qgen::matrix <double> M(n,m); // M is (n,m) matrix initialized by 0.
            
            res = M - 2.0; // res become (n,m) matrix where all elements are -2.0;
                           // M remains unchanged
    */

    matrix <T> C;

    if( !compact ){
        int status = C.allocate(numRow, numCol);
        if(status != 0) {
            C.failbit = true;
            throw std::string("Memory allocation error");
        }
        
        C.allocated = true;
        resizedElements = numRow*numCol;
        C.numCol = numCol;
        C.numRow = numRow;
    }
    else{
        int status = C.allocate(numRow);
        if(status != 0) {
            C.failbit = true;
            throw std::string("Memory allocation error");
        }
        
        C.allocated = true;
        resizedElements = (numRow*numRow+numRow)/2;
        C.numCol = numRow;
        C.numRow = numRow;
    }

    auto n_threads = std::thread::hardware_concurrency();
    auto block_size = static_cast<unsigned int> (C.size()/(n_threads));
    
    if(block_size < worksize){
        block_size = C.size();
        n_threads = 1;
    }

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
    for (size_t i = 0; i < C.size(); i++)
        C.A[i] = A[i] - val;

    return matrix(C);
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::gemmt_intrf (double *A, double *B, MKL_INT rowA, MKL_INT colA, MKL_INT colB){
    /*
        Interface to cblas_dgemmt routine which computes a matrix-matrix product with general matrices
        but updates only the upper or lower triangular part of the result matrix.
    */

    cblas_dgemmt (CblasRowMajor, CblasLower, CblasNoTrans, CblasTrans, rowA, colA, 1.0, A, colA, A, colA, 0.0, B, colB);

}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::gemmt_intrf (float *A, float *B, MKL_INT rowA, MKL_INT colA, MKL_INT colB){
    /*
        Interface to cblas_sgemmt routine which computes a matrix-matrix product with general matrices
        but updates only the upper or lower triangular part of the result matrix.
    */

    cblas_sgemmt (CblasRowMajor, CblasLower, CblasNoTrans, CblasTrans, rowA, colA, 1.0, A, colA, A, colA, 0.0, B, colB);

}
//===============================================================================================================

template <typename T>
qgen::matrix<T> qgen::matrix<T>::operator^(const int val) {
	/*
	    Overloaded matrix power operator which computes:

	    1)  matrix product: A * A' if val == 2.
            Return value: rectangular matrix.

        2)  inverse of matrix A: A^-1 if val == -1.
            Return value: (a) general rectangular matrix if A is not symmetric;
                          (b) matrix in compact form if A is symmetric.

        3)  inversion of matrix product: (A * A')^-1 if val == -2.
            Return value: rectangular matrix.
        
        Original matrix A remains unchanged.

        Example:
            qgen::matrix <double> M(n,n);
            qgen::matrix <double> res;
            
            for (auto i = 0; i < M.size(); i++){
                M[i] = static_cast <double> (i);
            
            res = M^(-2);
	*/

    matrix <T> C (*this);

    if(val == 2){

        matrix <T> tmpM;

        if ( C.compact )
            C.symtorec();

        int status = tmpM.allocate(numRow, numRow);
        if(status != 0) {
            C.failbit = true;
            throw std::string("Memory allocation error");
        }
        
        tmpM.allocated = true;
        tmpM.compact = false;
        tmpM.numRow = numRow;
        tmpM.numCol = numRow;
        tmpM.resizedElements = numRow*numRow;

        gemmt_intrf (C.A, tmpM.A, C.numRow, C.numCol, tmpM.numCol);

        /*
            Because C is symmetrical, it consists of values only for 'L' part, 'U' part are zeros.
            Therefore we restore and further return regular rectangular matrix.
        */
        for (size_t i = 0; i < tmpM.numRow; i++){
            for (size_t j = 0; j <= i; j++){
                tmpM(j,i) = tmpM(i,j);
            }
        }
        
        C = tmpM;

    }
    else if (val == -1){

        if( numRow == numCol ){
            
            if( !C.compact ){
                try{
                    inv_rec(C.A, numRow, numCol);
                }
                catch(std::string err){
                    C.failbit = true;
                    throw err;
                }
            }       
            else{
                try{
                    inv_sym(C.A, numCol);
                }
                catch(std::string err){
                    C.failbit = true;
                    throw err;
                }
            }
           
        }
        else{
            C.failbit = true;
            throw std::string("The matrix is not square");
        }
    }
    else if (val == -2){

        /* First do A*A' */

        matrix <T> tmpM;

        if ( C.compact )
            C.symtorec();

        int status = tmpM.allocate(numRow, numRow);
        if(status != 0) {
            C.failbit = true;
            throw std::string("Memory allocation error");
        }
        
        tmpM.allocated = true;
        tmpM.compact = false;
        tmpM.numRow = numRow;
        tmpM.numCol = numRow;
        tmpM.resizedElements = numRow*numRow;

        gemmt_intrf (C.A, tmpM.A, C.numRow, C.numCol, tmpM.numCol);

        /*
            Because C is symmetrical, it consists of values only for 'L' part, 'U' part are zeros.
            Therefore we restore and further return regular rectangular matrix.
        */
        for (size_t i = 0; i < tmpM.numRow; i++){
            for (size_t j = 0; j <= i; j++){
                tmpM(j,i) = tmpM(i,j);
            }
        }
        
        C = tmpM;
        tmpM.clear();

        /* Then do inv(A*A') */
            
        try{
            inv_rec(C.A, C.numRow, C.numCol);
        }
        catch(std::string err){
            C.failbit = true;
            throw err;
        }
    }
    else{
        C.failbit = true;
        throw std::string("Not supported operation");
    }

    return matrix(C);
}

//===============================================================================================================

template <typename T>
qgen::matrix<T> qgen::matrix<T>::operator^(const char *val) {
	/*
	    Overloaded matrix power operator which computes matrix transpose.
        Original matrix A remains unchanged.

        Return value: rectangular matrix.
        
        Example:
            qgen::matrix <double> M(n,m);
            qgen::matrix <double> res;
                       
            res = M^"T"; // res now is (m,n) matrix while M is (n,m) matrix.
	*/

    matrix <T> C (*this);
    if ( *val == 'T' || 't' ){
        C.transpose();
    }
    else {
        C.failbit = true;
        throw "Invalid operation";
    }
    return matrix(C);
}

//===============================================================================================================

template <typename T>
qgen::matrix<T> qgen::matrix<T>::operator*(const matrix<T>& rhs) {
	/*
	    Matrix dot product:
	    C = A * B;
        where A & B are rectangular matrices, C is always rectangular.

        Return value: rectangular matrix.
	*/
	qgen::matrix <T> C;

    /* Check if matrices can be multiplied. */
    if(numCol != rhs.numRow){
        C.failbit = true;
        throw std::string("Matrices are not consistent for multiplication");
    }
    
    int status = C.allocate(numRow, rhs.numCol);
    if(status != 0) {
        C.failbit = true;
        throw std::string("Memory allocation error");
    }
    
    C.allocated = true;
    C.numRow = numRow;
    C.numCol = rhs.numCol;
    C.resizedElements = C.numRow*C.numCol;

    dotprod(A, rhs.A, C.A, numRow, rhs.numRow, numCol, rhs.numCol);

    return matrix(C);
}

//===============================================================================================================

template <typename T>
bool qgen::matrix<T>::operator==(const matrix<T>& rhs) {
	/*
	    Checks matrices complete equality: A == B.

        Return value: logical TRUE if two matrices are equal and FALSE otherwise.
	*/
	
    if ( (rhs.numRow != numRow) || (rhs.numCol != numCol) )
        return false;

    if( (rhs.compact && compact) || (!rhs.compact && !compact) )
        return false;

    return std::equal(std::begin(A), std::end(A), std::begin(rhs.A));
}

//===============================================================================================================

template <typename T>
bool qgen::matrix<T>::eq(const matrix<T>& rhs) {
	/*
	    Checks matrices shapes equality: A == B.

        Return value: logical TRUE if shapes of two matrices are equal and FALSE otherwise.
                      Values of matrices elements are not compared.
	*/
	
    if ( (rhs.numRow != numRow) || (rhs.numCol != numCol) )
        return false;

    if( (rhs.compact && !compact) || (!rhs.compact && compact) )
        return false;

    return true;
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::dotprod (double *A, double *B, double *C, MKL_INT rowA, MKL_INT rowB, MKL_INT colA, MKL_INT colB) {

	/*
	 * Matrix product (rectangular matrices):
	 * C = A * B;
	 * 
	 * void cblas_dgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,...
	 * 					const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,...
	 * 					const MKL_INT k, const double alpha, const double *a,...
	 * 					const MKL_INT lda, const double *b, const MKL_INT ldb,...
	 * 					const double beta, double *c, const MKL_INT ldc);
	 */

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, colB, colA, 1.0, A, colA, B, colB, 0.0, C, colB);

}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::dotprod (float *A, float *B, float *C, MKL_INT rowA, MKL_INT rowB, MKL_INT colA, MKL_INT colB) {

	/*
	 * Matrix product (rectangular matrices):
	 * C = A * B;
	 * 
	 * void cblas_sgemm (const CBLAS_LAYOUT Layout, const CBLAS_TRANSPOSE transa,...
	 * 					const CBLAS_TRANSPOSE transb, const MKL_INT m, const MKL_INT n,...
	 * 					const MKL_INT k, const float alpha, const float *a,...
	 * 					const MKL_INT lda, const float *b, const MKL_INT ldb,...
	 * 					const float beta, float *c, const MKL_INT ldc);
	 */

    cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rowA, colB, colA, 1.0, A, colA, B, colB, 0.0, C, colB);

}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::clear() {
    /*
        Clears a memory allocated for the container A.

        Return value: none.
    */

    if(allocated) {
        mkl_free(A);
        allocated = false;
        failbit = false;
        failinfo = 0;
        resizedElements = 0;
        numRow = numCol = 0;
    }
}

//===============================================================================================================

template <typename T>
bool qgen::matrix<T>::empty() {
    /*
        Checks if a matrix is empty.

        Return value: logical TRUE if matrix is empty, otherwise returns FALSE.
    */
    
    if(allocated)
        return false;

    return true;
}

//===============================================================================================================

template <typename T>
qgen::matrix<T> & qgen::matrix<T>::operator=(const matrix<T>& rhs){
    /*
        Overloaded assignment operator.
    */
        
    qgen::matrix <T> tmpObj (rhs);

    std::swap(compact, tmpObj.compact);
    std::swap(rectangular, tmpObj.rectangular);
    std::swap(symetric, tmpObj.symetric);        
    std::swap(failbit, tmpObj.failbit);
    std::swap(failinfo, tmpObj.failinfo);        
    std::swap(numCol, tmpObj.numCol);
    std::swap(numRow, tmpObj.numRow);
    std::swap(resizedElements, tmpObj.resizedElements);
    std::swap(binFilename, tmpObj.binFilename);
    std::swap(allocated, tmpObj.allocated);
    std::swap(A, tmpObj.A);

    return *this;
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::fwrite(){
    /*
        Moves matrix to a DISK and clears memory allocated for the container A.

        Return value: none.
    */

    fA.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
    fA.open( binFilename, fA.binary | fA.trunc | fA.out );

    if (!fA.is_open()) {
        failbit = true;
        throw std::string("Error while opening a binary file");
    }
    
    size_t sz;
    if(!compact)
        sz = numRow*numCol;
    else
        sz = (numCol*numCol+numCol)/2;
    
    fA.write(reinterpret_cast<char*>(A), sz*sizeof( T ));

    fA.close();

    if (allocated){
        mkl_free(A);
        allocated = false;
    }
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::fread(){
    /*
        Moves saved on DISK matrix back to the memory.

        Return value: none.
    */

    fA.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
    fA.open( binFilename, fA.binary | fA.in );

    if (!fA.is_open()) {
        failbit = true;
        throw std::string("Error while opening a binary file");
    }
    
    size_t sz;
    if(!compact){

        sz = numRow*numCol;

        int status = allocate(numRow,numCol);
        if(status != 0) throw std::string("Memory allocation error");
        
        allocated = true;
    }
    else{

        sz = (numCol*numCol+numCol)/2;

        int status = allocate(numCol);
        if(status != 0) throw std::string("Memory allocation error");
        
        allocated = true;
    }
    
    fA.read(reinterpret_cast<char*>(A), sz*sizeof( T ));

    fA.close();
}

//===============================================================================================================

template <typename T>
void qgen::matrix <T>::inv_rec (double *A, MKL_INT rowA, MKL_INT colA){
    /*
        Matrix inversion. Interface to mkl routines.
        Here:
            A - matrix to invert;
            row - number of rows in A;
            col - number of columns in A;
            ipiv - empty vector of size lda.
        
        1. First, make LU factorizatioon of A: A_fact = P*L*U; ipiv = P.
        2. Second, make inversion of factorized A: A = inv(A_fact).
        The result of inversion is in matrix A.

        Return value: none.
    */

    lapack_int info = 0;
    lapack_int row = rowA;
    lapack_int col = colA;
    int matrix_order = LAPACK_ROW_MAJOR;

    lapack_int *ipiv;
    ipiv = (lapack_int *)mkl_malloc( row*sizeof( lapack_int ), sizeof(T)*8 );
    if (ipiv == NULL) {
        mkl_free(ipiv);
        failbit = true;
        throw std::string("Memory allocation error");
    }
    for (lapack_int i = 0; i < (row); i++)
        ipiv[i] = 1;


    info = LAPACKE_dgetrf (matrix_order, row, col, A, col, ipiv);
    if (info != 0){
        mkl_free(ipiv);
        throw std::string("Error during computation of the LU factorization of a general m-by-n matrix");
    }

    info = LAPACKE_dgetri (matrix_order, row, A, row, ipiv);
    if (info != 0){
        mkl_free(ipiv);
        throw std::string("Error during computation the inverse of an LU-factored general matrix");
    }
    
    mkl_free(ipiv);
}

//===============================================================================================================

template <typename T>
void qgen::matrix <T>::inv_rec (float *A, MKL_INT rowA, MKL_INT colA){
    /*
        Matrix inversion. Interface to mkl routines.
        Here:
            A - matrix to invert;
            row - number of rows in A;
            col - number of columns in A;
            ipiv - empty vector of size lda.
        
        1. First, make LU factorizatioon of A: A_fact = P*L*U; ipiv = P.
        2. Second, make inversion of factorized A: A = inv(A_fact).
        The result of inversion is in matrix A.

        Return value: none.
    */

    lapack_int info = 0;
    lapack_int row = rowA;
    lapack_int col = colA;
    int matrix_order = LAPACK_ROW_MAJOR;

    lapack_int *ipiv;
    ipiv = (lapack_int *)mkl_malloc( row*sizeof( lapack_int ), sizeof(T)*8 );
    if (ipiv == NULL) {
        mkl_free(ipiv);
        failbit = true;
        throw std::string("Memory allocation error");
    }
    for (lapack_int i = 0; i < (row); i++)
        ipiv[i] = 1;


    info = LAPACKE_sgetrf (matrix_order, row, col, A, col, ipiv);
    if (info != 0){
        mkl_free(ipiv);
        throw std::string("Error during computation of the LU factorization of a general m-by-n matrix");
    }

    info = LAPACKE_sgetri (matrix_order, row, A, row, ipiv);
    if (info != 0){
        mkl_free(ipiv);
        throw std::string("Error during computation the inverse of an LU-factored general matrix");
    }
    
    mkl_free(ipiv);
}

//===============================================================================================================

template <typename T>
void qgen::matrix <T>::inv_sym (double *A, MKL_INT colA){
    /*
        Symetric matrix inversion. Interface to mkl routines.

        Return value: none.
    */
    
    lapack_int info = 0;

    int matrix_order = LAPACK_ROW_MAJOR;

    info = LAPACKE_dpptrf (matrix_order, 'L', colA, A);
    if (info != 0) {
        failbit = true;
        failinfo = info;
        throw std::string("Error during computationof  the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage");
    }

    info = LAPACKE_dpptri (matrix_order, 'L', colA, A);
    if (info != 0) {
        failbit = true;
        failinfo = info;
        throw std::string("Error during computation the inverse of a packed symmetric (Hermitian) positive-definite matrix after the Cholesky factorization");
    }
}

//===============================================================================================================

template <typename T>
void qgen::matrix <T>::inv_sym (float *A, MKL_INT colA){
    /*
        Symetric matrix inversion. Interface to mkl routines.

        Return value: none.
    */
    
    lapack_int info = 0;

    int matrix_order = LAPACK_ROW_MAJOR;

    info = LAPACKE_spptrf (matrix_order, 'L', colA, A);
    if (info != 0) {
        throw std::string("Error during computationof  the Cholesky factorization of a symmetric (Hermitian) positive-definite matrix using packed storage");
    }

    info = LAPACKE_spptri (matrix_order, 'L', colA, A);
    if (info != 0) {
        throw std::string("Error during computation the inverse of a packed symmetric (Hermitian) positive-definite matrix after the Cholesky factorization");
    }
}

//===============================================================================================================

template <typename T>
void qgen::matrix <T>::invert(){
    /*
        Symetric/General matrix inversion.

        Return value: none.

        Example:

            qgen::matrix <double> M(n,n);
            qgen::matrix <double> res;    // empty matrix
            
            for (auto i = 0; i < M.size(); i++)
                M[i] = i;
            
            M.invert(); // now M is (n,n) inverted matrix
            res = M;       // now res is (n,n) inverted matrix

    */


    if( numRow == numCol ){
        if( !compact ){
           try{
               inv_rec(A, numRow, numCol);
           }
           catch(std::string err){
               failbit = true;
               throw err;
           }
        }       
        else{
           try{
               inv_sym(A, numCol);
           }
           catch(std::string err){
               failbit = true;
               throw err;
           }
        }
    }
    else{
        failbit = true;
        throw std::string("Matrix is not square");
    }
}

//===============================================================================================================

template <typename T>
size_t qgen::matrix<T>::capacity(){
    /*  
        Returns number of allocated elements for A.
    */

    return resizedElements;
}

//===============================================================================================================

template <typename T>
qgen::matrix <T>::~matrix(){
    /*
        Class destructor.
    */

    remove( binFilename.c_str() );
    if( allocated ) {
        mkl_free(A);
        allocated = false;
    }

}

//===============================================================================================================

template <typename T>
qgen::matrix<T>::matrix( const matrix<T>& obj){
    /*
        Copy constructor.
    */

    compact = obj.compact;
    rectangular = obj.rectangular;
    symetric = obj.symetric;   
    failbit = obj.failbit;
    failinfo = obj.failinfo;   
    numCol = obj.numCol;
    numRow = obj.numRow;
    resizedElements = obj.resizedElements;
    binFilename = obj.binFilename;  
    allocated = false;
    resize(numRow, numCol);
    std::memcpy( A, obj.A, sizeof(T)*this->size() );
}

//===============================================================================================================

template <typename T>
void qgen::matrix<T>::print(std::string whiichMatrix) {
    /*
        Prints part of a matrix into a LOG file.

        Return value: none.
    */

	FILE *dbgFile;
	dbgFile = fopen(debug_file.c_str(), "a");		
	
    int maxRows = 15;
    fprintf (dbgFile, "%s", whiichMatrix.c_str());
    //fprintf (dbgFile, "\n");
    if (rectangular){
        fprintf (dbgFile, "%s", ", Rectangular matrix");
        fprintf (dbgFile, "\n\n");
        for (auto i=0; i<_min(maxRows,numRow); i++) {
            for (auto j=0; j<_min(maxRows,numCol); j++) {
                fprintf (dbgFile, "%12.5G", A[i*numCol+j]);
            }
            fprintf (dbgFile, "\n");
        }
    }
    else if(symetric){
        fprintf (dbgFile, "%s", ", symetric matrix");
        fprintf (dbgFile, "\n\n");
        for (auto i=0; i<_min(maxRows,numRow); i++) {
            for (auto j=0; j<=i; j++) {
                fprintf (dbgFile, "%12.5G", A[i*(i+1)/2+j]);
            }
            fprintf (dbgFile, "\n");
        }
    }
    fprintf (dbgFile, "\n");
    fprintf (dbgFile, "\n");

    fclose(dbgFile);
}

//===============================================================================================================
