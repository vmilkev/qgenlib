/*
 * cs_apy.hpp
 *
 *  Created on: Aug 7, 2017
 *      Author: vimi
 */
#pragma once

#include "structures.hpp"
#include "cs_housekeeper.hpp"
#include "cs_matrixmath.hpp"

#include "cs_gmat.hpp"

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <sstream>

#include <algorithm>
#include <map>
#include <string>

#include <stdio.h>
#include <stdlib.h>
#include <cstdlib>

#include <omp.h>

#ifdef _WIN64
	#include <Windows.h>
#endif

#ifdef linux
	#include <sys/stat.h>
	#include <linux/limits.h>
	#include <unistd.h>
#endif


class APY : public HOUSEKEEPER {

public:

	APY(std::string f1, std::string f2, std::string f3, std::string f4, int m, bool in, double wt, std::string adj, std::string mark, double d_val);
	~APY();

	int run_apy (); /*  the main method called from the outside */

	std::map <dKey, double> r_ainv; /* A matrix container */

	/* Opptional;  the values red from file 'param' if it exists in working directory */

	bool isH; /* optional parameter: if true - produce full H_inv in addition to H_sub */
	bool isHtxt; /* optional parameter: if true - produce also ASCII file of H_sub.dat in addition to H_sub.bin file */
	double tol_a22inv; /* optional parameter: defines the threshold to zero out the elements in a22inv matrix */
	double tol_Hsub; /* optional parameter: defines the threshold to zero out the elements in H_sub matrix */
	unsigned int threadsAvail; /* the number of available threads set in OMP_NUM_THREADS */
	unsigned int n_threads;
	bool save; /* if TRUE - use compact storage form inversion; if FALSE - use full storage inversion */
	double gdiag_val;

private:

	int check_ID (); /* check if animal IDs are in G matrix and in pedigree */
	int write_toFile_h (std::string where, double *what, std::vector<int> id_list, bool isBinary);
	int write_toFile2 (std::string where, std::map <dKey, double> &what);

};
