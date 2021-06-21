/*
*
*
*
*/

#ifndef CS_GMAT_HPP
#define CS_GMAT_HPP

#pragma once

#include "structures.hpp"
#include "../../matrixclass/mathlib/src/cs_matrix.hpp"
#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <thread>
#include <cmath>

#ifdef _WIN64
	#include <Windows.h>
#endif

#ifdef linux
	#include <sys/stat.h>
	#include <linux/limits.h>
	#include <unistd.h>
#endif

#define workload 100000

class GMAT
{
    public:

        GMAT(const std::string& dataFile);
        ~GMAT();

        int makeG();
        void getValues(std::vector <int>&, std::vector <int>&, std::vector <float>&, double);
        void getIDs(std::vector <int>&);

        qgen::matrix <double> G;         /* G matrix */

        /*double ofd_mean = 0.0;
		double all_mean = 0.0;
		double diag_mean = 0.0;*/

    private:

        std::string marker_file;
        std::map<size_t, int> animID;    /* key: consecutive index, value: animal ID */
        std::map <int, std::string> SNP; /* initial markers data */
        qgen::matrix <double> Z;         /* Z (snp) matrix */
        size_t snpNum;
        double freq;

        int parseString(std::string& snpStr, std::vector<int>& markers);
        int getZ();
        int getG();
        int readSNP();
        //void getMeans();

        /* here are some helper methiods; because we do not want to have inheritance from HOUSEKEEPER class */
        void write_log( std::string location, int code );
        bool logfile_exists( std::string& name );
        std::string getFilePath( std::string &fileName );
};

#endif
