/*
 * cs_housekeeper.hpp
 *
 * Main class APY definition.
 *
 *  Created on: Aug 5, 2017
 *      Author: vimi
 */

#pragma once

#include <map>
#include <vector>
#include <sstream>
#include <chrono>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <thread>

#include "structures.hpp"


class HOUSEKEEPER {

public:

	// constructor
	HOUSEKEEPER(/*std::string ped_file, std::string typed_file, std::string gmat_file, std::string log_file, int whichMethod, bool inbreeding, double weight, std::string adjust*/);

	// destructor
	~HOUSEKEEPER();

	template <class X, class Y>
	unsigned long mapCapacity(const std::map<X,Y> &themap){
		unsigned long cap = sizeof(themap);
		unsigned long noElements = 0;
		for ( typename std::map<X,Y>::const_iterator i = themap.begin() ; i != themap.end() ; i++ )
			noElements++;
		cap += noElements * sizeof(X);
		cap += noElements * sizeof(Y);
		return cap;
	}

	template <class T>
	unsigned long arrCapacity(T *&a, int size);

	/* methods to get initial data from files */
	int read_ped ();
	int read_gtypd ();

	template <typename T, typename X>
	int read_genmT (std::vector<X> &g_row, std::vector<X> &g_col, std::vector<T> &g_val);

	/* methods to check if IDs are in pedigree and in G matrix */
	int check_gmatID ();
	int check_gtypedID ();

	void write_log (std::string location, int code);
	void runtimelog (bool wrTime, std::string message, double someValue, bool makeFileOpen, bool makeFileClose, bool showSize = false);
	void writeOut (bool wrTime, std::string message, double someValue, bool makeFileOpen, bool makeFileClose);

	int trace_back (std::map <dKey, pedPair> &from_ped, std::map <dKey, pedPair> &to_ped, std::vector <int> &gen_id);
	int trace_back_par (std::map <dKey, pedPair> &from_ped, std::map <dKey, pedPair> &to_ped, std::vector <int> &gen_id);

	int get_ainv (std::map <dKey, pedPair> &ped, std::map <dKey, double> &ai, bool inbreed);
	int get_A22 (std::map <dKey, pedPair> &ped, float *A22);
	int get_dinv (std::map <dKey, pedPair> &ped, std::vector <double> &dinv, bool inbreed);
	int findRecodedIdMap (std::map <int, int> &id_map, std::vector <int> &whereIidList, std::vector <int> &whatIdList);
	int getRecodedIdMap (std::map <int, int> &id_map, std::vector <int> &idVect);
	int getA22vector (std::vector <float> &w, std::vector <int> &v, std::vector<std::vector<int> > &Ped);

	std::string getFilePath (std::string &fileName);

	/* class data */
	std::map <dKey, pedPair> f_ped; /* original pedigree from red file. will be deleted after use */
	std::map <int, int> birth_id_map; /* the map which holds <animal_id, birth_day>, used when check correctness of pedigree */
	std::map <dKey, pedPair> pedigree; /* traced back original pedigree */
	std::map <dKey, pedPair> r_pedigree; /* reduced pedigree container, sorted by day order during initialization */
	std::map <dKey, double> ainv; /* A_inverse matrix container */
	std::map <dKey, double> r_ainv;
	std::map <dKey, float> a; /* A matrix container */
	std::vector <int> genotypedID;      /* container for genotyped IDs (from 'typed' file) */
	std::vector <int> coreID;           /* container for core IDs (from 'typed' file) */
	std::vector <int> gmatID;           /* container for the list of G matrix IDs */
	std::vector <int> pedID;            /* container for the list of pedigree IDs */
	std::vector<double> inbrF; /* inbreeding coeffecients */


	bool inbreeding;
	double wa;
	std::string adjust;

	std::map <int, int> rcode_map;
	std::map <int, int> fcode_map;

	/* variables to track the progress of G matrix file reading */
	int theLock;
	int modifFlag;
	bool tmpSatatus;
	size_t gFileSize;
	size_t logFilePos_bgn;
	size_t logFilePos_end;

protected:

	/* input data (files) should be used by class constructor */
	std::string gmat_file;
	std::string ped_file;
	std::string typed_file;
	std::string marker_file;
	std::string log_file;
	std::string runtimeLog_file = "RUNTIME_HMAT.LOG";
	std::fstream runtimelogStr;
	int whichMethod;

	void check (std::vector<int> &id, std::vector<int> &check_id, std::vector<int> &w_id);
	void write_log2 (std::string which_id, int code, std::vector<int> &miss_id); /* used only for missing IDs*/
	bool logfile_exists (const std::string& name); // LINUX version

	bool find_invect (std::vector <int> &where, int what);

	template <typename T>
	int find_invect2 (std::vector <T> &where, int what);

	bool find_inmap (std::map <int, int>  &where, int what);
	int pos_inped (std::map <int, int> &codemap, int id);
	bool is_unique(std::vector<int> &x);

};

//===============================================================================================================

template <class T>
unsigned long arrCapacity(T *&a, int size){
	return size*size*sizeof(T);
}

//===============================================================================================================

/*template <class X, class Y>
unsigned long mapCapacity(const std::map<X,Y> &themap){
	unsigned long cap = sizeof(themap);
	unsigned long noElements = 0;
	for ( typename std::map<X,Y>::const_iterator i = themap.begin() ; i != themap.end() ; i++ )
		noElements++;
	cap += noElements * sizeof(X);
	cap += noElements * sizeof(Y);
	return cap;
}*/

//===============================================================================================================

template <typename T>
int HOUSEKEEPER::find_invect2 (std::vector <T> &where, int what) {

	typename std::vector<T>::iterator it;
	it = find (where.begin(), where.end(), what);
	if (it != where.end()) {
		return it - where.begin();
	}
	else
		return -1;

}

//===============================================================================================================

template <typename T, typename X>
int HOUSEKEEPER::read_genmT (std::vector<X> &g_row, std::vector<X> &g_col, std::vector<T> &g_val) {

	/* reads G matrix from gmat.dat file into the linked list container */

	std::string where("Error during reading G matrix");
	std::string line;
	T val;
	dKey key;

	char* end; // new
	const char* p; // new
	std::vector<T> tmp_list;

	std::ifstream ped;
	ped.exceptions ( std::ifstream::failbit | std::ifstream::badbit );

	size_t diagonals = 0; /* for debugging */
	size_t fSize;
	theLock = 1;
	size_t curSize = 101;

	try {
		//id = 0;
		ped.open(getFilePath(gmat_file), std::fstream::in);

		ped.seekg (0, ped.end);
		fSize = ped.tellg();
		ped.seekg (0, ped.beg);

		while (getline (ped,line)) {
			p = line.c_str();
			for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end)) {
				p = end;
				if (errno == ERANGE){
					write_log ("Range error during reading G matrix", 10);
					errno = 0;
				}
				tmp_list.push_back(static_cast<T> (f));
			}
			g_row.push_back(static_cast<X> (tmp_list[0]));
			g_col.push_back(static_cast<X> (tmp_list[1]));
			g_val.push_back(tmp_list[2]);
			gmatID.push_back(int(tmp_list[0]));
			gmatID.push_back(int(tmp_list[1]));

			if (static_cast<X> (tmp_list[0]) == static_cast<X> (tmp_list[1])) diagonals++;

			tmp_list.erase(tmp_list.begin(), tmp_list.end());

			/* Change the size of read data, this will be used to track the read status operation */
			if ( (100*ped.tellg()/fSize)%1 == 0 && theLock && (100*ped.tellg()/fSize) != curSize) {
				theLock = 0;
				gFileSize = 100*ped.tellg()/fSize;
				curSize = gFileSize;
				modifFlag = 1;
				theLock = 1;
			}
		}
		ped.close();

		if(!is_unique(gmatID)) gmatID.erase( unique( gmatID.begin(), gmatID.end() ), gmatID.end() ); // here the vector should be sorted and unique

		return 0;
	}
	catch (std::exception const& e) {
		if (ped.eof()) {
			ped.close();

			//std::cout<<"diagonals = "<<diagonals<<std::endl;

			if(!is_unique(gmatID)) gmatID.erase( unique( gmatID.begin(), gmatID.end() ), gmatID.end() ); // here the vector should be sorted and unique

			if (diagonals != gmatID.size()) {
				write_log ("There are missing diagonals in G-matrix file.", 2);
				return 2;
			}
			else
				return 0;
		}
		else {
			write_log (e.what(), 10);
			//std::cerr << "Exception opening/reading file: "<< e.what() << std::endl;
			return 10;
		}
	}
	catch (...)
	{
		write_log (where, 1);
		return 1;
	}

	return 0;
}

//===============================================================================================================
