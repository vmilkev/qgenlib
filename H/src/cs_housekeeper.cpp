/*
 * cs_housekeeper.cpp
 *
 * Description of methods of main class APY.
 *
 *  Created on: Aug 5, 2017
 *      Author: vimi
 */


#include <fstream>
#include <stdexcept>
#include <sstream>
#include <math.h>

#ifdef _WIN32
#include <Windows.h>
#elif defined __linux__
#include <sys/stat.h>
#include <linux/limits.h>
#include <unistd.h> // LINUX API
#endif

#include <stdlib.h>
#include <string.h>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <array>
#include "cs_housekeeper.hpp"

//===============================================================================================================

HOUSEKEEPER::HOUSEKEEPER/*(std::string f1, std::string f2, std::string f3, std::string f4, int meth, bool inb, double w, std::string adj):\
ped_file(f1), typed_file(f2), gmat_file(f3), log_file(f4), whichMethod(meth), inbreeding(inb), wa(w), adjust(adj)*/() {

}

//===============================================================================================================

HOUSEKEEPER::~HOUSEKEEPER(){

}

//===============================================================================================================

bool HOUSEKEEPER::logfile_exists (const std::string& name) {
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

//===============================================================================================================

void HOUSEKEEPER::write_log (std::string location, int code) {

	std::fstream log;
	bool exists = logfile_exists (log_file);
	if (exists)
		log.open (log_file, std::fstream::out | std::fstream::app);
	else log.open (log_file, std::fstream::out);

	auto end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	if (log) log << std::ctime(&end_time) << " :  Exception -> "<<location<<";  code -> "<<code<<std::endl;
	log.close();
}

//===============================================================================================================

void HOUSEKEEPER::runtimelog (bool wrTime, std::string message, double someValue, bool makeFileOpen, bool makeFileClose, bool showSize) {

	char *ctime_no_newline;

	if (makeFileOpen) {
		runtimelogStr.open (runtimeLog_file, std::fstream::out);

		auto timeNow = std::chrono::system_clock::now();
		auto outTime = std::chrono::system_clock::to_time_t(timeNow);

		if (runtimelogStr) {
			runtimelogStr<<"HMAT PROGRAM RUNTIME LOG FILE   "<<std::ctime(&outTime)<<std::endl;

			runtimelogStr<<"INPUT DATA:"<<std::endl;
			runtimelogStr<<"          pedigree is in the file:                "<<ped_file<<std::endl;
			runtimelogStr<<"          genotyped IDs are in the file:          "<<typed_file<<std::endl;
			if ( !marker_file.empty() )
				runtimelogStr<<"          genomic markers data is in the file: "<<marker_file<<std::endl;
			else
				runtimelogStr<<"          genomic covariance data is in the file: "<<gmat_file<<std::endl;
			runtimelogStr<<"          apy method:                             "<<whichMethod<<std::endl;
			runtimelogStr<<"          inbreeding:                             "<<inbreeding<<std::endl;
			runtimelogStr<<"          weight:                                 "<<wa<<std::endl;
			runtimelogStr<<"          adjust:                                 "<<adjust<<std::endl;
		}
	}

	auto timeNow = std::chrono::system_clock::now();
	auto outTime = std::chrono::system_clock::to_time_t(timeNow);

	ctime_no_newline = strtok(std::ctime(&outTime), "\n");

	if (wrTime && !showSize) {
		if (runtimelogStr) {
			runtimelogStr <<std::endl;
			runtimelogStr << ctime_no_newline << "   "<<message<<std::endl;
			logFilePos_end = runtimelogStr.tellg();
		}
	}
	else if (runtimelogStr && !showSize)  {
		runtimelogStr << "                           "<<message<<"  "<<someValue<<std::endl;
		logFilePos_end = runtimelogStr.tellg();
	}
	else if (showSize) {
		runtimelogStr.seekg (logFilePos_end, runtimelogStr.beg);
		runtimelogStr << ctime_no_newline << "             "<<message<<"  "<<someValue<<std::endl;
	}

	if (makeFileClose) runtimelogStr.close();

}
//===============================================================================================================

void HOUSEKEEPER::writeOut (bool wrTime, std::string message, double someValue, bool makeFileOpen, bool makeFileClose) {

	if (theLock && modifFlag) {
		theLock = 0;
		runtimelog (false, "Completed, % ", gFileSize, false, false, 1);
		modifFlag = 0;
		theLock = 1;
	}
}

//===============================================================================================================

void HOUSEKEEPER::write_log2 (std::string which_id, int code, std::vector<int> &miss_id) {

	std::fstream log;
	bool exists = logfile_exists (log_file);
	if (exists)
		log.open (log_file, std::fstream::out | std::fstream::app);
	else log.open (log_file, std::fstream::out);
	if (log) {

		auto end = std::chrono::system_clock::now();
		std::time_t end_time = std::chrono::system_clock::to_time_t(end);

		log << std::ctime(&end_time)<< " :  Exception code -> "<<code<<std::endl;
		auto count = 1;
		for (auto const& elem : miss_id) {
			log <<"			"<<count<<". "<<which_id<<", missing id: "<<elem<<std::endl;
			count++;
		}
	}
	log.close();
}

//===============================================================================================================

int HOUSEKEEPER::read_ped () {

	std::string where("Error during reading pedigree file");
	std::string line;
	//int id, sir, dam, dayy;
	dKey key;
	pedPair p_pair;

	char* end; // new
	const char* p; // new
	std::vector<double> tmp_list;
	//tmp_list.resize(4);

	std::ifstream ped;
	ped.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
	
	/*
	std::cout<<"ped file1: "<<ped_file<<std::endl;
	std::cout<<"ped file2: "<<getFilePath(ped_file)<<", short name: "<<ped_file<<std::endl;
	*/

	try {
		ped.open(getFilePath(ped_file), std::fstream::in);
		if (ped) {
			while (getline (ped,line)) {
				p = line.c_str();
				for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end)) {
					p = end;
					if (errno == ERANGE){
						write_log ("Range error during reading pedigree file", 10);
						errno = 0;
					}
					tmp_list.push_back(f);
				}
				key.day = static_cast<int>(tmp_list[3]);
				key.id = static_cast<int>(tmp_list[0]);
				p_pair.id_1 = static_cast<int>(tmp_list[1]);
				p_pair.id_2 = static_cast<int>(tmp_list[2]);
				f_ped[key] = p_pair;
				birth_id_map[static_cast<int>(tmp_list[0])] = static_cast<int>(tmp_list[3]);
				pedID.push_back(static_cast<int>(tmp_list[0]));

				tmp_list.erase(tmp_list.begin(), tmp_list.end());

			}
			ped.close();
			if(!is_unique(pedID)) pedID.erase( unique( pedID.begin(), pedID.end() ), pedID.end() ); // here the vector should be sorted and unique
			return 0;

		}
		else throw 11;
	}
	catch (std::exception const& e) {
		if (ped.eof()) {
			ped.close();
			if(!is_unique(pedID)) pedID.erase( unique( pedID.begin(), pedID.end() ), pedID.end() ); // here the vector should be sorted and unique
			return 0;
		}
		else {
			write_log (e.what(), 10);
			return 10;
		}
	}
	catch (int e_int) {
		write_log (where, e_int); // error on opening file
		return e_int;
	}
	catch (...)
	{
		write_log (where, 1);
		return 1;
	}

	return 0;
}

//===============================================================================================================

int HOUSEKEEPER::read_gtypd () {

	/* reads genotyped and core IDs from typed file into the vectors: genotyped, core */

	std::string where("Error during reading genotyped IDs file");
	std::string line;
	int t_gtyp, t_gcor;
	t_gtyp = t_gcor = 0;

	char* end; // new
	const char* p; // new
	std::vector<double> tmp_list;

	std::ifstream ped;
	ped.exceptions ( std::ifstream::failbit | std::ifstream::badbit );

	try {
		ped.open(getFilePath(typed_file), std::fstream::in);
		while (getline (ped,line)) {
			p = line.c_str();
			for (double f = std::strtod(p, &end); p != end; f = std::strtod(p, &end)) {
				p = end;
				if (errno == ERANGE){
					write_log ("Range error during reading genotyped IDs file", 10);
					errno = 0;
				}
				tmp_list.push_back(f);
			}
			tmp_list.erase(tmp_list.begin(), tmp_list.end());

			if (static_cast<int>(tmp_list[0]) != 0 && static_cast<int>(tmp_list[0]) != t_gtyp) genotypedID.push_back(static_cast<int>(tmp_list[0]));
			if (static_cast<int>(tmp_list[1]) != 0 && static_cast<int>(tmp_list[0]) != t_gtyp) coreID.push_back(static_cast<int>(tmp_list[0]));
			t_gtyp = static_cast<int>(tmp_list[0]);
			t_gcor = static_cast<int>(tmp_list[1]);

		}
		ped.close();
		if(!is_unique(genotypedID)) genotypedID.erase( unique( genotypedID.begin(), genotypedID.end() ), genotypedID.end() ); // here the vector should be sorted as well
		if(!is_unique(coreID)) coreID.erase( unique( coreID.begin(), coreID.end() ), coreID.end() ); // here the vector should be sorted as well
		return 0;
	}
	catch (std::exception const& e) {
		if (ped.eof()) {
			if(!is_unique(genotypedID)) genotypedID.erase( unique( genotypedID.begin(), genotypedID.end() ), genotypedID.end() ); // here the vector should be sorted as well
			if(!is_unique(coreID)) coreID.erase( unique( coreID.begin(), coreID.end() ), coreID.end() ); // here the vector should be sorted as well
			return 0;
		}
		else {
			write_log (e.what(), 10);
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

void HOUSEKEEPER::check (std::vector<int> &id, std::vector<int> &check_id, std::vector<int> &w_id) {
	// id - vector where to check
	// check_id - vector of what to check
	// w_id - vector of missing values

	for (size_t i = 0; i < check_id.size(); i++)
		if( !std::binary_search(id.begin(), id.end(), check_id[i]) )
			w_id.push_back(check_id[i]);

}

//===============================================================================================================

int HOUSEKEEPER::check_gmatID () {

	// are all GENOTYPES (from G matrix) in pedigree ?

	std::string where("check_gmatID");
	int status = 0;
	try
	{
		std::vector<int> w_ID;
		check(pedID, gmatID, w_ID);
		if (!w_ID.empty()) {
			std::string where("G matrix IDs is not in pedigree");
			write_log2 (where, 20, w_ID);
			status = 20;
			w_ID.clear();
		}
	}
	catch (std::exception const& e)
	{
		write_log (e.what(), 1);
		write_log (where, 1);
		status = 1;
		return status;
	}

	return status;
}

//===============================================================================================================

int HOUSEKEEPER::check_gtypedID () {

	// are all GENOTYPED IDs (from typed file) in pedigree and in G matrix ?

	std::string where("check_gtypedID");
	int status = 0;
	try
	{
		std::vector<int> w_ID;
		check(gmatID, genotypedID, w_ID);
		if (!w_ID.empty()) {
			std::string where("typed IDs is not in G matrix");
			write_log2 (where, 20, w_ID);
			status = 20;
			w_ID.clear();
		}
	}
	catch (std::exception const& e)
	{
		write_log (e.what(), 1);
		write_log (where, 1);
		status = 1;
		return status;
	}
	return status;
}

//===============================================================================================================

bool HOUSEKEEPER::find_invect (std::vector <int> &where, int what) {

	std::vector<int>::iterator it;
	it = find (where.begin(), where.end(), what);
	if (it != where.end())
		return true;
	else
		return false;

}

//===============================================================================================================

int HOUSEKEEPER::findRecodedIdMap (std::map <int, int> &id_map, std::vector <int> &whereIidList, std::vector <int> &whatIdList) {

	try
	{
		for (size_t i = 0; i < whatIdList.size(); i++) {
			id_map[whatIdList[i]] = find_invect2(whereIidList, whatIdList[i])+1;
		}
	}
	catch (std::exception const& e) {
		write_log (e.what(), 1);
		write_log ("Exception in findRecodedIdMap", 1);
		return 1;
	}

	return 0;
}

//===============================================================================================================

int HOUSEKEEPER::getRecodedIdMap (std::map <int, int> &id_map, std::vector <int> &idVect) {

	try
	{
		size_t code_id = 1;
		for (auto const& elem : idVect) {
			id_map[elem] = code_id;
			code_id++;
		}
	}
	catch (std::exception const& e) {
		write_log (e.what(), 1);
		write_log ("Exception in getRecodedIdMap", 1);
		return 1;
	}

	return 0;
}

//===============================================================================================================

int HOUSEKEEPER::trace_back (std::map <dKey, pedPair> &from_ped, std::map <dKey, pedPair> &to_ped, std::vector <int> &gen_id) {

	try
	{
		pedPair p_pair;
		dKey key;
		std::vector <int> gen_pedID(gen_id);
		int elem_v, iter_v = 0;
		bool exists;

		std::vector <int> days;
		std::vector <int> ids;
		std::vector <int> sire;
		std::vector <int> dame;

		for (auto const& elem_m : from_ped) {
			days.push_back(elem_m.first.day);
			ids.push_back(elem_m.first.id);
			sire.push_back(elem_m.second.id_1);
			dame.push_back(elem_m.second.id_2);
		}

		do
		{
			elem_v = gen_pedID[iter_v];

			exists = false;

			for (size_t i = 0; i < days.size(); i++/*auto const& elem_m : from_ped*/) {

				if (ids[i]/*elem_m.first.id*/ == elem_v) {

					exists = true;
					p_pair.id_1 = sire[i]/*elem_m.second.id_1*/;
					p_pair.id_2 = dame[i]/*elem_m.second.id_2*/;

					if ( (birth_id_map[ sire[i] ] >= days[i]) || (birth_id_map[ dame[i] ] >= days[i]) ) throw 30;

					key.day = days[i];
					key.id = elem_v;
					to_ped[key/*elem_m.first*/] = p_pair;

					if (p_pair.id_1 != 0) {
						if (!find_invect (gen_pedID, p_pair.id_1))
							gen_pedID.push_back(p_pair.id_1);
					}
					if (p_pair.id_2 != 0) {
						if (!find_invect (gen_pedID, p_pair.id_2))
							gen_pedID.push_back(p_pair.id_2);
					}

					break;
				}

			}

			if (!exists) {
				key.day = -1;
				key.id = elem_v;
				p_pair.id_1 = 0;
				p_pair.id_2 = 0;
				to_ped[key] = p_pair;
				exists = false;
			}

			iter_v++;

		} while (iter_v < gen_pedID.size());

	}
	catch (int ex) {
		write_log("Tracing Pedigree", ex);
		return ex;
	}
	catch (std::exception const& e)
	{
		write_log (e.what(), 1);
		write_log ("Exception in trace_back", 1);
		return 1;
	}

	return 0;
}

//===============================================================================================================

int HOUSEKEEPER::trace_back_par (std::map <dKey, pedPair> &from_ped, std::map <dKey, pedPair> &to_ped, std::vector <int> &gen_id) {

	auto sub_size = ceil(gen_id.size()/20);

	std::vector<int> pedID1;
	std::vector<int> pedID2;
	std::vector<int> pedID3;
	std::vector<int> pedID4;
	std::vector<int> pedID5;
	std::vector<int> pedID6;
	std::vector<int> pedID7;
	std::vector<int> pedID8;
	std::vector<int> pedID9;
	std::vector<int> pedID10;

	std::vector<int> pedID11;
	std::vector<int> pedID12;
	std::vector<int> pedID13;
	std::vector<int> pedID14;
	std::vector<int> pedID15;
	std::vector<int> pedID16;
	std::vector<int> pedID17;
	std::vector<int> pedID18;
	std::vector<int> pedID19;
	std::vector<int> pedID20;

	std::map<dKey, pedPair> pedigree1;
	std::map<dKey, pedPair> pedigree2;
	std::map<dKey, pedPair> pedigree3;
	std::map<dKey, pedPair> pedigree4;
	std::map<dKey, pedPair> pedigree5;
	std::map<dKey, pedPair> pedigree6;
	std::map<dKey, pedPair> pedigree7;
	std::map<dKey, pedPair> pedigree8;
	std::map<dKey, pedPair> pedigree9;
	std::map<dKey, pedPair> pedigree10;

	std::map<dKey, pedPair> pedigree11;
	std::map<dKey, pedPair> pedigree12;
	std::map<dKey, pedPair> pedigree13;
	std::map<dKey, pedPair> pedigree14;
	std::map<dKey, pedPair> pedigree15;
	std::map<dKey, pedPair> pedigree16;
	std::map<dKey, pedPair> pedigree17;
	std::map<dKey, pedPair> pedigree18;
	std::map<dKey, pedPair> pedigree19;
	std::map<dKey, pedPair> pedigree20;

	int status, status1, status2, status3, status4, status5, status6, status7, status8, status9, status10;

	status = 0;

#pragma omp parallel sections //num_threads (10)
	{
#pragma omp section
		{
			for (auto i = 0; i < 1*sub_size; i++) {
				pedID1.push_back(gen_id[i]);
			}
			status1 = 0;
			status1 = trace_back (from_ped, pedigree1, pedID1);
		}
#pragma omp section
		{
			for (auto i = 1*sub_size; i < 2*sub_size; i++) {
				pedID2.push_back(gen_id[i]);
			}
			status2 = 0;
			status2 = trace_back (from_ped, pedigree2, pedID2);
		}
#pragma omp section
		{
			for (auto i = 2*sub_size; i < 3*sub_size; i++) {
				pedID3.push_back(gen_id[i]);
			}
			status3 = 0;
			status3 = trace_back (from_ped, pedigree3, pedID3);
		}
#pragma omp section
		{
			for (auto i = 3*sub_size; i < 4*sub_size; i++) {
				pedID4.push_back(gen_id[i]);
			}
			status4 = 0;
			status4 = trace_back (from_ped, pedigree4, pedID4);
		}
#pragma omp section
		{
			for (auto i = 4*sub_size; i < 5*sub_size; i++) {
				pedID5.push_back(gen_id[i]);
			}
			status5 = 0;
			status5 = trace_back (from_ped, pedigree5, pedID5);
		}
#pragma omp section
		{
			for (auto i = 5*sub_size; i < 6*sub_size; i++) {
				pedID6.push_back(gen_id[i]);
			}
			status6 = 0;
			status6 = trace_back (from_ped, pedigree6, pedID6);
		}
#pragma omp section
		{
			for (auto i = 6*sub_size; i < 7*sub_size; i++) {
				pedID7.push_back(gen_id[i]);
			}
			status7 = 0;
			status7 = trace_back (from_ped, pedigree7, pedID7);
		}
#pragma omp section
		{
			for (auto i = 7*sub_size; i < 8*sub_size; i++) {
				pedID8.push_back(gen_id[i]);
			}
			status8 = 0;
			status8 = trace_back (from_ped, pedigree8, pedID8);
		}
#pragma omp section
		{
			for (auto i = 8*sub_size; i < 9*sub_size; i++) {
				pedID9.push_back(gen_id[i]);
			}
			status9 = 0;
			status9 = trace_back (from_ped, pedigree9, pedID9);
		}
#pragma omp section
		{
			for (auto i = 9; i < 10*sub_size; i++) {
				pedID10.push_back(gen_id[i]);
			}
			status10 = 0;
			status10 = trace_back (from_ped, pedigree10, pedID10);
		}
	}

	status = status1 + status2 + status3 + status4 + status5 + status6 + status7 + status8 + status9 + status10;

	if (status) {
		write_log ("Exception during the call of trace_back_par", 1);
		return status;
	}

	to_ped.insert(pedigree1.begin(), pedigree1.end()); pedigree1.clear();
	to_ped.insert(pedigree2.begin(), pedigree2.end()); pedigree2.clear();
	to_ped.insert(pedigree3.begin(), pedigree3.end()); pedigree3.clear();
	to_ped.insert(pedigree4.begin(), pedigree4.end()); pedigree4.clear();
	to_ped.insert(pedigree5.begin(), pedigree5.end()); pedigree5.clear();
	to_ped.insert(pedigree6.begin(), pedigree6.end()); pedigree6.clear();
	to_ped.insert(pedigree7.begin(), pedigree7.end()); pedigree7.clear();
	to_ped.insert(pedigree8.begin(), pedigree8.end()); pedigree8.clear();
	to_ped.insert(pedigree9.begin(), pedigree9.end()); pedigree9.clear();
	to_ped.insert(pedigree10.begin(), pedigree10.end()); pedigree10.clear();


#pragma omp parallel sections //num_threads (10)
	{
#pragma omp section
		{
			for (auto i = 10*sub_size; i < 11*sub_size; i++) {
				pedID11.push_back(gen_id[i]);
			}
			status1 = 0;
			status1 = trace_back (from_ped, pedigree11, pedID11);
		}
#pragma omp section
		{
			for (auto i = 11*sub_size; i < 12*sub_size; i++) {
				pedID12.push_back(gen_id[i]);
			}
			status2 = 0;
			status2 = trace_back (from_ped, pedigree12, pedID12);
		}
#pragma omp section
		{
			for (auto i = 12*sub_size; i < 13*sub_size; i++) {
				pedID13.push_back(gen_id[i]);
			}
			status3 = 0;
			status3 = trace_back (from_ped, pedigree13, pedID13);
		}
#pragma omp section
		{
			for (auto i = 13*sub_size; i < 14*sub_size; i++) {
				pedID14.push_back(gen_id[i]);
			}
			status4 = 0;
			status4 = trace_back (from_ped, pedigree14, pedID14);
		}
#pragma omp section
		{
			for (auto i = 14*sub_size; i < 15*sub_size; i++) {
				pedID15.push_back(gen_id[i]);
			}
			status5 = 0;
			status5 = trace_back (from_ped, pedigree15, pedID15);
		}
#pragma omp section
		{
			for (auto i = 15*sub_size; i < 16*sub_size; i++) {
				pedID16.push_back(gen_id[i]);
			}
			status6 = 0;
			status6 = trace_back (from_ped, pedigree16, pedID16);
		}
#pragma omp section
		{
			for (auto i = 16*sub_size; i < 17*sub_size; i++) {
				pedID17.push_back(gen_id[i]);
			}
			status7 = 0;
			status7 = trace_back (from_ped, pedigree17, pedID17);
		}
#pragma omp section
		{
			for (auto i = 17*sub_size; i < 18*sub_size; i++) {
				pedID18.push_back(gen_id[i]);
			}
			status8 = 0;
			status8 = trace_back (from_ped, pedigree18, pedID18);
		}
#pragma omp section
		{
			for (auto i = 18*sub_size; i < 19*sub_size; i++) {
				pedID19.push_back(gen_id[i]);
			}
			status9 = 0;
			status9 = trace_back (from_ped, pedigree19, pedID19);
		}
#pragma omp section
		{
			for (auto i = 19*sub_size; i < gen_id.size(); i++) {
				pedID20.push_back(gen_id[i]);
			}
			status10 = 0;
			status10 = trace_back (from_ped, pedigree20, pedID20);
		}
	}

	status = status1 + status2 + status3 + status4 + status5 + status6 + status7 + status8 + status9 + status10;

	if (status) {
		write_log ("Exception during the call of trace_back_par", 1);
		return status;
	}

	to_ped.insert(pedigree11.begin(), pedigree11.end()); pedigree11.clear();
	to_ped.insert(pedigree12.begin(), pedigree12.end()); pedigree12.clear();
	to_ped.insert(pedigree13.begin(), pedigree13.end()); pedigree13.clear();
	to_ped.insert(pedigree14.begin(), pedigree14.end()); pedigree14.clear();
	to_ped.insert(pedigree15.begin(), pedigree15.end()); pedigree15.clear();
	to_ped.insert(pedigree16.begin(), pedigree16.end()); pedigree16.clear();
	to_ped.insert(pedigree17.begin(), pedigree17.end()); pedigree17.clear();
	to_ped.insert(pedigree18.begin(), pedigree18.end()); pedigree18.clear();
	to_ped.insert(pedigree19.begin(), pedigree19.end()); pedigree19.clear();
	to_ped.insert(pedigree20.begin(), pedigree20.end()); pedigree20.clear();

	return 0;
}
//===============================================================================================================

int HOUSEKEEPER::get_ainv (std::map <dKey, pedPair> &ped, std::map <dKey, double> &ai, bool inbreed) {

	int status = 0;
	dKey akey;
	std::vector <double> di; // store D(-1)

	try
	{
		status = get_dinv (ped, di, inbreed);
		if (status != 0) throw status;

		int s, d, id;
		double dinv;
		for (auto const& elem : ped) {

			dinv = di.front();
			s = elem.second.id_1;
			d = elem.second.id_2;
			id = elem.first.id;

			if (s && d) {

				akey.day = id;
				akey.id = id;
				ai[akey] = ai[akey] + dinv;

				akey.day = id;
				akey.id = s;
				if (s > id) {
					akey.day = s;
					akey.id = id;
				}
				ai[akey] = ai[akey] - dinv/2.0;

				akey.day = id;
				akey.id = d;
				if (d > id) {
					akey.day = d;
					akey.id = id;
				}
				ai[akey] = ai[akey] - dinv/2.0;

				akey.day = s;
				akey.id = s;
				ai[akey] = ai[akey] + dinv/4.0;

				akey.day = d;
				akey.id = d;
				ai[akey] = ai[akey] + dinv/4.0;

				if (s >= d) {
					akey.day = s;
					akey.id = d;
				}
				else {
					akey.day = d;
					akey.id = s;
				}
				ai[akey] = ai[akey] + dinv/4.0;

			}
			else if (!s && !d) {
				akey.day = id;
				akey.id = id;
				ai[akey] = ai[akey] + dinv;
			}
			else {
				akey.day = id;
				akey.id = id;
				ai[akey] = ai[akey] + dinv;

				akey.day = id;
				if (s) {
					akey.id = s;
					if (s > id) {
						akey.day = s;
						akey.id = id;
					}
				}
				else {
					akey.id = d;
					if (d > id) {
						akey.day = d;
						akey.id = id;
					}
				}
				ai[akey] = ai[akey] - dinv/2.0;

				if (s) {
					akey.day = s;
					akey.id = s;
				}
				else {
					akey.day = d;
					akey.id = d;
				}
				ai[akey] = ai[akey] + dinv/4.0;
			}
			di.erase(di.begin());

		}

	}
	catch (std::exception const& e)
	{
		status = 1;
		write_log (e.what(), status);
		write_log ("Exception during the call of get_dinv", status);
		return status;
	}
	catch (int ex)
	{
		write_log ("Exception in get_ainv", ex);
		return ex;
	}

	return 0;
}

//===============================================================================================================

int HOUSEKEEPER::get_A22 (std::map <dKey, pedPair> &ped, float *A22) {

	try
	{
		/*

		// allocate space in a disc for full 'A' and initialize it by '0.0'
		// we will build full 'A' directly on hard drive due to its big size
		std::string filename1 = "A.bin";
		std::fstream fA(filename1, fA.binary | fA.trunc | fA.in | fA.out);
		fA.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fA.is_open()) throw 10;

		size_t sz = (ped.size()*ped.size() + ped.size())/2;

		try {
			float *a; a = (float *)malloc( 1*sizeof( float ));
			for (auto i = 0; i < sz; i++) {
				a[0] = 0.0f;
				fA.write(reinterpret_cast<char*> (a), 1*sizeof( float ));
			}

			free(a);
		}
		catch (std::exception const& e) {
			std::cerr << "Exception writing file: "<< e.what() << std::endl;
			return 10;
		}



		float inbr = 0.0f;
		if (inbreed) inbr = 1.0f;

		std::unordered_map <int, float>::iterator it, it2;
		std::unordered_map<int, float> a_map;
		std::cout<<"map size : "<<a_map.size()<<std::endl;


		//float *a; a = (float *)malloc( 1*sizeof( float )); // temporal array to read/write in binary file

		for (auto e_i = ped.begin(); e_i != std::next(ped.end(),2); e_i++) {

			int row = code_map[e_i->first.id];

			for (auto e_j = ped.begin(); e_j != std::next(e_i,1); e_j++) {

				int col = code_map[e_j->first.id];
				int i = row*(row-1)/2 + col-1;
				int si = e_i->second.id_1;
				int di = e_i->second.id_2;

				int s = code_map[e_i->second.id_1];
				int d = code_map[e_i->second.id_2];

				if (row == col) {
					int ind;
					if (si && di) {
						ind = s*(s-1)/2 + d-1;
						if (d > s) ind = d*(d-1)/2 + s-1;

						//do operation: Agen[row*(row-1)/2 + col-1] = 1.0 + 0.5 * Agen[ind];
						it = a_map.find(ind);
						if (it == a_map.end()) a_map[i] = 1.0f;
						else a_map[i] = 1.0 + 0.5 * a_map[ind];


						fA.seekg(sizeof(float)*ind, std::ios::beg);
						fA.read(reinterpret_cast<char*> (a), 1*sizeof( float ));
						a[0] = 1.0 + 0.5 * a[0];
						fA.seekp(sizeof(float)*i, std::ios::beg);
						fA.write(reinterpret_cast<char*> (a), 1*sizeof( float ));

					}
					else {
						//do operation: Agen[row*(row-1)/2 + col-1] = 1.0;
						a_map[i] = 1.0;


						a[0] = 1.0;
						fA.seekp(sizeof(float)*i, std::ios::beg);
						fA.write(reinterpret_cast<char*> (a), 1*sizeof( float ));

					}
				}
				else {
					int ind_a, ind_b;
					if (si && di) {
						ind_a = s*(s-1)/2 + col-1;
						if (col > s) ind_a = col*(col-1)/2 + s-1;
						ind_b = d*(d-1)/2 + col-1;
						if (col > d) ind_b = col*(col-1)/2 + d-1;

						//do operation: Agen[row*(row-1)/2 + col-1] = 0.5 * (Agen[ind_a] + Agen[ind_b]);
						it = a_map.find(ind_a);
						it2 = a_map.find(ind_b);
						if (it != a_map.end() && it2 != a_map.end()) a_map[i] = 0.5 * (a_map[ind_a] + a_map[ind_b]);
						else if (it != a_map.end() && it2 == a_map.end()) a_map[i] = 0.5 * (a_map[ind_a]);
						else if (it == a_map.end() && it2 != a_map.end()) a_map[i] = 0.5 * (a_map[ind_b]);


						fA.seekg(sizeof(float)*ind_a, std::ios::beg);
						fA.read(reinterpret_cast<char*> (a), 1*sizeof( float ));
						auto r_a = a[0];
						fA.seekg(sizeof(float)*ind_b, std::ios::beg);
						fA.read(reinterpret_cast<char*> (a), 1*sizeof( float ));
						auto r_b = a[0];
						a[0] = 0.5 * (r_a + r_b);
						fA.seekp(sizeof(float)*i, std::ios::beg);
						fA.write(reinterpret_cast<char*> (a), 1*sizeof( float ));

					}
					else if (!si && !di) {
						//do operation: Agen[row*(row-1)/2 + col-1] = 0.0;


						a[0] = 0.0;
						fA.seekp(sizeof(float)*i, std::ios::beg);
						fA.write(reinterpret_cast<char*> (a), 1*sizeof( float ));

					}
					else {
						if (si && !di) {
							ind_a = s*(s-1)/2 + col-1;
							if (col > s) ind_a = col*(col-1)/2 + s-1;
							//do operation: Agen[row*(row-1)/2 + col-1] = 0.5 * Agen[ind_a];
							it = a_map.find(ind_a);
							if (it == a_map.end()) a_map[i] = 0.0f;
							else a_map[i] = 0.5 * a_map[ind_a];


							fA.seekg(sizeof(float)*ind_a, std::ios::beg);
							fA.read(reinterpret_cast<char*> (a), 1*sizeof( float ));
							a[0] = 0.5 * a[0];
							fA.seekp(sizeof(float)*i, std::ios::beg);
							fA.write(reinterpret_cast<char*> (a), 1*sizeof( float ));

						}
						else {
							ind_b = d*(d-1)/2 + col-1;
							if (col > d) ind_b = col*(col-1)/2 + d-1;
							//do operation: Agen[row*(row-1)/2 + col-1] = 0.5 * Agen[ind_b];
							it = a_map.find(ind_b);
							if (it == a_map.end()) a_map[i] = 0.0f;
							else a_map[i] = 0.5 * a_map[ind_b];


							fA.seekg(sizeof(float)*ind_b, std::ios::beg);
							fA.read(reinterpret_cast<char*> (a), 1*sizeof( float ));
							a[0] = 0.5 * a[0];
							fA.seekp(sizeof(float)*i, std::ios::beg);
							fA.write(reinterpret_cast<char*> (a), 1*sizeof( float ));

						}
					}
				}
			}
		}

		std::cout<<"map size : "<<a_map.size()<<std::endl;

		// fill in 'Agen' from 'A':
		std::vector<int> posinGen;
		std::vector<int> posinA;


		for (auto e_i = ped.begin(); e_i != std::next(ped.end(),2); e_i++) {

			int pos = find_invect2(genotypedID, e_i->first.id);
			if (pos != -1) {
				posinA.push_back(code_map[e_i->first.id]);
				posinGen.push_back(pos+1);
			}
		}

		for (auto i = 0; i < posinA.size(); i++) {
			for (auto j = 0; j <= i; j++) {
				int ind_m, ind_f;
				ind_m = posinGen[i]*(posinGen[i]-1)/2 + posinGen[j] - 1;
				if (posinGen[i] < posinGen[j]) ind_m = posinGen[j]*(posinGen[j]-1)/2 + posinGen[i] - 1;

				ind_f = (posinA[i]*(posinA[i]-1)/2 + posinA[j] - 1);
				if (posinA[i] < posinA[j]) ind_f = (posinA[j]*(posinA[j]-1)/2 + posinA[i] - 1);
				Agen[ind_m] = a_map[ind_f];


				ind_f = sizeof(float)*(posinA[i]*(posinA[i]-1)/2 + posinA[j] - 1);
				if (posinA[i] < posinA[j]) ind_f = sizeof(float)*(posinA[j]*(posinA[j]-1)/2 + posinA[i] - 1);
				fA.seekg(ind_f, std::ios::beg);
				fA.read(reinterpret_cast<char*> (a), 1*sizeof( float ));
				Agen[ind_m] = a[0];

			}
		}

		//free(a);

		 */
		//--------------------------------------------------------------------------------------------------

		int n = ped.size(); // total amount of animals in pedigree
		std::vector<std::vector<int> > Ped(n+1, std::vector<int>(2, 0.0));
		std::vector<int> GenID; // list of genotyped IDs
		//std::vector<float> w(n+1, 0.0);
		//std::vector<int> v(n+1, 0);

		std::map <int, int> code_map;
		std::map <int, int> gen_map;

		int code_id = 1;
		int t_k = 0;
		for (auto const& elem : ped) {

			code_map[elem.first.id] = code_id;

			Ped[code_id][0] = pos_inped (code_map, elem.second.id_1);
			Ped[code_id][1] = pos_inped (code_map, elem.second.id_2);

			/*
			if (find_invect(genotypedID, elem.first.id)) {
				GenID.push_back(code_id);
			}
			 */

			int pos = find_invect2(genotypedID, elem.first.id);
			if (pos != -1) {
				GenID.push_back(code_id);
				gen_map[code_id] = pos;

				/*
				std::cout<<"i = "<<t_k<<"; code_id = "<<code_id<<"; pos = "<<pos<<std::endl;
				t_k++;
				 */
			}

			code_id++;
		}

		int m = GenID.size(); // number of genotyped IDs

#pragma omp parallel for
		for (size_t i = 0; i < m; i++) {

			std::vector<float> w(n+1, 0.0f);
			std::vector<int> v(n+1, 0);

			size_t ii = gen_map[GenID[i]];

			v[GenID[i]] = 1;

			auto  status = getA22vector (w, v, Ped);
			if (status != 0) throw status;

			for (size_t j = 0; j < m; j++) {
				//A22[ii*(ii+1)/2+j] = w[ GenID[ j ] ];
				if (ii >= gen_map[GenID[j]]) {
					/*
					size_t ind = ii*(ii+1)/2+gen_map[GenID[j]];
					if (ind >= (genotypedID.size()*genotypedID.size()+genotypedID.size())/2)
						std::cout<<"ii, jj: "<<ii<<"; "<< gen_map[GenID[j]]<<"; ind = "<<ind<< std::endl;
					 */

					A22[ii*(ii+1)/2+gen_map[GenID[j]]] = w[ GenID[j] ];
				}
			}

			//v[GenID[i]] = 0;
		}

	}
	catch (std::exception const& e)
	{
		write_log (e.what(), 1);
		write_log ("Exception in get_A22", 1);
		return 1;
	}
	catch (int ex)
	{
		write_log ("Exception during the call of getA22vector", ex);
		return ex;
	}

	return 0;

}

//===============================================================================================================

int HOUSEKEEPER::getA22vector (std::vector <float> &w, std::vector <int> &v, std::vector<std::vector<int> > &Ped) {

	try
	{
		int n = w.size()-1;
		std::vector<float> q(n+1, 0.0);

		for (auto i = n; i >= 1; i--) {
			q[i] += v[i];
			auto s = Ped[i][0];
			auto d = Ped[i][1];
			if (s) q[s] += q[i]*0.5;
			if (d) q[d] += q[i]*0.5;
		}

		for (auto i = 1; i <= n; i++) {
			auto s = Ped[i][0];
			auto d = Ped[i][1];
			auto di = (std::count (Ped[i].begin(), Ped[i].end(), 0)+2.0)/4.0 - 0.25 * (inbrF[s] + inbrF[d]);
			float temp = 0.0;
			if (s) temp += w[s];
			if (d) temp += w[d];
			w[i] = 0.5*temp;
			w[i] += di*q[i];
		}
	}
	catch (std::exception const& e)
	{
		write_log (e.what(), 1);
		write_log ("Exception in getA22vector", 1);
		return 1;
	}
	return 0;
}

//===============================================================================================================

int HOUSEKEEPER::pos_inped (std::map <int, int> &codemap, int id) {

	auto pos = codemap.find(id);
	if (pos == codemap.end()) {
		return 0;
	} else {
		return pos->second;
	}
}

//===============================================================================================================

int HOUSEKEEPER::get_dinv (std::map <dKey, pedPair> &ped, std::vector <double> &dinv, bool inbreed) {

	try
	{

		if (inbreed)
		{
			/*
			// This algorithm is not efficient for LARGE data since it require full memory allocation for 2D vector 'l'

			int n = ped.size();
			int pd[n+1][2];

			std::map <int, int> code_map;
			int code_id = 1;

			for (auto const& elem : ped) {

				code_map[elem.first.id] = code_id; // local recoding (get consecutive IDs)

				pd[code_id][0] = pos_inped (code_map, elem.second.id_1);
				pd[code_id][1] = pos_inped (code_map, elem.second.id_2);

				std::cout<<"elem.first.id = "<< elem.first.id<<"; code_id = " << code_id <<std::endl;
				std::cout<< "sire & its pos: " << elem.second.id_1 << "; " << pos_inped (code_map, elem.second.id_1) << std::endl;
				std::cout<< "dame & its pos: " << elem.second.id_2 << "; " << pos_inped (code_map, elem.second.id_2) << std::endl;

				code_id++;
			}

			std::vector<int> as;
			std::vector<int> ad;
			std::vector<double> dii(n+1, 1.0f);

			//start = std::chrono::high_resolution_clock::now();

			std::vector<std::vector<double>> l(n+1, std::vector<double>(n+1,0.0f));

			std::vector<double> F(n+1);

			F[0] = -1.0f;
			for (auto i = 1; i <= n; i++) {

				F[i] = 0;
				auto si = pd[i][0];
				auto di = pd[i][1];

				if (si) {
					as.push_back(si);
					l[si][si] = 1.0;
				}
				if (di) {
					ad.push_back(di);
					l[di][di] = 1.0;
				}

				while ( !as.empty() && !ad.empty() ) {

					auto j_iter = std::max_element(std::begin(as), std::end(as));
					auto j = *j_iter;
					auto k_iter = std::max_element(std::begin(ad), std::end(ad));
					auto k = *k_iter;

					auto sj = pd[j][0];
					auto dj = pd[j][1];
					auto sk = pd[k][0];
					auto dk = pd[k][1];

					if ( j > k) {

						if (sj) {
							as.push_back(sj);
							l[si][sj] = l[si][sj] + 0.5 * l[si][j];
						}
						if (dj) {
							as.push_back(dj);
							l[si][dj] = l[si][dj] + 0.5 * l[si][j];
						}
						//as.delete(j)
						as.erase(std::remove(as.begin(), as.end(), j), as.end());
					}
					else if ( k > j) {

						if (sk) {
							ad.push_back(sk);
							l[di][sk] = l[di][sk] + 0.5 * l[di][k];
						}
						if (dk) {
							ad.push_back(dk);
							l[di][dk] = l[di][dk] + 0.5 * l[di][k];
						}
						//ad.delete(k)
						ad.erase(std::remove(ad.begin(), ad.end(), k), ad.end());
					}
					else {
						if (sj) {
							as.push_back(sj);
							ad.push_back(sj);
							l[si][sj] = l[si][sj] + 0.5 * l[si][j];
							l[di][sj] = l[di][sj] + 0.5 * l[di][j];
						}
						if (dj) {
							as.push_back(dj);
							ad.push_back(dj);
							l[si][dj] = l[si][dj] + 0.5 * l[si][j];
							l[di][dj] = l[di][dj] + 0.5 * l[di][j];
						}

						//std::cout<<"j = " <<j<<std::endl;
						F[i] = F[i] + l[si][j] * l[di][j] * 0.5 * dii[j];
						//inbrF[i] = F[i];

						as.erase(std::remove(as.begin(), as.end(), j), as.end());
						ad.erase(std::remove(ad.begin(), ad.end(), j), ad.end());
					}
											
				}

				if (si && di) {
					dii[i] = 0.5 - 0.25 * (F[si] + F[di]);
					double g = 1 / dii[i];
					dinv.push_back(g);
				}
				else if (!si && !di) {
					dii[i] = 1.0;
					dinv.push_back(dii[i]);
				}
				else {
					if (si) dii[i] = 0.75 - 0.25 * (F[si]);
					else dii[i] = 0.75 - 0.25 * (F[di]);
					double g = 1 / dii[i];
					dinv.push_back(g);
				}

			}
			*/

			// Another algorithm: M.Sargolzaei et al. 'A fast algorithm for computing inbreeding ...' (works correct!)

			/**/
			int n = ped.size();
			int m = n + 1;
			int Ped[n+1][2];

			std::map <int, int> code_map;
			int code_id = 1;
			inbrF.push_back(0.0); // we shall be able to get the coeff. starting from index '1'
			for (auto const& elem : ped) {

				code_map[elem.first.id] = code_id;

				Ped[code_id][0] = pos_inped (code_map, elem.second.id_1);
				Ped[code_id][1] = pos_inped (code_map, elem.second.id_2);

				inbrF.push_back(0.0); // initialize vector of inbreeding coeff.

/* 				std::cout<<"elem.first.id = "<< elem.first.id<<"; code_id = " << code_id <<std::endl;
				std::cout<< "sire & its pos: " << elem.second.id_1 << "; " << pos_inped (code_map, elem.second.id_1) << std::endl;
				std::cout<< "dame & its pos: " << elem.second.id_2 << "; " << pos_inped (code_map, elem.second.id_2) << std::endl;
 */
				code_id++;
			}

			//n = number of animals in total; m = number of sires and dams in total

			//Variables declaration
			int i, j, k, rN, rS, S, D, MIP;

			//int Ped[n + 1][2]; // main pedigree
			int rPed[m + 1][2]; // reduced pedigree
			int SId[n + 1]; // will contain the sorted animals ID based on the ID of their sires
			int Link[n + 1]; // will contain new ID of ancestors at position of their original ID
			//int MaxIdP[m + 1]; // will contain maximum new ID of parents for each paternal group at position of the new ID of each sire
			std::vector <int> MaxIdP(m + 1, 0.0);
			double F[n + 1]; // inbreeding coefficients
			double B[m + 1]; // within family segregation variances
			double x[m + 1]; // x arrays

			//Kernel of the modified algorithm

			// set values for the unknown parent
			F[0] = -1.0;
			x[0] = 0.0;
			Link[0] = 0;
			rN = 1;
			i = 1;

			for(; i <= n; ++i){ // extract and recode ancestors

				// initialization
				SId[i] = i;
				Link[i] = 0;
				if(i <= m) x[i] = 0.0;

				S = Ped[i][0];
				D = Ped[i][1];

				if(S && !Link[S]){
					MaxIdP[rN] = Link[S] = rN;
					rPed[rN][0] = Link[Ped[S][0]];
					rPed[rN++][1] = Link[Ped[S][1]];
				}
				if(D && !Link[D]){
					/*MaxIdP[rN] =*/ Link[D] = rN;
					rPed[rN][0] = Link[Ped[D][0]];
					rPed[rN++][1] = Link[Ped[D][1]];
				}
				if(MaxIdP[Link[S]] < Link[D]) MaxIdP[Link[S]] = Link[D]; // determine maximum ID of parents for each paternal group
			}

			/*for (auto i = 1; i <= m; i++)
				std::cout<<"MaxIdP[i] = "<<MaxIdP[i]<<std::endl;*/

			// sort animals according to ID of their sires into SId
			// The pedigree is already sorted so the oldest animal IDs appears first
			for (auto i = 1; i <= n; i++) {
				SId[i] = i;
			}

			k = 1;
			i = 1;
			for(; i <= n;){ // do from the oldest sire to the youngest sire

				int t = Ped[SId[i]][0];
				if(!t) {
					F[SId[i++]] = 0.0; // sire is not known
				}
				else{

					S = Ped[SId[i]][0];
					rS = Link[S];
					MIP = MaxIdP[rS];
					x[rS] = 1.0;

					//std::cout<<"i = "<<i<<"; S = "<<S<<"; rS = "<<rS<<"; MIP = "<< MIP <<std::endl;

					for(; k <= S; ++k) {// compute within family segregation variances
						if(Link[k]) B[Link[k]] = 0.5 - 0.25 * (F[Ped[k][0]] + F[Ped[k][1]]);
					}

					for(j = rS; j; --j){ // trace back the reduced pedigree
						if(x[j]){ // consider only ancestors of the sire

							if(rPed[j][0]) x[rPed[j][0]] += x[j] * 0.5;
							if(rPed[j][1]) x[rPed[j][1]] += x[j] * 0.5;
							x[j] *= B[j]; //std::cout<<"j, x[j]: "<<j<<", "<<x[j]<<std::endl;
						}
					}

					for(j = 1; j <= MIP; ++j) {// trace forth the reduced pedigree
						//std::cout<<"j, x[j]: "<<j<<", "<<x[j]<<std::endl;
						x[j] += (x[rPed[j][0]] + x[rPed[j][1]]) * 0.5;
						//std::cout<<"j, x[j]: "<<j<<", "<<x[j]<<std::endl;
					}

					for(; i <= n; ++i) // obtain F for progeny of the current sire
						if(S != Ped[SId[i]][0]) break;
						else F[SId[i]] = x[Link[Ped[SId[i]][1]]] * 0.5;

					for(j = 1; j <= MIP; ++j) {
						x[j] = 0.0; // set to 0 for next evaluation of sire
					}
				}
			}

			for (auto i = 1; i <= n; i++) {

				auto s = Ped[i][0];
				auto d = Ped[i][1];

				inbrF[i] = F[i]; // fill global vector of inbr. coeff.

				if (s && d) {
					double g = 1.0 / (0.5 - 0.25 * (F[s] + F[d]));
					dinv.push_back(g);
				}
				else if (!s && !d) {
					double g = 1.0;
					dinv.push_back(g);
				}
				else {
					double g;
					if (s) g = 1.0/(0.75 - 0.25 * (F[s]));
					else g = 1.0/(0.75 - 0.25 * (F[d]));
					dinv.push_back(g);
				}

			}
			

		}
		else {
			for (auto i = 0; i <= ped.size(); i++)
				inbrF.push_back(0.0);

			int s, d;
			for (auto const& elem : ped) {
				s = elem.second.id_1;
				d = elem.second.id_2;
				if (s && d) dinv.push_back(2.0);
				else if (!s && !d) dinv.push_back(1.0);
				else dinv.push_back(4.0/3.0);
			}
		}
		

	}
	catch (std::exception const& e)
	{
		write_log (e.what(), 1);
		write_log ("Exception in get_dinv", 1);
		return 1;
	}

	return 0;

}

//===============================================================================================================

bool HOUSEKEEPER::is_unique(std::vector<int> &x) {
	sort( x.begin(), x.end() );
	return adjacent_find( x.begin(), x.end() ) == x.end();
}

//===============================================================================================================

std::string HOUSEKEEPER::getFilePath (std::string &fileName) {

	/*
	std::array<char, 128> buffer;
	std::string result;
	std::string cmd("realpath ");
	cmd.append(fileName);
	std::cout<<"before pipe"<<std::endl;
	std::shared_ptr<FILE> pipe(popen(cmd.c_str(), "r"), pclose);
	std::cout<<"after pipe"<<std::endl;
	if (!pipe) throw std::runtime_error("popen() failed!");
	while (!feof(pipe.get())) {
		if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
			result += buffer.data();
	}
	 */
	/*
	std::string result;
	char buffer[128];
	std::string cmd("realpath ");
	cmd.append(fileName);
	std::cout<<"command: "<<cmd.c_str()<<std::endl;
	FILE* pipe = popen(cmd.c_str(), "r");
	if (!pipe) {
		std::cout<<"ERROR in pipeopen"<<std::endl;
		throw std::runtime_error("popen() failed!");
	}
	try {
		while (!feof(pipe)) {
			if (fgets(buffer, 128, pipe) != NULL)
				result += buffer;
		}
	} catch (...) {
		pclose(pipe);
		throw;
	}
	pclose(pipe);

	result.erase(remove_if(result.begin(), result.end(), isspace), result.end());
	return result;
	 */

	//char buf[PATH_MAX + 1]; /* not sure about the "+ 1" */
	char *res = realpath(fileName.c_str(), NULL);
	if (!res) {
		perror("realpath");
		exit(EXIT_FAILURE);
	}
	return fileName = res;
}

//===============================================================================================================
