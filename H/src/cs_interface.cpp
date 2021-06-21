/*
*
*
*
*/

#include "cs_interface.hpp"

#ifdef _WIN64
	#include <Windows.h>
#endif

#ifdef linux
	#include <sys/stat.h>
	#include <linux/limits.h>
	#include <unistd.h> // LINUX API
#endif

//===============================================================================================================

HMINTERFACE::HMINTERFACE (){
}

//===============================================================================================================

void HMINTERFACE::get_ParamFileName (const std::string& data_file){
	apyParFile = data_file;
}

//===============================================================================================================

int HMINTERFACE::getHmatr(){

	APY a(apyDatStr.file_ped, apyDatStr.file_typd, apyDatStr.file_gen, apyDatStr.file_log, apyDatStr.method, apyDatStr.inbreed, apyDatStr.weight, apyDatStr.adjust, apyDatStr.file_marker, apyDatStr.gdiag);

	if( a.run_apy() ){
		a.write_log ("ERROR while running APY", -1);
		return 1;
	}

	return 0;
}

//===============================================================================================================

int HMINTERFACE::getGmatr(){

	std::string dataFile = apyDatStr.file_marker;

	if( !dataFile.empty() ){

		GMAT g(dataFile);
		g.makeG();
		g.G.print("G matrix with spaces");

		//std::cout<<"Means: all, diag, ofd = "<<g.all_mean<<", "<<g.diag_mean<<", "<<g.ofd_mean<<std::endl;
		
	}

	return 0;
}

//===============================================================================================================

int HMINTERFACE::read_apy_input_data () { // !!! Redundant method (can be deleted?)

	std::string line;
	std::ifstream ped;

    std::vector<std::string> apy_data_list;

	ped.exceptions ( std::ifstream::failbit | std::ifstream::badbit );

	/* some default values */
	apyDatStr.file_log = "LOG_HMAT.log";
	apyDatStr.adjust = " ";
	apyDatStr.method = 0;
	apyDatStr.inbreed = 0;
	apyDatStr.weight = 0.0f;

	try {

		char *res_file = realpath(apyParFile.c_str(), NULL);
		if (!res_file) {
			write_log ("realpath error", EXIT_FAILURE);
			exit(EXIT_FAILURE);
		}

		ped.open(res_file, std::fstream::in);

		while (getline (ped,line)) {

			std::string delimiter = " ";

			size_t pos = 0;
			std::string token;
			while ((pos = line.find(delimiter)) != std::string::npos) {
				token = line.substr(0, pos);
				apy_data_list.push_back(token);
				line.erase(0, pos + delimiter.length());
			}
			/* get the last element of the string */
			apy_data_list.push_back(line);

            if (apy_data_list.size() < 3) {
                write_log (" There is not enough parameters in a data file !", -1);
                exit(-1);
            }
            else if (apy_data_list.size() == 3) {
                apyDatStr.file_ped = apy_data_list[0];
                apyDatStr.file_typd = apy_data_list[1];
                apyDatStr.file_gen = apy_data_list[2];
            }
            else if (apy_data_list.size() == 4) {
                apyDatStr.file_ped = apy_data_list[0];
                apyDatStr.file_typd = apy_data_list[1];
                apyDatStr.file_gen = apy_data_list[2];
                apyDatStr.inbreed = std::stoi(apy_data_list[3]);
            }
            else if (apy_data_list.size() == 5) {
                apyDatStr.file_ped = apy_data_list[0];
                apyDatStr.file_typd = apy_data_list[1];
                apyDatStr.file_gen = apy_data_list[2];
                apyDatStr.inbreed = std::stoi(apy_data_list[3]);
                apyDatStr.weight = std::stod(apy_data_list[4]);
            }
            else if (apy_data_list.size() == 6) {
                apyDatStr.file_ped = apy_data_list[0];
                apyDatStr.file_typd = apy_data_list[1];
                apyDatStr.file_gen = apy_data_list[2];
                apyDatStr.inbreed = std::stoi(apy_data_list[3]);
                apyDatStr.weight = std::stod(apy_data_list[4]);
                apyDatStr.adjust = apy_data_list[5];
            }
            else {
                apyDatStr.file_ped = apy_data_list[0];
                apyDatStr.file_typd = apy_data_list[1];
                apyDatStr.file_gen = apy_data_list[2];
                apyDatStr.inbreed = std::stoi(apy_data_list[3]);
                apyDatStr.weight = std::stod(apy_data_list[4]);
                apyDatStr.adjust = apy_data_list[5];
				apyDatStr.file_marker = apy_data_list[6];
            }

            apy_data_list.erase(apy_data_list.begin(), apy_data_list.end());
			
		}
		ped.close();

		return 0;
	}
	catch (std::exception const& e) {
		if (ped.eof())
			ped.close();
		else {
			write_log (e.what(), 1);
			return 10;
		}
	}
	catch (...) {
		return 1;
	}

	return 0;
}

//===============================================================================================================

int HMINTERFACE::getParameters() {

	std::string line;
	std::ifstream ped;

    std::vector<std::string> apy_data_list;

	ped.exceptions ( std::ifstream::failbit | std::ifstream::badbit );

	/* some default values */
	apyDatStr.file_log = "LOG_HMAT.log";
	apyDatStr.adjust = " ";
	apyDatStr.method = 0;
	apyDatStr.inbreed = 0;
	apyDatStr.weight = 0.0f;
	apyDatStr.gdiag = 0.0f;

	/* Keywords */
	std::string pedigree("$PEDIGREE");
	std::string genid("$GENID");
	std::string gmatrix("$GMATRIX");
	std::string marker("$MARKER");
	std::string adjust("$ADJUST");
	std::string weight("$WEIGHT");
	std::string gdiag("$GDIAG");
	std::string inbreeding("$INBREEDING");

	try {

		char *res_file = realpath(apyParFile.c_str(), NULL);
		if (!res_file) {
			write_log ("realpath error", EXIT_FAILURE);
			exit(EXIT_FAILURE);
		}

		ped.open(res_file, std::fstream::in);

		while (getline (ped,line)) {

			if ( line.compare(pedigree) == 0 ) {
				getline (ped,line);
				apyDatStr.file_ped = line;
				continue;
			}
			else if (line.compare(genid) == 0 ) {
				getline (ped,line);
				apyDatStr.file_typd = line;
				continue;
			}
			else if (line.compare(gmatrix) == 0 ) {
				getline (ped,line);
				apyDatStr.file_gen = line;
				continue;
			}
			else if (line.compare(marker) == 0 ) {
				getline (ped,line);
				apyDatStr.file_marker = line;
				continue;
			}
			else if (line.compare(adjust) == 0 ) {
				getline (ped,line);
				apyDatStr.adjust = line;
				continue;
			}
			else if (line.compare(weight) == 0 ) {
				getline (ped,line);
				apyDatStr.weight = std::stod( line );
				continue;
			}
			else if (line.compare(gdiag) == 0 ) {
				getline (ped,line);
				apyDatStr.gdiag = std::stod( line );
				continue;
			}
			else if (line.compare(inbreeding) == 0 ) {
				getline (ped,line);
				apyDatStr.inbreed = std::stoi( line );
				continue;
			}
			
		}
		ped.close();

		return 0;
	}
	catch (std::exception const& e) {
		if (ped.eof()){
			ped.close();

			/*std::cout<<"file_ped = "<< apyDatStr.file_ped << std::endl;
			std::cout<<"file_typd = "<< apyDatStr.file_typd << std::endl;
			std::cout<<"file_gen = "<< apyDatStr.file_gen << std::endl;
			std::cout<<"file_log = "<< apyDatStr.file_log << std::endl;
			std::cout<<"adjust = "<< apyDatStr.adjust << std::endl;
			std::cout<<"data_file = "<< apyParFile << std::endl;
			std::cout<<"file_marker = "<< apyDatStr.file_marker << std::endl;
			std::cout<<"method = "<< apyDatStr.method << std::endl;
			std::cout<<"inbreed = "<< apyDatStr.inbreed << std::endl;
			std::cout<<"weight = "<< apyDatStr.weight << std::endl;*/

			return 0;

		}
		else {
			write_log (e.what(), 1);
			return 10;
		}
	}
	catch (...) {
		return 1;
	}

	return 0;
}
//===============================================================================================================

bool HMINTERFACE::logfile_exists (const std::string& name) {
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

//===============================================================================================================

void HMINTERFACE::write_log (std::string location, int code) {

	std::fstream log;
	std::string log_file = "LOG_HMAT.log";
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
