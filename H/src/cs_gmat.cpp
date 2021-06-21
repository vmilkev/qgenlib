/*
*
*
*
*/

#include "cs_gmat.hpp"

//===============================================================================================================

GMAT::GMAT(const std::string& dataFile){
    marker_file = dataFile;
}

//===============================================================================================================

GMAT::~GMAT(){

    G.clear();
    
}

//===============================================================================================================

int GMAT::readSNP (){

    try {
        std::ifstream snpF;
        std::string line;
        std::vector<std::string> data_list;

        //snpF.exceptions(std::ifstream::failbit | std::ifstream::badbit);
        snpF.open( getFilePath( marker_file ), std::fstream::in );
        if (!snpF.good())
            throw 10;

        while (getline(snpF, line)){

            std::string delimiter = " ";
            size_t pos = 0;
			std::string token;
/*             while(  pos = line.find(delimiter) != std::string::npos ){
                
                std::cout<<"pos = "<<pos<<std::endl;
                pos++;
            }
 */
			while ((pos = line.find(delimiter)) != std::string::npos) {
                
                
                if (pos == 0)
                    token = " ";
                else
                    token = line.substr(0, pos);
				
                //std::cout<<"token = "<<token<<"; pos = "<<pos<<std::endl;
                line.erase(0, pos + delimiter.length());
                //std::cout<<"after erase pos = "<<pos<<std::endl;
            
                if (token.compare(delimiter) == 0){
                    //std::cout<<"(ins)"<<std::endl;                  
                    continue;
                }

				data_list.push_back(token);
                
				//line.erase(0, pos + delimiter.length());

                if ( data_list.size() == 2 )
                    break;
			}
			/* get the last element of the string */
			data_list.push_back(line);

            /* now we have got the SNP data for one ID */
            SNP[stoi(data_list[0])] = data_list[2];
            data_list.erase(data_list.begin(), data_list.end());
        }

        /* build the map: <index, ID> */
        size_t tmpInd = 0;
        for(const auto& snp: SNP){
            animID[tmpInd] = snp.first;
            tmpInd++;
        }

        // debug
/*         for (const auto& snp: SNP)
            std::cout<<"id, string = "<<snp.first<<"|"<<snp.second<<std::endl;
 */
        //exit(1);/**/

    }
    catch (int err) {
        write_log( "Error while opening SNP data file", 10 );
        return 1;
    }
    catch (...) {
        write_log( "Error while reading SNP data file", 10 );
        return 1;
    }

    return 0;
}

//===============================================================================================================

#ifdef _WIN64
    bool GMAT::logfile_exists( std::string& name ) {
        // ???
        return true;
    }
#endif

#ifdef linux
    bool GMAT::logfile_exists( std::string& name ) {
        struct stat buffer;
        return (stat (name.c_str(), &buffer) == 0);
    }
#endif

//===============================================================================================================

std::string GMAT::getFilePath( std::string &fileName ) {

	char *res = realpath(fileName.c_str(), NULL);
	if (!res) {
		perror("realpath");
		exit(EXIT_FAILURE);
	}
	return fileName = res;
}

//===============================================================================================================

void GMAT::write_log( std::string location, int code ) {

	std::fstream log;
    std::string log_file = "LOG_HMAT.log";
	bool exists = logfile_exists( log_file );
	if (exists)
		log.open (log_file, std::fstream::out | std::fstream::app);
	else log.open (log_file, std::fstream::out);

	auto end = std::chrono::system_clock::now();
	std::time_t end_time = std::chrono::system_clock::to_time_t(end);

	if (log) log << std::ctime(&end_time) << " :  Exception -> "<<location<<";  code -> "<<code<<std::endl;
	log.close();
}
//===============================================================================================================

int GMAT::makeG(){

    /*
        reading SNPs from the marker file and store it in a map as a string,
        where key is a animal ID and a value is a snp string (needs to be parsed)
    */
    if( readSNP() ){
        write_log( "Error while reading SNP from data file", -1 );
        return 1;
    }

    if( getZ() ){
        write_log( "Error while making Z matrix", -1 );
        return 1;
    }

    if( getG() ){
        write_log( "Error while making G matrix", -1 );
        return 1;
    }

    return 0;
}

//===============================================================================================================

int GMAT::parseString( std::string& snpStr, std::vector<int>& markers ){

    size_t sz = snpStr.length()+1;
    char * cstr = new char [snpStr.length()+1];
    std::strcpy (cstr, snpStr.c_str());
    for(size_t i = 0; i < sz; i++){
        if( isdigit(cstr[i]) ){
            markers.push_back( (int)cstr[i] - (int)48 );
            //std::cout<<"parsed = "<< (int)cstr[i] - (int)48 <<std::endl;
        }
    }

    return 0;
}

//===============================================================================================================

int GMAT::getZ(){

    try {

        /* get the number of SNPs */
        std::vector<int> tmpVect;
        auto it = SNP.begin();
        std::string tmpStr = it->second;
        parseString( tmpStr, tmpVect );
        snpNum = tmpVect.size();
        tmpVect.clear(); tmpVect.shrink_to_fit();
        tmpStr.clear(); tmpStr.shrink_to_fit();

        // debug
        //std::cout<<"snpNum = "<<snpNum<<std::endl;
        //exit(1);
        // end debug

        /* declare matrix M*/
        qgen::matrix <double> M( SNP.size(), snpNum );

        /* vector of SNPs frequences and missed values */
        std::vector <double> P(snpNum, 0.0);
        std::vector <int> missed(snpNum, 0);
        std::vector <double> missed2pq(SNP.size(), 0.0);

        /* map of missed values locations */
        std::vector < std::vector<int> > missedLocation;
        for(size_t i = 0; i < snpNum; i++)
            missedLocation.push_back(std::vector<int>());

        /* parse SNPs and fill matrix M*/
        size_t rowI = 0;
        for (auto const& e: SNP){

            std::vector<int> parsedMarkers;
            std::string strToParse = e.second;
            parseString( strToParse, parsedMarkers );

            for(size_t i = 0; i < parsedMarkers.size(); i++){
                M( rowI,i ) = static_cast <double>( parsedMarkers[i] );
                if ( parsedMarkers[i] != 0 && parsedMarkers[i] != 1 && parsedMarkers[i] != 2 ){
                    missed[i] += 1;
                    missedLocation[i].push_back(rowI);
                }
                else
                    P[i] += static_cast <double>( parsedMarkers[i] );
            }
            rowI++;

        }

        /* finish to calculate allele frequences, additionally accounting missing values */

        auto n_threads = std::thread::hardware_concurrency();
        auto block_size = static_cast<unsigned int> (P.size()/(n_threads));
        
        if(block_size < workload){
            block_size = P.size();
            n_threads = 1;
        }

    #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
            for(size_t i = 0; i < P.size(); i++){
                P[i] = P[i] / ( 2*(SNP.size() - missed[i]) );
            }

        Z.resize( SNP.size(),snpNum );

        for(size_t i = 0; i < SNP.size(); i++){
            for(size_t j = 0; j < snpNum; j++){
                Z(i,j) = M(i,j) - 2*P[j];
            }
        }

        /* modify Z matrix, so instead of missing values we put population average (0.0) */

    #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for(size_t i = 0; i < missedLocation.size(); i++){
            for(size_t j = 0; j < missedLocation[i].size(); j++){
                Z( missedLocation[i][j],i ) = 0.0;;
                missed2pq[ missedLocation[i][j] ] = missed2pq[ missedLocation[i][j] ] + 2*P[i]*(1.0-P[i]);
            }
        }

        freq = 0.0;
    #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for(size_t j = 0; j < P.size(); j++){
            freq += P[j]*(1.0-P[j]);
        }
        freq = 2*freq;

        /*
        double missedFreq = 0.0;
        for(size_t i = 0; i < missed2pq.size(); i++){
            missedFreq += missed2pq[i];
            std::cout<<"not adjusted missed 2pq = "<<missed2pq[i]<<std::endl;
        }
        missedFreq = missedFreq/SNP.size();
        std::cout<<"2pq = "<<freq<<"; av. missed 2pq = "<<missedFreq<<std::endl;       
        double missed_av = 0.0;
        for(auto i = 0; i < missed.size(); i++)
            missed_av += missed[i];
        std::cout<<"aver number of missed = "<<missed_av/SNP.size()<<std::endl;
        */

    #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for(size_t i = 0; i < SNP.size(); i++){
            missed2pq[i] = sqrt( freq/(freq - missed2pq[i]) );
        }

        /*
            After centering, adjust for missing markers for each animal;
            adjust for sqrt[sum of 2pq over all loci /sum of 2pq over non-missing loci.
        */
    #pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
        for(size_t i = 0; i < SNP.size(); i++){
            for(size_t j = 0; j < snpNum; j++){
                Z(i,j) = Z(i,j) * missed2pq[i];
            }
        }

        /*
        for(size_t i = 0; i < SNP.size(); i++){
            std::cout<<"missed2pq[i] = "<<missed2pq[i]<<std::endl;
        }
        for(size_t i = 0; i < SNP.size(); i++){
            for(size_t j = 0; j < snpNum; j++){
                std::cout<<"Z[i,j] adjusted = "<<Z(i,j)<<std::endl;
            }
            std::cout<<std::endl;
        }
        */

    }
    catch (...){
        write_log( "Error while making Z matrix", -1 );
        return 1;
    }
    return 0;

}

//===============================================================================================================

int GMAT::getG(){

    try {

        G = ( Z^2 ) * (1 / freq);

    }
    catch (...){
        write_log( "Error while making G matrix", -1 );
        return 1;
    }
    return 0;

}

//===============================================================================================================
/*
void GMAT::getMeans(){

    size_t ids = SNP.size();
    double _diag_mean = 0.0;
    double _ofd_mean = 0.0;

    #pragma omp parallel for reduction (+: _diag_mean, _ofd_mean)
        for(size_t i = 0; i < ids; i++) {
            for(size_t j = 0; j <= i; j++) {
                if (i == j)
                    _diag_mean = _diag_mean + G(i,j);
                else
                    _ofd_mean = _ofd_mean + G(i,j);
            }
        }

        all_mean = (_diag_mean + 2*_ofd_mean)/(ids*ids);
        diag_mean = _diag_mean/ids;
        ofd_mean = 2*_ofd_mean/(ids*ids-ids);

}
*/
//===============================================================================================================

void GMAT::getValues(std::vector <int>& row, std::vector <int>& col, std::vector <float>& val, double diag_val){

    for(size_t i = 0; i < animID.size(); i++) {
        //row.push_back(animID[i]);
        for(size_t j = 0; j <= i; j++) {
            col.push_back(animID[j]);
            row.push_back(animID[i]);
            if ( i == j ) {
                val.push_back( static_cast <float> ( G(i,j) + diag_val ) );
            }
            else {
                val.push_back( static_cast <float> ( G(i,j) ) );
            }
        }
    }

}

//===============================================================================================================

void GMAT::getIDs(std::vector <int>& IDs){

    for(const auto& id: animID)
        IDs.push_back(id.second);

}

//===============================================================================================================