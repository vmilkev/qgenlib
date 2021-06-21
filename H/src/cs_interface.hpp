/*
*
*
*
*/

#pragma once

#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include "cs_apy.hpp"
#include "cs_gmat.hpp"

class HMINTERFACE {

    private:

        struct INDATA_APY {
            std::string file_ped;
            std::string file_typd;
            std::string file_gen;
            std::string file_log;
            std::string adjust;
            std::string data_file;
            std::string file_marker;
            int method;
            int inbreed;
            double weight;
            double gdiag;
        };

        INDATA_APY apyDatStr;
        std::string apyParFile;
        bool logfile_exists (const std::string& name);
        

    public:

        HMINTERFACE();
        int read_apy_input_data();
        void get_ParamFileName(const std::string& data_file);
        int getHmatr();
        void write_log (std::string location, int code);
        int getGmatr();

        int getParameters();

};