/*
*
*
*
*/

#include "cs_interface.hpp"

int main(int argc, char* argv[]) {

    HMINTERFACE apy;

    if (argc < 2) {            
        apy.write_log ("There are not enough arguments passed to the executable !", -1);
        return 1;
    }

    apy.get_ParamFileName(argv[1]);
    		
	if ( apy.getParameters()/*read_apy_input_data()*/ != 0 ) {
        apy.write_log ("ERROR while reading the input data from parameters file !", -1);
        return 1;
    }

    //apy.getGmatr();
    //exit(1);

	if ( apy.getHmatr() != 0 ) {
        apy.write_log ("ERROR while calling getHmatr !", -1);
        return 1;
    }

    return 0;
}