/* catchmain.cpp

 This is the file dedicated to compile the
 source code of Catch itself.
*/

/* Catch provides main(): */
#define CATCH_CONFIG_MAIN

#include "catch.hpp"

TEST_CASE( "1: All test cases reside in other .cpp files (empty)", "[multi-file:1]" ) {
}
