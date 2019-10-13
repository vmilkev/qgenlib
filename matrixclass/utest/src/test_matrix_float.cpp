/* test_matrix_mytype.cpp 

 Tests of the matrix class from a static library.
 Class instance of type 'mytype'.
*/

#include "catch.hpp"
#include "../../mathlib/src/cs_matrix.hpp"

typedef float mytype;

TEST_CASE("Checking class constructors, type = float"){
    
    SECTION("Default constructor"){
        
        mklmatr::matrix <mytype> M;

        CHECK(M.size() == 0);
        CHECK(M.failbit == false);
        REQUIRE_FALSE(M.failbit);

        SECTION("Resize & Clear"){
            M.resize(10,10);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 100);
            M.resize(1,3);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 3);
            M.clear();
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 0);
            CHECK(M.failbit == false);
        }
    }

    SECTION("Constructor: M(val1, val2)"){
        
        mklmatr::matrix <mytype> M(3,5);

        REQUIRE_FALSE(M.failbit);
        
        CHECK(M.size() == 15);
        CHECK(M.failbit == false);
        
        SECTION("Resize & Clear"){
            M.resize(10,10);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 100);
            M.resize(0,0);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 0);
            M.clear();
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 0);
            CHECK(M.failbit == false);
        }
    }

    SECTION("Constructor: M(val1)"){
        
        mklmatr::matrix <mytype> M(6);

        REQUIRE_FALSE(M.failbit);
        
        CHECK(M.size() == 21);
        CHECK(M.failbit == false);

        SECTION("Resize & Clear"){
            
            M.resize(10,10);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 55);

            M.resize(1,1);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 1);

            M.clear();
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 0);

            M.resize(25);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 325);
            CHECK(M.failbit == false);
        }
    }

    SECTION("Constructor: M(char_var)"){
        
        mklmatr::matrix <mytype> M("s");

        REQUIRE_FALSE(M.failbit);

        CHECK(M.size() == 0);
        CHECK(M.failbit == false);

        SECTION("Resize & Clear"){
            M.resize(100,100);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 5050);
            M.resize(1,1);
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 1);
            M.clear();
            REQUIRE_FALSE(M.failbit);
            CHECK(M.size() == 0);
            CHECK(M.failbit == false);
        }
    }

    SECTION("Constructors & overloaded assignment operator"){
        
        mklmatr::matrix <mytype> M1;
        mklmatr::matrix <mytype> M2(5,3);
        mklmatr::matrix <mytype> M3(4,5);

        CHECK(M1.size() == 0);
        CHECK(M2.size() == 15);
        CHECK(M3.size() == 20);
        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(M3.failbit);

        SECTION("Assignment to uninitialised matrix"){
            M1 = M2;
            CHECK(M1.size() == 15);
            CHECK(M2.size() == 15);
            REQUIRE_FALSE(M1.failbit);
            REQUIRE_FALSE(M2.failbit);
        }

        SECTION("Assignment to initialised matrix"){
            M3 = M2;
            CHECK(M3.size() == 15);
            CHECK(M2.size() == 15);
            REQUIRE_FALSE(M3.failbit);
            REQUIRE_FALSE(M2.failbit);
        }
    }

    SECTION("Constructors & overloaded assignment operator: matrices with values"){
        
        mklmatr::matrix <mytype> M1;
        mklmatr::matrix <mytype> M2(5,3);
        mklmatr::matrix <mytype> M3(4,5);

        CHECK(M1.size() == 0);
        CHECK(M2.size() == 15);
        CHECK(M3.size() == 20);
        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(M3.failbit);

        SECTION("Assignment to uninitialised matrix"){
            M1 = M2;
            CHECK(M1.size() == 15);
            CHECK(M2.size() == 15);
            REQUIRE_FALSE(M1.failbit);
            REQUIRE_FALSE(M2.failbit);

            for (auto i = 0; i < M1.size(); i++)
                M1[i] = 5.0;

            M2 = M1;

            REQUIRE_FALSE(M1.failbit);
            REQUIRE_FALSE(M2.failbit);

            for (auto i = 0; i < M1.size(); i++)
                CHECK(M1[i] == M2[i]);

        }

        SECTION("Assignment to initialised matrix"){
            M3 = M2;
            CHECK(M3.size() == 15);
            CHECK(M2.size() == 15);
            REQUIRE_FALSE(M3.failbit);
            REQUIRE_FALSE(M2.failbit);

            for (auto i = 0; i < M1.size(); i++)
                M3[i] = 15.0;

            M2 = M3;

            REQUIRE_FALSE(M3.failbit);
            REQUIRE_FALSE(M2.failbit);

            for (auto i = 0; i < M1.size(); i++)
                CHECK(M3[i] == M2[i]);

        }

    }

    SECTION("Self-assignment"){

        mklmatr::matrix <mytype> M2(5,3);
        REQUIRE_FALSE(M2.failbit);
        M2 = M2;
        REQUIRE_FALSE(M2.failbit);
    }

}

TEST_CASE("Check access operators, type = float"){
    
    SECTION("In case of the default constructor"){
        
        mklmatr::matrix <mytype> M;

        SECTION("Checking at()"){
            
            M.at(10,10) = 0.8;

            CHECK(M.at(10,10) == Approx(0.8));
            CHECK(M.size() == 121);

            M.resize(1,3);
            M.at(0,2) = 0.2;

            CHECK(M.at(0,2) == Approx(0.2));
            CHECK(M.size() == 3);

            M.clear();
            M.at(1000, 200) = 1.2;

            CHECK(M.at(1000,200) == Approx(1.2));
            CHECK(M.size() == 201201);
            CHECK(M.capacity() == 201201);
            CHECK(M.failbit == false);           
        }

        SECTION("Checking operator()"){
            
            M.resize(10,10);
            M(9,9) = 0.9;
            CHECK(M(9,9) == Approx(0.9));
            CHECK(M(0,0) == Approx(0.0));
            CHECK(M.failbit == false);
        }
    }

    SECTION("In case of the constructor: M(val1, val2)"){
        
        mklmatr::matrix <mytype> M(11,11);

        SECTION("Checking at()"){
            
            CHECK(M.size() == 121);

            M.at(10,10) = 0.8;

            CHECK(M.at(10,10) == Approx(0.8));
            
            M.resize(1,3);
            M.at(0,2) = 0.2;

            CHECK(M.at(0,2) == Approx(0.2));
            CHECK(M.size() == 3);

            M.clear();
            M.at(1000, 200) = 1.2;

            CHECK(M.at(1000,200) == Approx(1.2));
            CHECK(M.size() == 201201);
            CHECK(M.capacity() == 201201);
            CHECK(M.failbit == false);           
        }

        SECTION("Checking operator()"){
            
            M(9,9) = 0.9;

            CHECK(M(9,9) == Approx(0.9));
            CHECK(M(0,0) == Approx(0.0));
            CHECK(M.failbit == false);
        }

    }

    SECTION("In case of the constructor: M(val1)"){
        
        mklmatr::matrix <mytype> M(11);

        SECTION("Checking at()"){
            
            CHECK(M.size() == 66);

            M.at(10,10) = 0.8;

            CHECK(M.at(10,10) == Approx(0.8));
            
            M.resize(3,1);
            M.at(0,0) = 0.5;

            CHECK(M.at(0,0) == Approx(0.5));
            CHECK(M.size() == 6);

            M.resize(1);

            CHECK(M.size() == 1);

            M.at(2,0) = 0.2;

            CHECK(M.at(2,0) == Approx(0.2));
            CHECK(M.size() == 6);

            M.clear();
            M.at(1000, 200) = 1.2;

            CHECK(M.at(1000,200) == Approx(1.2));
            CHECK(M.size() == 501501);
            CHECK(M.failbit == false);

            M.at(3, 2) = 1.2;

            CHECK(M.at(3,2) == Approx(1.2));
            CHECK(M.failbit == false);                    
        }

        SECTION("Checking operator()"){
            
            M(9,9) = 0.9;

            CHECK(M(9,9) == Approx(0.9));
            CHECK(M(0,0) == Approx(0.0));
            CHECK(M.failbit == false);
        }

    }

}

TEST_CASE("Testing scale() method, type = float"){
    
    SECTION("Scaling to 0.5 a square matrix"){

        mklmatr::matrix <mytype> M(3,3);
        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
        }

        M.scale(0.5);
        
        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            mytype j = i*0.5;
            CHECK(M[i] == Approx(j));
        }

    }

    SECTION("Scaling to -2.5 a square matrix"){

        mklmatr::matrix <mytype> M(3,3);
        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
        }

        M.scale(-2.5);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            mytype j = i*(-2.5);
            CHECK(M[i] == Approx(j));
        }

    }

    SECTION("Scaling to 0.5 a rectangular matrix"){

        mklmatr::matrix <mytype> M(3,5);
        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
        }

        M.scale(0.5);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            mytype j = i*0.5;
            CHECK(M[i] == Approx(j));
        }

    }

    SECTION("Scaling to 1.5 a rectangular matrix"){

        mklmatr::matrix <mytype> M(5,3);
        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
        }

        M.scale(1.5);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            mytype j = i*1.5;
            CHECK(M[i] == Approx(j));
        }

    }

    SECTION("Scaling to -2.5 a rectangular matrix"){

        mklmatr::matrix <mytype> M(13,5);
        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
        }

        M.scale(-2.5);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            mytype j = i*(-2.5);
            CHECK(M[i] == Approx(j));
        }

    }

}

TEST_CASE("Testing overloaded operator = , type = float"){
    
    SECTION("Checking an initialised square matrix"){
        
        mklmatr::matrix <mytype> M1(3,3);
        mklmatr::matrix <mytype> M2(3,3);

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 1;
        }
        M1.scale(5.0);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] > M2[i]);
        }

        M2 = M1;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }
    }

    SECTION("Checking an uninitialised square matrix"){
        
        mklmatr::matrix <mytype> M1(3,3);
        mklmatr::matrix <mytype> M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 1;
        }
        M1.scale(5.0);

        M2 = M1;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }
    }

    SECTION("Checking an initialised rectangular matrix (3,5)"){
        
        mklmatr::matrix <mytype> M1(3,5);
        mklmatr::matrix <mytype> M2(3,5);

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 1;
        }
        M1.scale(5.0);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] > M2[i]);
        }

        M2 = M1;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }
    }

    SECTION("Checking an uninitialised rectangular matrix (3,5)"){
        
        mklmatr::matrix <mytype> M1(3,5);
        mklmatr::matrix <mytype> M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 1;
        }

        M1.scale(5.0);

        M2 = M1;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }
    }

    SECTION("Checking a rectangular matrix (5,3)"){
        
        mklmatr::matrix <mytype> M1(5,3);
        mklmatr::matrix <mytype> M2(5,3);

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 1;
        }
        M1.scale(5.0);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] > M2[i]);
        }

        M2 = M1;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }
    }

    SECTION("Checking an initialised symmetrical (compact) matrix"){
        
        mklmatr::matrix <mytype> M1(3);
        mklmatr::matrix <mytype> M2(3);

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 1;
        }
        M1.scale(-5.0);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] < M2[i]);
        }

        M2 = M1;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }
    }

    SECTION("Checking an uninitialised symmetrical (compact) matrix"){
        
        mklmatr::matrix <mytype> M1(3);
        mklmatr::matrix <mytype> M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 1;
        }
        M1.scale(-5.0);

        M2 = M1;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for (auto i = 0; i < M1.size(); i++){
            CHECK(M1[i] == M2[i]);
        }
    }

    SECTION("Checking non equivalent matrices"){
        
        mklmatr::matrix <mytype> M1(3,3);
        mklmatr::matrix <mytype> M2(3);

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        CHECK(M1.size() > M2.size());

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 1;
        }
        M1.scale(-5.0);

        M2[0] = 1.0;

        for (auto i = 0; i < M2.size(); i++){
            CHECK(M1[i] != M2[i]);
        }

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
    }

}


TEST_CASE("Testing symtorec() method, type = float"){

    mklmatr::matrix <mytype> M(3);
    mytype m[] = {0, 1, 3, 1, 2, 4, 3, 4, 5};
    mytype m2[] = {-0, -1, -3, -1, -2, -4, -3, -4, -5};

    REQUIRE_FALSE(M.failbit);

    CHECK(M.size() == 6);

    SECTION("Positive matrix"){

        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
        }
        M.symtorec();
        
        REQUIRE_FALSE(M.failbit);

        CHECK(M.size() == 9);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }
    }

    SECTION("Negative matrix"){       

        for(auto i = 0; i < M.size(); i++){
            M[i] = -1*static_cast<mytype> (i);
        }
        M.symtorec();

        REQUIRE_FALSE(M.failbit);

        CHECK(M.size() == 9);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m2[i]));
        }
    }

}

TEST_CASE("Check rectosym() method, type = float"){
    
    mklmatr::matrix <mytype> M(3,3);
    CHECK(M.size() == 9);
    REQUIRE_FALSE(M.failbit);

    SECTION("Positive matrix"){
        
        mytype m[] = {0, 3, 4, 6, 7, 8};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
        }

        M.rectosym();

        REQUIRE_FALSE(M.failbit);

        CHECK(M.size() == 6);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }
    }

    SECTION("Negative matrix"){
        
        mytype m[] = {-0.0, -0.3, -0.4, -0.6, -0.7, -0.8};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i * static_cast <mytype> (-0.1);
        }

        M.rectosym();

        REQUIRE_FALSE(M.failbit);

        CHECK(M.size() == 6);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }
    }

}

TEST_CASE("Check transpose() method, type = float"){

    SECTION("Square matrix"){
        
        size_t dim = 3;
        mklmatr::matrix <mytype> M(dim,dim);

        REQUIRE_FALSE(M.failbit);

        mytype m[] = {-0.0, -0.3, -0.6, -0.1, -0.4, -0.7, -0.2, -0.5, -0.8};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i * static_cast <mytype> (-0.1);
        }

        for(auto i = 0; i < M.size(); i++){
            if( i%(dim-1) )
                CHECK(M[i] != Approx(m[i]));
        }

        M.transpose();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }
    }

    SECTION("Rectangular matrix (3,5)"){
        
        size_t dim1 = 3;
        size_t dim2 = 5;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        mytype m[] = {-0.0, -0.5, -1.0, -0.1, -0.6, -1.1, -0.2, -0.7, -1.2, -0.3, -0.8, -1.3, -0.4, -0.9, -1.4};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i * static_cast <mytype> (-0.1);
        }

        M.transpose();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }
    }

    SECTION("Rectangular matrix (5,3)"){
        
        size_t dim1 = 5;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        mytype m[] = {-0.0, -0.3, -0.6, -0.9, -1.2, -0.1, -0.4, -0.7, -1.0, -1.3, -0.2, -0.5, -0.8, -1.1, -1.4};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i * static_cast <mytype> (-0.1);
        }

        M.transpose();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }
    }

    SECTION("Symmetrical matrix"){
        
        size_t dim = 3;
        mklmatr::matrix <mytype> M(dim);

        REQUIRE_FALSE(M.failbit);

        mytype m[] = {-0.0, -0.1, -0.3, -0.1, -0.2, -0.4, -0.3, -0.4, -0.5};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i * static_cast <mytype> (-0.1);
        }

        M.transpose();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }
    }

}

TEST_CASE("Chacking fread()/fwrite() methods, type = float"){
    
    SECTION("Rectangular matrix (5,3)"){
        
        size_t dim1 = 5;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        mytype m[] = {-0.0, -0.3, -0.6, -0.9, -1.2, -0.1, -0.4, -0.7, -1.0, -1.3, -0.2, -0.5, -0.8, -1.1, -1.4};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i * static_cast <mytype> (-0.1);
        }

        M.transpose();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }

        M.fwrite();

        REQUIRE_FALSE(M.failbit);
        CHECK(M.failbit == false);
        REQUIRE(M.empty() == true);

        M.fread();

        REQUIRE_FALSE(M.failbit);
        CHECK(M.failbit == false);
        CHECK(M.size() == 15);
        CHECK(M.capacity() == 15);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }

    }

    SECTION("Rectangular matrix (3,5)"){
        
        size_t dim1 = 3;
        size_t dim2 = 5;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        mytype m[] = {-0.0, -0.5, -1.0, -0.1, -0.6, -1.1, -0.2, -0.7, -1.2, -0.3, -0.8, -1.3, -0.4, -0.9, -1.4};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i * static_cast <mytype> (-0.1);
        }

        M.transpose();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }

        M.fwrite();

        REQUIRE_FALSE(M.failbit);
        REQUIRE(M.empty() == true);
        
        M.fread();

        REQUIRE_FALSE(M.failbit);
        CHECK(M.size() == 15);
        CHECK(M.capacity() == 15);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }

    }

    SECTION("Symmetrical matrix (3,3)"){
        
        size_t dim = 3;
        mklmatr::matrix <mytype> M(dim);

        REQUIRE_FALSE(M.failbit);

        mytype m[] = {-0.0, -0.1, -0.3, -0.1, -0.2, -0.4, -0.3, -0.4, -0.5};

        for(auto i = 0; i < M.size(); i++){
            M[i] = i * static_cast <mytype> (-0.1);
        }

        M.transpose();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }

        M.fwrite();

        REQUIRE_FALSE(M.failbit);
        REQUIRE(M.empty() == true);

        M.fread();

        REQUIRE_FALSE(M.failbit);
        CHECK(M.size() == 9);
        CHECK(M.capacity() == 9);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }

    }
}

TEST_CASE("Checking the invert() method, type = float"){
    
    SECTION("Positive square matrix (3,3)"){
        
        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = i;
        }
        M[0] = 3;

        mytype m[] = {0.333333, -0.66666, 0.333333, -0.66666, -1.333333, 1.0, 0.333333, 1.66666, -1.0};

        M.invert();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK( M[i] == Approx( m[i]) );
        }

        CHECK_FALSE(M.failbit);
    }

    SECTION("Negative square matrix (3,3)"){
        
        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = -1.0*i;
        }
        M[0] = -3;

        mytype m[] = {-0.333333, 0.66666, -0.333333, 0.66666, 1.333333, -1.0, -0.333333, -1.66666, 1.0};

        M.invert();

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK( M[i] == Approx( m[i]) );
        }

        CHECK_FALSE(M.failbit);
    }

    SECTION("Symmetrical matrix (3,3)"){
        
        size_t dim = 3;
        mklmatr::matrix <mytype> M(dim);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = i;
        }
        M[0] = 3.0; M[2] = 6.0;

        mytype m[] = {1.99999, 0.99999, 0.857142, -1.99999, -1.2857142, 2.4285714};

        M.invert();

        REQUIRE_FALSE(M.failbit);

        CHECK(M.size() == 6);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }

        CHECK_FALSE(M.failbit);

    }

}

TEST_CASE("Checking overloaded operator ^, type = float"){
    
    SECTION("Positive square matrix (3,3) -> checking inversion with self-assignment"){
        
        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = i;
        }
        M[0] = 3;

        mytype m[] = {0.333333, -0.66666, 0.333333, -0.66666, -1.333333, 1.0, 0.333333, 1.66666, -1.0};

        M = M^(-1);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK( M[i] == Approx( m[i]) );
        }

        CHECK_FALSE(M.failbit);
    }
    SECTION("Positive square matrix (3,3) -> checking inversion without self-assignment"){
        
        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);
        mklmatr::matrix <mytype> M2;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = i;
        }
        M[0] = 3;

        mytype m[] = {0.333333, -0.66666, 0.333333, -0.66666, -1.333333, 1.0, 0.333333, 1.66666, -1.0};

        M2 = M^(-1);

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK( M2[i] == Approx( m[i]) );
        }

        CHECK_FALSE(M2.failbit);
    }

    SECTION("Negative square matrix (3,3) -> checking inversion"){
        
        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = -1.0*i;
        }
        M[0] = -3;

        mytype m[] = {-0.333333, 0.66666, -0.333333, 0.66666, 1.333333, -1.0, -0.333333, -1.66666, 1.0};

        M = M^(-1);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK( M[i] == Approx( m[i]) );
        }

        CHECK_FALSE(M.failbit);
    }

    SECTION("Symmetrical matrix (3,3) -> checking inversion"){
        
        size_t dim = 3;
        mklmatr::matrix <mytype> M(dim);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = i;
        }
        M[0] = 3.0; M[2] = 6.0;

        mytype m[] = {1.99999, 0.99999, 0.857142, -1.99999, -1.2857142, 2.4285714};

        M = M^(-1);

        REQUIRE_FALSE(M.failbit);

        CHECK(M.size() == 6);

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == Approx(m[i]));
        }

        CHECK_FALSE(M.failbit);

    }

    SECTION("Symmetrical matrix (3,3) -> checking A*A'"){
        
        size_t dim = 3;
        mklmatr::matrix <mytype> M(dim);
        mklmatr::matrix <mytype> res;

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = i;
        }
        M[0] = 3.0; M[2] = 6.0;

        mytype m[] = {19, 21, 28, 21, 53, 47, 28, 47, 50};
        mytype m0[] = {3,1,6,3,4,5};

        res = M^2;

        REQUIRE_FALSE(res.failbit);
        REQUIRE_FALSE(M.failbit);

        CHECK(res.size() == 9);

        for(auto i = 0; i < res.size(); i++){
            CHECK(res[i] == Approx(m[i]));
        }

        for(auto i = 0; i < M.size(); i++){
            CHECK(M[i] == m0[i]);
        }

    }

    SECTION("Negative square matrix (3,3) -> checking A*A'"){
        
        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = -1.0*i;
        }
        M[0] = -3;

        mytype m[] = {14, 23, 41, 23, 50, 86, 41, 86, 149};

        M = M^2;

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK( M[i] == Approx( m[i]) );
        }

        CHECK_FALSE(M.failbit);
    }

    SECTION("Part-negative square matrix (3,3) -> checking A*A'"){
        
        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = 1.0*i;
        }
        M[0] = -3;

        mytype m[] = {14, 5, 5, 5, 50, 86, 5, 86, 149};

        M = M^2;

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK( M[i] == Approx( m[i]) );
        }

        CHECK_FALSE(M.failbit);
    }

    SECTION("Part-negative square matrix (3,3) -> checking A*A' without self-assignment"){
        
        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M(dim1,dim2);
        mklmatr::matrix <mytype> M2;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 1; i < M.size(); i++){
            M[i] = 1.0*i;
        }
        M[0] = -3;

        mytype m[] = {14, 5, 5, 5, 50, 86, 5, 86, 149};

        M2 = M^2;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 0; i < M.size(); i++){
            CHECK( M2[i] == Approx( m[i]) );
        }

        CHECK_FALSE(M.failbit);
    }

}

TEST_CASE("Checking overloaded operator *, type = float"){

        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M1(dim1,dim2);
        mklmatr::matrix <mytype> M2(dim1,dim2);
        mklmatr::matrix <mytype> res;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        mytype m1[] = {42, 45, 48, 150, 162, 174, 258, 279, 300};
        mytype m2[] = {-42, -45, -48, -150, -162, -174, -258, -279, -300};

    SECTION("Positive square matrices: (3,3)*(3,3): result is in uninitialised matrix"){

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = i;
        }

        for(auto i = 0; i < M2.size(); i++){
            M2[i] = i+9;
        }

        res = M1*M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        if( !res.failbit ){
            CHECK(res.size() == 9);
            for(auto i = 0; i < res.size(); i++){
                CHECK( res[i] == Approx( m1[i]) );
            }
        }
    }

    SECTION("Positive square matrices: (3,3)*(3,3): result as self-assignment"){

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = i;
        }

        for(auto i = 0; i < M2.size(); i++){
            M2[i] = i+9;
        }

        M1 = M1*M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        if( !M1.failbit ){
            CHECK(M1.size() == 9);
            for(auto i = 0; i < M1.size(); i++){
                CHECK( M1[i] == Approx( m1[i]) );
            }
        }
    }
}

TEST_CASE("(2) Checking overloaded operator *, type = float"){

        size_t dim1 = 3;
        size_t dim2 = 3;
        mklmatr::matrix <mytype> M1(dim1,dim2);
        mklmatr::matrix <mytype> M2(dim1,dim2);
        mklmatr::matrix <mytype> res;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        mytype m2[] = {-42, -45, -48, -150, -162, -174, -258, -279, -300};

    SECTION("Square matrices: -(3,3)*(3,3): result is in uninitialised matrix"){
        
        for(auto i = 0; i < M1.size(); i++){
            M1[i] = -1.0*i;
        }

        for(auto i = 0; i < M2.size(); i++){
            M2[i] = i+9;
        }

        res = M1*M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        for(auto i = 0; i < res.size(); i++){
            CHECK( res[i] == Approx( m2[i]) );
        }
    }

    SECTION("Square matrices: -(3,3)*(3,3): result as self-assignment"){
        
        for(auto i = 0; i < M1.size(); i++){
            M1[i] = -1.0*i;
        }

        for(auto i = 0; i < M2.size(); i++){
            M2[i] = i+9;
        }

        M2 = M1*M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 0; i < M2.size(); i++){
            CHECK( M2[i] == Approx( m2[i]) );
        }
    }
}


TEST_CASE("(3) Checking overloaded operator *, type = float"){

        size_t dim1 = 3;
        size_t dim2 = 5;
        mklmatr::matrix <mytype> M1(dim1,dim2);
        mklmatr::matrix <mytype> M2;
        mklmatr::matrix <mytype> res;
        mklmatr::matrix <mytype> M(5,3);

        mytype m[] = {30, 80, 130, 80, 255, 430, 130, 430, 730};

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

    SECTION("Rectangular square matrices: A*A'"){

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = i;
        }

        M2 = M1;
        M2.transpose();

        res = M1*M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        CHECK(res.size() == 9);

        for(auto i = 0; i < res.size(); i++){
            CHECK( res[i] == Approx( m[i]) );
        }

        CHECK_FALSE(res.failbit);

        res.clear();
        M2.clear();
        M2 = M1^"T";
        res = M1*M2;

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        CHECK(res.size() == 9);

        for(auto i = 0; i < res.size(); i++){
            CHECK( res[i] == Approx( m[i]) );
        }

        CHECK_FALSE(res.failbit);

        res.clear();
        M2.clear();
        res = M1*(M1^"T");

        REQUIRE_FALSE(M1.failbit);
        REQUIRE_FALSE(res.failbit);

        CHECK(res.size() == 9);

        for(auto i = 0; i < res.size(); i++){
            CHECK( res[i] == Approx( m[i]) );
        }

        CHECK_FALSE(res.failbit);

    }

    SECTION("Rectangular matrix multiplied by positive scalar"){
        
        mytype val = 1.5;
        mytype m2[15];

        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
            m2[i] = i*val;
        }

        res = M * val;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(res.failbit);

        M = M * val;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(res.failbit);

        for(auto i = 0; i < res.size(); i++){
            CHECK( res[i] == Approx( m2[i]) );
        }
       
        for(auto i = 0; i < res.size(); i++){
            CHECK( M[i] == Approx( m2[i]) );
        }

        CHECK_FALSE(res.failbit);
    }

    SECTION("Rectangular matrix multiplied by negative scalar"){
        
        mytype val = -1.5;
        mytype m3[15];

        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
            m3[i] = i*val;
        }

        res = M * val;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(res.failbit);

        M = M * val;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(res.failbit);

        for(auto i = 0; i < res.size(); i++){
            CHECK( res[i] == Approx( m3[i]) );
        }
       
        for(auto i = 0; i < res.size(); i++){
            CHECK( M[i] == Approx( m3[i]) );
        }

        CHECK_FALSE(res.failbit);
    }

}

TEST_CASE("Checking overloaded +/- operators, type = float"){

        mklmatr::matrix <mytype> M(3,5);
        mklmatr::matrix <mytype> M2;
        mklmatr::matrix <mytype> res;
        std::vector <mytype> m(15, 0.0);

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

    SECTION("Adding two positive matrices"){

        for(auto i = 0; i < M.size(); i++){
            M[i] = i;
            m[i] = i*2.0;
        }
        
        M2 = M;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        CHECK(M.size() == 15);
        CHECK(M2.size() == 15);
        CHECK(res.size() == 0);

        res = M + M2;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        if ( !res.failbit ) {
            CHECK(res.size() == 15);
            for(auto i = 0; i < res.size(); i++){
                CHECK( res[i] == Approx( m[i]) );
            }
        }
        else {
            CAPTURE(res.failinfo);
        }
    }
}

TEST_CASE("(2) Checking overloaded +/- operators, type = float"){

        mklmatr::matrix <mytype> M(3,5);
        mklmatr::matrix <mytype> M2;
        mklmatr::matrix <mytype> res;
        std::vector <mytype> m(15, 0.0);
        mklmatr::matrix <mytype> M3(3,5);

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

    SECTION("Adding two negative matrices"){
        
        M.resize(5,3);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            M[i] = -1.0*i;
            m[i] = -1.0*i*2.0;
        }

        REQUIRE_FALSE(M.failbit);

        res = M + M;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(res.failbit);

        if ( !res.failbit ) {
            CHECK(res.size() == 15);
            for(auto i = 0; i < res.size(); i++){
                CHECK( res[i] == Approx( m[i]) );
            }
        }
        else {
            CAPTURE(res.failinfo);
        }
    }

    SECTION("Adding negative and positive matrices"){
        
        M.resize(5,3);
        M2.resize(5,3);

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);

        for(auto i = 0; i < M.size(); i++){
            M[i] = 1.0*i;
            m[i] = -1.0*i;
        }

        M2 = M*(-2.0);

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);

        res = M + M2;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        if ( !res.failbit ) {
            CHECK(res.size() == 15);
            for(auto i = 0; i < res.size(); i++){
                CHECK( res[i] == Approx( m[i]) );
            }
        }
        else {
            CAPTURE(res.failinfo);
        }
        
    }

    SECTION("Substitution two positive matrices"){
        
        for(auto i = 0; i < M.size(); i++){
            M[i] = i*2.0;
            m[i] = i;
        }

        M2 = M*0.5;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        res = M - M2;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);
        REQUIRE_FALSE(res.failbit);

        if ( !res.failbit ) {
            CHECK(res.size() == 15);
            for(auto i = 0; i < res.size(); i++){
                CHECK( res[i] == Approx( m[i]) );
            }
        }
        else {
            CAPTURE(res.failinfo);
        }
    }

    SECTION("Substitution two negative matrices"){
        
        M.resize(5,3);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            M[i] = -1.0*i;
            m[i] = i;
        }

        M2 = M*2.0;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);

        res = M - M2;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(res.failbit);

        if ( !res.failbit ) {
            CHECK(res.size() == 15);
            for(auto i = 0; i < res.size(); i++){
                CHECK( res[i] == Approx( m[i]) );
            }
        }
        else {
            CAPTURE(res.failinfo);
        }
    }

    SECTION("Substitution negative and positive matrices"){
        
        M.resize(5,3);

        REQUIRE_FALSE(M.failbit);

        for(auto i = 0; i < M.size(); i++){
            M[i] = -1.0*i;
            m[i] = -1.0*i*3.0;
        }

        M2 = M*(-2.0);

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);

        res = M - M2;

        REQUIRE_FALSE(M.failbit);
        REQUIRE_FALSE(M2.failbit);
        CHECK_FALSE(res.failbit);
        
        if ( !res.failbit ) {
            CHECK(res.size() == 15);
            for(auto i = 0; i < res.size(); i++){
                CHECK( res[i] == Approx( m[i]) );
            }
        }
        else {
            CAPTURE(res.failinfo);
        }
    }

}

TEST_CASE ("Checking dual direction overloaded operator +, type = float"){
    
    mklmatr::matrix <mytype> M1(3,3);
    mklmatr::matrix <mytype> M2(3,3);
    mklmatr::matrix <mytype> res;

    SECTION("Left/right addition"){

        for(auto i = 0; i < M1.size(); i++)
            M1[i] = 2.0;
        
        res = 2.0+M2;

        REQUIRE_FALSE(res.failbit);
        
        for (auto i = 0; i < res.size(); i++)
            CHECK(res[i] == M1[i]);

        res.clear();

        REQUIRE_FALSE(res.failbit);

        res = M2 + 2.0;

        REQUIRE_FALSE(res.failbit);

        for (auto i = 0; i < res.size(); i++)
            CHECK(res[i] == M1[i]);

    }
}

TEST_CASE ("Checking dual direction overloaded operators * , type = float"){
    
    mklmatr::matrix <mytype> M1(3,3);
    mklmatr::matrix <mytype> M2(3,3);
    mklmatr::matrix <mytype> res;

    SECTION("Positive left/right multiplication"){

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 2.0;
            M2[i] = 1.0;
        }
        
        res = 2.0 * M2;

        REQUIRE_FALSE(res.failbit);
        
        for (auto i = 0; i < res.size(); i++)
            CHECK(res[i] == M1[i]);

        res.clear();

        REQUIRE_FALSE(res.failbit);

        res = M2 * 2.0;

        REQUIRE_FALSE(res.failbit);

        for (auto i = 0; i < res.size(); i++)
            CHECK(res[i] == M1[i]);

    }

    SECTION("Negative left/right multiplication"){

        for(auto i = 0; i < M1.size(); i++){
            M1[i] = 2.0;
            M2[i] = -1.0;
        }
        
        res = -2.0 * M2;

        REQUIRE_FALSE(res.failbit);
        
        for (auto i = 0; i < res.size(); i++)
            CHECK(res[i] == M1[i]);

        res.clear();

        REQUIRE_FALSE(res.failbit);

        M2 = M2*(-1);

        REQUIRE_FALSE(M2.failbit);

        res = -1.0 * M2 * (-2.0);

        REQUIRE_FALSE(res.failbit);

        for (auto i = 0; i < res.size(); i++)
            CHECK(res[i] == M1[i]);

    }
}

TEST_CASE("Final integration test, type = float"){
    mklmatr::matrix <mytype> M;
    mklmatr::matrix <mytype> res;
    mklmatr::matrix <mytype> res_tr;
    mklmatr::matrix <mytype> res2;
    mklmatr::matrix <mytype> tmp_res;
    mklmatr::matrix <mytype> tmp_res2;

    SECTION("Rectangular matrix"){

        mytype m_res[] = {-2.08, -7.76, -3.08, -6.76};
        mytype m_sqv[] = {91, 217, 217, 559};
        mytype m0[] = {1,2,3,4,5,6,7,8,9,10,11,12};
        mytype m_inv[] = {0.147883, -0.057407, -0.057407, 0.024074};

        M.resize(2,6);

        for (auto i = 0; i < M.size(); i++){
            M[i] = static_cast <mytype> (i+1);
            CHECK (M[i] == m0[i]);
        }

        tmp_res = M^2;

        REQUIRE_FALSE(tmp_res.failbit);

        for (auto i = 0; i < tmp_res.size(); i++)
            CHECK (tmp_res[i] == m_sqv[i]);

        tmp_res.transpose();

        REQUIRE_FALSE(tmp_res.failbit);

        //tmp_res2 = (M^2)^-1;
        tmp_res2 = M^2;

        REQUIRE_FALSE(tmp_res2.failbit);

        for (auto i = 0; i < tmp_res2.size(); i++)
            CHECK (tmp_res2[i] == m_sqv[i]);

        tmp_res2.invert();

        REQUIRE_FALSE(tmp_res2.failbit);

        for (auto i = 0; i < tmp_res2.size(); i++)
            CHECK( Approx(tmp_res2[i]) == (m_inv[i]) );    

        res = ( ( (M^2)^(-1) ) - 0.01 )*tmp_res;

        REQUIRE_FALSE(res.failbit);

        res_tr = (M^2)^"T";

        REQUIRE_FALSE(res_tr.failbit);

        for (auto i = 0; i < res_tr.size(); i++)
            CHECK (res_tr[i] == m_sqv[i]);

        res2 = ( (M^-2) - 0.01 )*( (M^2)^"T" );

        REQUIRE_FALSE(res2.failbit);

        for (auto i = 0; i < res.size(); i++)
            CHECK( res[i] == Approx(m_res[i]) );

        for (auto i = 0; i < res.size(); i++)
            CHECK( res[i] == res2[i] );
    }

    SECTION("Square symmetrical matrix"){

        M.resize(4);
        REQUIRE_FALSE(M.failbit);

        REQUIRE( M.size() == 10 );

        mklmatr::matrix <mytype> b(4,2);
        mytype m0[] = {1,2,3,4,5,6,7,8,9,10};
        mytype m0_tr[] = {1,2,4,7,2,3,5,8,4,5,6,9,7,8,9,10};
        mytype b0[] = {1,2,3,4,5,6,7,8};
        mytype m_res[] = {-3.149999,-2.909999,3.539999,3.269999,-0.249999,-0.229999,-0.399999,-0.369999};
        mytype m_sqv[] = {70,84,101,129,84,102,125,163,101,125,158,212,129,163,212,294};
        mytype m_inv[] = {290.99,-327.9999,24.4999,36.4999,-327.9999,370.9999,-29.4999,
                          -40.4999,24.4999,-29.4999,4.99999,1.99999,36.4999,-40.4999,1.99999,4.99999};
        
        for (auto i = 0; i < M.size(); i++){
            M[i] = static_cast <mytype> (i+1);
            CHECK (M[i] == m0[i]);
        }

        for (auto i = 0; i < b.size(); i++){
            b[i] = static_cast <mytype> (i+1);
            CHECK (b[i] == b0[i]);
        }

        tmp_res = M^2;

        REQUIRE_FALSE(tmp_res.failbit);

        for (auto i = 0; i < tmp_res.size(); i++)
            CHECK (tmp_res[i] == m_sqv[i]);

        tmp_res.invert();
        tmp_res = tmp_res;

        REQUIRE_FALSE(tmp_res.failbit);

        for (auto i = 0; i < tmp_res.size(); i++)
            CHECK ( (tmp_res[i]) == Approx(m_inv[i]).margin(0.1) );

        res = (M^-2)*b*0.01;

        REQUIRE_FALSE(res.failbit);

        for (auto i = 0; i < res.size(); i++)
            CHECK( res[i] == Approx(m_res[i]).margin(0.01) );

        for (auto i = 0; i < M.size(); i++){
            CHECK (M[i] == m0[i]);
        }

        res2.clear();

        REQUIRE_FALSE(res2.failbit);

        res2 = M^"T";

        REQUIRE_FALSE(res2.failbit);

        CHECK(res2.size() == 16);

        for (auto i = 0; i < res2.size(); i++)
            CHECK( res2[i] == Approx(m0_tr[i]) );

        M.symtorec();

        res2 = M*(M^"T");

        REQUIRE_FALSE(res2.failbit);

        for (auto i = 0; i < res2.size(); i++)
            CHECK( res2[i] == Approx(m_sqv[i]) );

        res2 = ( ( M*(M^"T") )^-1 )*b*0.01;

        REQUIRE_FALSE(res2.failbit);

        for (auto i = 0; i < res2.size(); i++)
            CHECK( res2[i] == Approx(m_res[i]).margin(0.01) );

    }
}