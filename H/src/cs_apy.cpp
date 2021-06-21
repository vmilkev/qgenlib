/*
 * cs_apy.cpp
 *
 *  Created on: Aug 7, 2017
 *      Author: vimi
 */

#include "cs_apy.hpp"


//===============================================================================================================

APY::APY(std::string f1, std::string f2, std::string f3, std::string f4, int m, bool in, double wt, std::string adj, std::string mark, double d_val)/* : HOUSEKEEPER(f1, f2, f3, f4, m, in, wt, adj) */{

	ped_file = f1;
	typed_file = f2;
	gmat_file = f3;
	log_file = f4;
	whichMethod = m;
	inbreeding = in;
	wa = wt;
	adjust = adj;
	marker_file = mark;

	gdiag_val = d_val;

	std::fstream log;
	std::string log_file("param"); /* optional parameters file */

	struct stat buffer;
	bool exists = (stat (log_file.c_str(), &buffer) == 0); /* check if the file exists */

	if (exists) {
		log.open (log_file, std::fstream::in);
		if (log) {
			log >> isH;
			log >> isHtxt;
			log >> tol_a22inv;
			log >> tol_Hsub;
			log >> save;
		}
		log.close();
	}
	else {
		isH = false;
		isHtxt = false;
		tol_a22inv = 0.0;
		tol_Hsub = 0.0;
		save = false;
	}

	if( const char* env_p = std::getenv("OMP_NUM_THREADS") ) {
		threadsAvail = std::max(atoi(env_p), 1);
		omp_set_num_threads(threadsAvail);
		n_threads = threadsAvail;
	}
	else {
		n_threads = std::thread::hardware_concurrency();
		omp_set_num_threads(n_threads);
	}

}

//===============================================================================================================

APY::~APY(){

}

//===============================================================================================================

int APY::run_apy () {

	int status = 0;
	std::string where("run_apy");

	try {

		std::string str_adjust("G-ADJUST");

		bool if_adjust = true;

		if (str_adjust.compare(adjust) != 0) if_adjust = false;

		std::string filename1 = "A22.bin";
		std::fstream fA22(filename1, fA22.binary | fA22.trunc | fA22.in | fA22.out);
		fA22.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fA22.is_open()) throw 51;

		std::string filename2 = "A22_inv.bin";
		std::fstream fA22i(filename2, fA22i.binary | fA22i.trunc | fA22i.in | fA22i.out);
		fA22i.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fA22i.is_open()) throw 52;

		std::string filename3 = "Af_inv.bin";
		std::fstream fAfi(filename3, fAfi.binary | fAfi.trunc | fAfi.in | fAfi.out);
		fAfi.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fAfi.is_open()) throw 53;

		std::string filename4 = "Ar_inv.bin";
		std::fstream fAri(filename4, fAri.binary | fAri.trunc | fAri.in | fAri.out);
		fAri.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fAri.is_open()) throw 54;

		std::string filename5 = "Gcc_inv.bin";
		std::fstream fGcci(filename5, fGcci.binary | fGcci.trunc | fGcci.in | fGcci.out);
		fGcci.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fGcci.is_open()) throw 55;

		std::string filename6 = "Gnc.bin";
		std::fstream fGnc(filename6, fGnc.binary | fGnc.trunc | fGnc.in | fGnc.out);
		fGnc.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fGnc.is_open()) throw 56;

		std::string filename7 = "Gnn_inv.bin";
		std::fstream fGnni(filename7, fGnni.binary | fGnni.trunc | fGnni.in | fGnni.out);
		fGnni.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fGnni.is_open()) throw 57;

		std::string filename8 = "Gcn.bin";
		std::fstream fGcn(filename8, fGcn.binary | fGcn.trunc | fGcn.in | fGcn.out);
		fGcn.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fGcn.is_open()) throw 58;

		std::string filename9 = "G12.bin";
		std::fstream fG12(filename9, fG12.binary | fG12.trunc | fG12.in | fG12.out);
		fG12.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fG12.is_open()) throw 59;

		std::string filename10 = "G_inv.bin";
		std::fstream fGi(filename10, fGi.binary | fGi.trunc | fGi.in | fGi.out);
		fGi.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fGi.is_open()) throw 61;

		std::string filename11 = "a11.bin";
		std::fstream fa11(filename11, fa11.binary | fa11.trunc | fa11.in | fa11.out);
		fa11.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fa11.is_open()) throw 62;

		std::string filename12 = "a21.bin";
		std::fstream fa21(filename12, fa21.binary | fa21.trunc | fa21.in | fa21.out);
		fa21.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fa21.is_open()) throw 63;

		std::string filename13 = "a22.bin";
		std::fstream fa22(filename13, fa22.binary | fa22.trunc | fa22.in | fa22.out);
		fa22.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fa22.is_open()) throw 64;

		std::string filename14 = "a22i.bin";
		std::fstream fa22i(filename14, fa22i.binary | fa22i.trunc | fa22i.in | fa22i.out);
		fa22i.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fa22i.is_open()) throw 65;

		/* for debugging */
		/*
		std::string filename14_prm = "a22i_prm.bin";
		std::fstream fa22i_prm(filename14_prm, fa22i_prm.binary | fa22i_prm.trunc | fa22i_prm.in | fa22i_prm.out);
		fa22i_prm.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fa22i_prm.is_open()) throw 652;
		 */
		/* -------------*/

		std::string filename15 = "a12.bin";
		std::fstream fa12(filename15, fa12.binary | fa12.trunc | fa12.in | fa12.out);
		fa12.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
		if (!fa12.is_open()) throw 66;


		MKL_INT64 AllocatedBytes;
		int N_AllocatedBuffers;

		runtimelog (true, "Data from the optional parameters file: ", 0, true, false);
		runtimelog (false, "NOTE, for the logical parameters, if 0 - FALSE, if 1 - TRUE", static_cast<double>(save), false, false);
		runtimelog (false, "1) produce full H_inv in addition to H_sub", static_cast<double>(isH), false, false);
		runtimelog (false, "2) produce ASCII file (H_sub.dat) in addition to H_sub.bin", static_cast<double>(isHtxt), false, false);
		runtimelog (false, "3) threshold to zero out the elements in a22inv matrix", static_cast<double>(tol_a22inv), false, false);
		runtimelog (false, "4) threshold to zero out the elements in H_sub matrix", static_cast<double>(tol_Hsub), false, false);
		runtimelog (false, "5) use compact matrix storage for inversion ('save' option)", static_cast<double>(save), false, false);

		runtimelog (false, "6) the value added to a diagonal elements of G matrix", static_cast<double>(gdiag_val), false, false);

		runtimelog (true, "Reading pedigree from file ...", 0, false, false);

		status = read_ped();
		if (status != 0) throw status;

		runtimelog (false, "number of unique pedIDs", static_cast<double>(pedID.size()), false, false);
		runtimelog (false, "number of elements of full pedigree", static_cast<double>(f_ped.size()), false, false);

		runtimelog (true, "Reading genotyped file ...", 0, false, false);

		status = read_gtypd();
		if (status != 0) throw status;

		runtimelog (false, "number of unique genotyped IDs", static_cast<double>(genotypedID.size()), false, false);
		runtimelog (false, "number of unique core IDs", static_cast<double>(coreID.size()), false, false);

		runtimelog (true, "Tracing-back the pedigree ...", 0, false, false);

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		  std::ofstream myfile ("pedigree_out.txt");
		  if (myfile.is_open())
		  {
				for(auto const& elem : f_ped) {
					myfile <<elem.first.day <<" "<<elem.first.id<< " "<<elem.second.id_1<<" "<<elem.second.id_2<<"\n";
				}

		    myfile.close();
		  }
		  else std::cout << "Unable to open file";
		  //----- End DEBUG PRINT ---------------------------------------------------
		 */

		if (pedID.size() > 10000.0) {

			if (isH) {
				status = trace_back_par (f_ped, pedigree, pedID); // f_ped is original pedigree. pedigree is full traced pedigree
				if (status != 0) throw status;

				runtimelog (false, "size of pedigree (in memory), B", static_cast<double>(mapCapacity(pedigree)), false, false);
				runtimelog (false, "elements in pedigree", static_cast<double>(pedigree.size()), false, false);
			}

			status = trace_back_par (f_ped, r_pedigree, genotypedID); // f_ped is original pedigree. pedigree is full traced pedigree
			if (status != 0) throw status;

			runtimelog (false, "size of reduced pedigree (in memory), B", static_cast<double>(mapCapacity(r_pedigree)), false, false);
			runtimelog (false, "elements in reduced pedigree", static_cast<double>(r_pedigree.size()), false, false);

		}
		else {

			if (isH) {
				status = trace_back (f_ped, pedigree, pedID); // f_ped is original pedigree. pedigree is full traced pedigree
				if (status != 0) throw status;

				runtimelog (false, "size of pedigree (in memory), B", static_cast<double>(mapCapacity(pedigree)), false, false);
				runtimelog (false, "elements in pedigree", static_cast<double>(pedigree.size()), false, false);
			}

			status = trace_back (f_ped, r_pedigree, genotypedID); // f_ped is original pedigree. pedigree is full traced pedigree
			if (status != 0) throw status;

			runtimelog (false, "size of reduced pedigree (in memory), B", static_cast<double>(mapCapacity(r_pedigree)), false, false);
			runtimelog (false, "elements in reduced pedigree", static_cast<double>(r_pedigree.size()), false, false);

		}

		std::vector<int> nongenotID; // make list of non genotyped IDs

		for (auto const& e : r_pedigree) {
			if (! find_invect(genotypedID, e.first.id))
				nongenotID.push_back(e.first.id);
		}
		std::sort(nongenotID.begin(), nongenotID.end());

		runtimelog (false, "elements of non-genotyped IDs", static_cast<double>(nongenotID.size()), false, false);

		f_ped.clear(); // we do not need it any more
		birth_id_map.clear();

		/* Optional part related to full H_inv*/

		size_t ainv_sz;

		if (isH) {

			runtimelog (true, "Build full pedigree A(-1) and write it to a file ...", 0, false, false);

			status = get_ainv (pedigree, ainv, inbreeding); // here ainv consists of non-zero elements associated to [row, col] which are not yet recoded IDs
			if (status != 0) throw status;

			runtimelog (false, "size of ainv (in memory), B", static_cast<double>(mapCapacity(ainv)), false, false);
			runtimelog (false, "elements in ainv", static_cast<double>(ainv.size()), false, false);

			ainv_sz = ainv.size();

			// write out A_inv based on full pedigree
			try {
				strAinv *a; a = (strAinv *)malloc( sizeof( strAinv ));
				if (a == NULL) throw 40;
				for (auto const &e : ainv) {
					a->id_1 = e.first.day;
					a->id_2 = e.first.id;
					a->val = e.second;
					fAfi.write(reinterpret_cast<char*> (a), sizeof( strAinv ));
				}
				free(a);
			}
			catch (int ex) {
				write_log("Memory for a", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception in writing ainv to the file", 3);
				throw 3;
			}

			/*
			//----- Start DEBUG PRINT -----------------------------------------------------------
			std::cout<<"fcode_map:"<<std::endl;
			auto debug_count = 0;
			for(auto const& elem : fcode_map) {
				if (debug_count < 40) { std::cout << "key: "<<elem.first <<" "<< ",  values: "<<elem.second<<std::endl; debug_count++; }
			}
			std::cout<<std::endl;

			std::cout<<"full pedigree:"<<std::endl;
			debug_count = 0;
			for(auto const& elem : pedigree) {
				if (debug_count < 100) { std::cout << "key: "<<elem.first.day <<" "<<elem.first.id<< ",  values: "<<elem.second.id_1<<" "<<elem.second.id_2<<"\n"; debug_count++; }
			}
			std::cout<<std::endl;

			std::cout<<"reduced pedigree:"<<std::endl;
			debug_count = 0;
			for(auto const& elem : r_pedigree) {
				if (debug_count < 40) { std::cout << "key: "<<elem.first.day <<" "<<elem.first.id<< ",  values: "<<elem.second.id_1<<" "<<elem.second.id_2<<"\n"; debug_count++; }
			}
			std::cout<<std::endl;
			debug_count = 0;
			for(auto const& elem : pedID) {
				if (debug_count < 20) { std::cout<<"ped IDs : "<<elem<<std::endl; debug_count++; }
			}
			std::cout<<std::endl;
			debug_count = 0;
			for(auto const& elem : genotypedID) {
				if (debug_count < 20) { std::cout<<"genotyped IDs : "<<elem<<std::endl; debug_count++; }
			}
			std::cout<<std::endl;
			debug_count = 0;
			for(auto const& elem : coreID) {
				if (debug_count < 20) { std::cout<<"core IDs : "<<elem<<std::endl; debug_count++; }
			}
			std::cout<<"full A inv:"<<std::endl;
			debug_count = 0;
			for(auto const& elem : ainv) {
				if (debug_count < 200) { std::cout << "key: "<<elem.first.day <<" "<<elem.first.id<< ",  values: "<<elem.second<<"\n"; debug_count++; }
			}
			AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
			printf("\nprogram uses %li bytes in %d buffers",AllocatedBytes,N_AllocatedBuffers);
			std::cout<<std::endl<<std::endl;
			//----- End DEBUG PRINT -------------------------------------------------------------
			 */

			pedigree.clear();
			ainv.clear();
			fcode_map.clear();

		}

		runtimelog (true, "Build reduced pedigree A(-1) and write it to a file ...", 0, false, false);

		status = get_ainv (r_pedigree, r_ainv, inbreeding); // here ainv consists of non-zero elements associated to [row, col] which are not yet recoded IDs
		if (status != 0) throw status;

		runtimelog (false, "size of r_ainv (in memory), B", static_cast<double>(mapCapacity(r_ainv)), false, false);
		runtimelog (false, "elements in r_ainv", static_cast<double>(r_ainv.size()), false, false);

		size_t r_ainv_size = r_ainv.size();

		// write out A_inv based on reduced pedigree
		try {
			strAinv *a; a = (strAinv *)malloc( sizeof( strAinv ));
			if (a == NULL) throw 40;
			for (auto const &e : r_ainv) {
				a->id_1 = e.first.day;
				a->id_2 = e.first.id;
				a->val = e.second;
				fAri.write(reinterpret_cast<char*> (a), sizeof( strAinv ));
			}
			free(a);
		}
		catch (int ex) {
			write_log("Memory for a", ex);
			throw ex;
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception in writing r_ainv to the file", 3);
			throw 3;
		}

		if (inbreeding) {
			/* record inbreeding coefficients to the file */

			std::ofstream fHtxt;
			fHtxt.open ("inbreeding_coef", fHtxt.out);
			fHtxt.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
			if (!fHtxt.is_open()) throw 3;

			size_t eNums = 0;
			size_t animId = 1;
			for (auto const& e: r_pedigree) {
				fHtxt << e.first.id <<" " << inbrF[animId] << std::endl;
				animId++;
			}
			fHtxt.close();
		}

		/*
		//----- Start DEBUG PRINT -----------------------------------------------------------
		std::cout<<"fcode_map:"<<std::endl;
		auto debug_count = 0;
		for(auto const& elem : fcode_map) {
			if (debug_count < 40) { std::cout << "key: "<<elem.first <<" "<< ",  values: "<<elem.second<<std::endl; debug_count++; }
		}
		std::cout<<std::endl;

		std::cout<<"full pedigree:"<<std::endl;
		debug_count = 0;
		for(auto const& elem : pedigree) {
			if (debug_count < 100) { std::cout << "key: "<<elem.first.day <<" "<<elem.first.id<< ",  values: "<<elem.second.id_1<<" "<<elem.second.id_2<<"\n"; debug_count++; }
		}
		std::cout<<std::endl;

		std::cout<<"reduced pedigree:"<<std::endl;
		debug_count = 0;
		for(auto const& elem : r_pedigree) {
			if (debug_count < 40) { std::cout << "key: "<<elem.first.day <<" "<<elem.first.id<< ",  values: "<<elem.second.id_1<<" "<<elem.second.id_2<<"\n"; debug_count++; }
		}
		std::cout<<std::endl;
		debug_count = 0;
		for(auto const& elem : pedID) {
			if (debug_count < 20) { std::cout<<"ped IDs : "<<elem<<std::endl; debug_count++; }
		}
		std::cout<<std::endl;
		debug_count = 0;
		for(auto const& elem : genotypedID) {
			if (debug_count < 20) { std::cout<<"genotyped IDs : "<<elem<<std::endl; debug_count++; }
		}
		std::cout<<std::endl;
		debug_count = 0;
		for(auto const& elem : coreID) {
			if (debug_count < 20) { std::cout<<"core IDs : "<<elem<<std::endl; debug_count++; }
		}
		std::cout<<"reduced A inv:"<<std::endl;
		debug_count = 0;
		for(auto const& elem : r_ainv) {
			if (debug_count < 200) { std::cout << "key: "<<elem.first.day <<" "<<elem.first.id<< ",  values: "<<elem.second<<"\n"; debug_count++; }
		}
		AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
		printf("\nprogram uses %li bytes in %d buffers",AllocatedBytes,N_AllocatedBuffers);
		std::cout<<std::endl<<std::endl;
		//----- End DEBUG PRINT -------------------------------------------------------------
		 */

		r_ainv.clear();

		runtimelog (true, "Make A22, based on reduced pedigree ...", 0, false, false);

		// Make A22, based on reduced pedigree:

		// 1) Get 'A' matrix:

		MATRIXMATH sstepmat;

		runtimelog (false, "Allocate memory for A22 (half-stored, 'L' part) ...", 0, false, false);

		// Get A22 (half-stored, 'L' part):
		float *Agen; // genotyped
		status = sstepmat.make_res_array_h (Agen, genotypedID.size());
		if (status != 0) {
			write_log ("Memory for A22", status);
			throw status;
		}

		runtimelog (false, "Memory allocated for A22, GiB", (static_cast<double> ( sizeof( float[(genotypedID.size()*genotypedID.size()+genotypedID.size())/2] ) ))/(1024*1024*1024), false, false);

		runtimelog (true, "call get_A22 () ...", 0, false, false);

		/* Matrix A22 (Agen in the code) is aranged according to 'genotypedID' vector */
		status = get_A22 (r_pedigree, Agen);
		if (status != 0) throw status;

		r_pedigree.clear();

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		std::cout<<std::endl;
		printf ("\n directly generated Agen: \n");
		for (auto i=1; i<=min(genotypedID.size(),26); i++) {
			for (auto j=1; j<=i; j++) {
				printf ("%12.5G", Agen[i*(i - 1)/2 + j - 1]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT ---------------------------------------------------
		 */

		runtimelog (true, "Get alpha and betha for A22 ...", 0, false, false);

		// Get mean values
		double a_ofd_mean = 0.0;
		double a_all_mean = 0.0;
		double a_diag_mean = 0.0;

		if (if_adjust) {

#pragma omp parallel for reduction (+: a_diag_mean)
			for (size_t i = 1; i <= genotypedID.size(); i++) {
				a_diag_mean += double(Agen[i*(i - 1)/2 + i - 1]);
			}
			a_all_mean = a_diag_mean;
			a_diag_mean = a_diag_mean/(genotypedID.size());

#pragma omp parallel for reduction (+: a_ofd_mean)
			for (size_t i = 1; i <= genotypedID.size(); i++) {
				for (size_t j = 1; j <= i; j++) {
					if (i != j) {
						a_ofd_mean += double(Agen[i*(i - 1)/2 + j - 1]);
					}
				}
			}
			a_all_mean = (a_all_mean + 2*a_ofd_mean)/(genotypedID.size()*genotypedID.size());
			a_ofd_mean = 2*a_ofd_mean/(genotypedID.size()*genotypedID.size() - genotypedID.size());

		}

		//----- Start DEBUG PRINT -----------------------------------------------------------
		/*
		std::cout<<"a_diag_mean : "<<a_diag_mean<<std::endl;
		std::cout<<"a_all_mean : "<<a_all_mean<<std::endl;
		std::cout<<"a_ofd_mean : "<<a_ofd_mean<<std::endl;
		 */

		runtimelog (false, "program buffers uses, GiB", static_cast<double> (mkl_mem_stat(&N_AllocatedBuffers))/(1024*1024*1024), false, false);

		a.clear();

		runtimelog (true, "Writing A22 to the file ...", 0, false, false);

		try
		{
			fA22.write(reinterpret_cast<char*>(Agen), ((genotypedID.size()*genotypedID.size() + genotypedID.size())/2)*sizeof( float ));
		}
		catch (std::exception const& e)
		{
			write_log (e.what(), 3);
			write_log ("Exception in writing A22 to the file", 3);
			throw 3;
		}

		runtimelog (true, "Completed of writing A22 to the file", 0, false, false);

		mkl_free(Agen);

		/* Produce A22(-1) based on matrrix partition (but not through the direct inversion) */

		/*Working with reduced A_inv :*/

		runtimelog (true, "Making A22(-1). Prepare vectors and restore r_ainv from the file ...", 0, false, false);

		size_t r_ginv = genotypedID.size();
		size_t c_ginv = genotypedID.size();

		// Restore r_ainv
		try {
			dKey k;
			fAri.seekg(0, std::ios::beg);
			strAinv *a; a = (strAinv *)malloc( sizeof( strAinv ));
			if (a == NULL) throw 40;
			for (auto i = 0; i < r_ainv_size; i++) {
				fAri.read(reinterpret_cast<char*> (a), sizeof( strAinv ));
				k.day = a->id_1;
				k.id = a->id_2;
				r_ainv[k] = a->val;
			}
			free(a);
		}
		catch (int ex) {
			write_log("Memory for a", ex);
			throw ex;
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception reading r_ainv from the file", 3);
			throw 3;
		}

		fAri.close(); remove(filename4.c_str());

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		std::cout<<"r_ainv:"<<std::endl;
		for (auto const& elem: r_ainv) {
			std::cout<<elem.first.day<<" "<<elem.first.id<<" "<<elem.second<<std::endl;
		}
		//----- End DEBUG PRINT ---------------------------------------------------
		 */

		runtimelog (true, "Build A11 and write it to a file ...", 0, false, false);

		double *A11; // non-genotyped IDs matrix

		size_t r_a11 = nongenotID.size();
		size_t c_a11 = nongenotID.size();

		/*Give the range of IDs to define the submatrix (USE the ORIGINAL IDs from the pedigree)*/

		//sstepmat.make_res_array (A11, r_a11, c_a11); // full-store matrix
		status = sstepmat.make_res_array_h (A11, r_a11); //half-store matrix
		if (status != 0) {
			write_log("Memory for A11", status);
			throw status;
		}

		//unsigned int n_threads = std::thread::hardware_concurrency();
		unsigned int block_size = static_cast<unsigned int> (nongenotID.size()/(1*n_threads));

		//size_t a11_el = 0;
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
		for (size_t i = 1; i <= nongenotID.size(); i++ ) {
			dKey k;
			k.day = nongenotID[i-1];
			for (size_t j = 1; j <= i; j++ ) {
				k.id = nongenotID[j-1];
				if (r_ainv.count(k)) {
					A11[i*(i - 1)/2 + j - 1] = r_ainv[k];
					//a11_el++;
				}
			}
		}
		//std::cout << "a11_el = " <<a11_el << "\n";

		try {
			fa11.write(reinterpret_cast<char*>(A11), ((r_a11*r_a11 + r_a11)/2)*sizeof( double )); //half-store matrix
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception writing A11 to the file", 3);
			throw 3;
		}

		runtimelog (false, "The size of A11 in memory, GiB", static_cast<double> ( sizeof (double[(r_a11*r_a11 + r_a11)/2]) ) /(1024*1024*1024), false, false);

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		printf ("\n Top left corner of matrix A11: \n");
		for (auto i=1; i<=min(nongenotID.size(),20); i++) {
			for (auto j=1; j<=i; j++) {
				printf ("%12.5G", A11[i*(i - 1)/2 + j - 1]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT ---------------------------------------------------
		 */

		mkl_free(A11);

		runtimelog (true, "Build A22 and write it to a file ...", 0, false, false);

		double *A22; // genotyped

		size_t r_a22 = genotypedID.size();
		size_t c_a22 = genotypedID.size();

		status = sstepmat.make_res_array_h (A22, r_a22); // half-store matrix
		if (status != 0) {
			write_log("Memory for A22", status);
			throw status;
		}

		//n_threads = std::thread::hardware_concurrency();
		block_size = static_cast<unsigned int> (genotypedID.size()/(1*n_threads));

		//size_t a22_el = 0;
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
		for (size_t i = 1; i <= genotypedID.size(); i++ ) {
			dKey k;
			k.day = genotypedID[i-1];
			for (size_t j = 1; j <= i; j++ ) {
				k.id = genotypedID[j-1];
				//A22[(i - 1)*genotypedID.size() + j - 1] = A22[(j - 1)*genotypedID.size() + i - 1] = r_ainv[k]; // full-store matrix
				if (r_ainv.count(k)) {
					A22[i*(i - 1)/2 + j - 1] = r_ainv[k];
					//a22_el++;
				}
			}
		}
		//std::cout << "a22_el = " <<a22_el << "\n";

		try {
			fa22.write(reinterpret_cast<char*>(A22), ((r_a22*r_a22 + r_a22)/2)*sizeof( double ));
			if (fa22.fail()) throw 4;
		}
		catch (int ex) {
			write_log ("Writing A22 to binary file", ex);
			throw ex;
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception writing A22 to the file", 3);
			throw 3;
		}

		runtimelog (false, "The size of A22 in mmemory, GiB", static_cast<double> ( sizeof (double[(r_a22*r_a22 + r_a22)/2]) ) /(1024*1024*1024), false, false);

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		std::cout<<std::endl;
		printf ("\n Top left corner of matrix A22: \n");
		for (auto i=1; i<=min(genotypedID.size(),29); i++) {
			for (auto j=1; j<=i; j++) {
				printf ("%12.5G", A22[i*(i - 1)/2 + j - 1]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT ---------------------------------------------------
		 */

		mkl_free(A22);

		runtimelog (true, "Build A21 and write it to a file ...", 0, false, false);

		double *A21; // genotyped - non-genotyped
		size_t c_a21 = nongenotID.size();
		size_t r_a21 = genotypedID.size();

		status = sstepmat.make_res_array (A21, r_a21, c_a21);
		if (status != 0) {
			write_log("Memory for A21", status);
			throw status;
		}

		//n_threads = std::thread::hardware_concurrency();
		block_size = static_cast<unsigned int> (genotypedID.size()/(1*n_threads));

		//size_t elm_inA21 = 0;
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
		for (size_t i = 1; i <= genotypedID.size(); i++ ) {
			dKey k;
			k.day = genotypedID[i-1];
			bool changed_key = false;
			for (size_t j = 1; j <= nongenotID.size(); j++ ) {
				k.id = nongenotID[j-1];
				/* if k.day < k.id -> the key does not exists!
				 * Hence we will have zero element in a full matrix.*/
				if (k.day < k.id) {
					size_t t_id = k.id;
					k.id = k.day;
					k.day = t_id;
					changed_key = true;
				}
				if (r_ainv.count(k)) {
					A21[(i - 1)*nongenotID.size() + j - 1] = r_ainv[k];
					//elm_inA21++;
				}
				if (changed_key){
					k.day = genotypedID[i-1];
					changed_key = false;
				}
			}
		}

		//std::cout<<"elements in A21 = "<< elm_inA21<<std::endl;

		r_ainv.clear();

		try {
			fa21.write(reinterpret_cast<char*>(A21), r_a21*c_a21*sizeof( double ));
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception writing A21 to the file", 3);
			throw 3;
		}

		runtimelog (false, "The size of A21 in mmemory, GiB", static_cast<double> ( sizeof (double[(r_a21*c_a21)]) ) /(1024*1024*1024), false, false);

		runtimelog (true, "Transpose A21 to get A12 ...", 0, false, false);

		// transpose A21 to get A12
		double *A12;

		status = sstepmat.make_res_array (A12, c_a21, r_a21);
		if (status != 0) {
			write_log("Memory for A12", status);
			throw status;
		}

		//n_threads = std::thread::hardware_concurrency();
		block_size = static_cast<unsigned int> (r_a21/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
		for (size_t i = 0; i < r_a21; i++) {
			for (size_t j = 0; j < c_a21; j++) {
				A12[j*r_a21+i] = A21[i*c_a21+j];
			}
		}

		try {
			fa12.write(reinterpret_cast<char*>(A12), r_a21*c_a21*sizeof( double ));
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception writing A12 to the file", 3);
			throw 3;
		}

		size_t r_a12 = c_a21;
		size_t c_a12 = r_a21;

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		std::cout<<std::endl;
		printf ("\n A21: \n");
		for (auto i=0; i<min(r_a21,19); i++) {
			for (auto j=0; j<min(c_a21,19); j++) {
				printf ("%12.5G", A21[j+i*c_a21]);
			}
			printf ("\n");
		}
		std::cout<<std::endl;
		printf ("\n A12: \n");
		for (auto i=0; i<min(c_a21,26); i++) {
			for (auto j=0; j<min(r_a21,20); j++) {
				printf ("%12.5G", A12[j+i*r_a21]);
			}
			printf ("\n");
		}
		std::cout<<std::endl;

		double ma_el = 0.0;
		double mi_el = 0.0;
		size_t ma_num = 0;
		size_t mi_num = 0;
		for (size_t i = 0; i<r_a21*c_a21; i++){
			if (A12[i] != 0.0){
				if(A12[i] < mi_el){
					mi_num++;
					mi_el = A12[i];
				}
				if(A12[i] > ma_el){
					ma_num++;
					ma_el = A12[i];
				}
			}
		}
		std::cout<<"mi_num, mi_el:"<<mi_num<<", "<<mi_el<<std::endl;
		std::cout<<"ma_num, ma_el:"<<ma_num<<", "<<ma_el<<std::endl;
		std::cout<<std::endl;

		//----- End DEBUG PRINT ---------------------------------------------------
		 */

		mkl_free(A21);
		mkl_free(A12);

		runtimelog (true, "Restore A11 and inverting it ...", 0, false, false);

		try {
			status = sstepmat.make_res_array_h (A11, r_a11);
			if (status != 0) throw status;
			fa11.seekg(0, std::ios::beg);
			fa11.read(reinterpret_cast<char*>(A11), ((r_a11*r_a11 + r_a11)/2)*sizeof( double )); // half-store
		}
		catch (int ex)
		{
			write_log("Memory for A11 when restore it", ex);
			throw ex;
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception reading A11 from the file", 3);
			throw 3;
		}

		fa11.close(); remove(filename11.c_str());

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		printf ("\n A11: \n");
		for (auto i=1; i<=min(nongenotID.size(),29); i++) {
			for (auto j=1; j<=i; j++) {
				printf ("%12.5G", A11[i*(i - 1)/2 + j - 1]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT ---------------------------------------------------
		 */

		double *A11f;

		if (save) {

			runtimelog (true, "SAVE option is sat to TRUE !..", 0, false, false);

			runtimelog (true, "Do halve-stored A11 inverse ...", 0, false, false);

			status = sstepmat.matr_inv_h(A11, r_a11);
			if (status != 0) {
				write_log("Inversion of A11", status);
				throw status;
			}

			/*
		//----- Start DEBUG PRINT -------------------------------------------------
		printf ("\n inverse of matrix A11: \n");
		for (auto i=1; i<=min(nongenotID.size(),20); i++) {
			for (auto j=1; j<=i; j++) {
				printf ("%12.5G", A11[i*(i - 1)/2 + j - 1]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT ---------------------------------------------------
			 */

			// A_11t_12 = A11_inv * A12. Rows = r_a11; columns = c_a12.

			runtimelog (true, "Make A11_i full-stored to do multiplication further on ...", 0, false, false);

			// make A11_inv full to do multiplication further on
			//double *A11f;
			status = sstepmat.make_res_array (A11f, r_a11, c_a11);
			if (status != 0) {
				write_log("Memory for A11f", status);
				throw status;
			}

#pragma omp parallel for
			for (size_t i = 1; i <= nongenotID.size(); i++) {
				for (size_t j = 1; j <= i; j++) {
					A11f[(i - 1)*nongenotID.size() + j - 1] = A11f[(j - 1)*nongenotID.size() + i - 1] = A11[i*(i - 1)/2 + j - 1];
				}
			}

			mkl_free(A11); // we do not need A11 anymore

			runtimelog (false, "The size of A11_full in memory, GiB", static_cast<double> ( sizeof (double[(r_a11*c_a11)]) ) /(1024*1024*1024), false, false);

			/*
		//----- Start DEBUG PRINT -------------------------------------------------
		printf ("\n Top left corner of A11f: \n");
		for (auto i=0; i<min(r_a11,19); i++) {
			for (auto j=0; j<min(c_a11,19); j++) {
				printf ("%12.5G", A11f[j+i*c_a11]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT ---------------------------------------------------
			 */

		}
		else {

			runtimelog (true, "SAVE option is sat to FALSE !..", 0, false, false);

			runtimelog (true, "Do full-stored A11 inverse ...", 0, false, false);

			runtimelog (true, "Make A11 full-stored to do full-stored inversion and multiplication further on ...", 0, false, false);

			// make A11 full to do full-stored inversion and multiplication further on
			//double *A11f;
			status = sstepmat.make_res_array (A11f, r_a11, c_a11);
			if (status != 0) {
				write_log("Memory for A11f", status);
				throw status;
			}

#pragma omp parallel for
			for (size_t i = 1; i <= nongenotID.size(); i++) {
				for (size_t j = 1; j <= i; j++) {
					A11f[(i - 1)*nongenotID.size() + j - 1] = A11f[(j - 1)*nongenotID.size() + i - 1] = A11[i*(i - 1)/2 + j - 1];
				}
			}

			mkl_free(A11); // we do not need A11 anymore

			runtimelog (false, "The size of A11_full in memory, GiB", static_cast<double> ( sizeof (double[(r_a11*c_a11)]) ) /(1024*1024*1024), false, false);

			/* Do A11f inversion*/
			runtimelog (true, "Inverting A11 ...", 0, false, false);

			lapack_int *pvt;
			status = sstepmat.make_pivot_array(pvt, r_a11);
			if (status != 0) {
				write_log("Memory for pvt (A11 inv part)", status);
				throw status;
			}

			status = sstepmat.matr_inv(A11f, pvt, r_a11, r_a11);
			if (status != 0) {
				write_log("Inversion of A11f (inv part)", status);
				throw status;
			}

			mkl_free(pvt);

		}

		runtimelog (true, "A11(-1) * A12 ...", 0, false, false);

		double *A_11t_12;

		status = sstepmat.make_res_array (A_11t_12, r_a11, c_a12);
		if (status != 0) {
			write_log("Memory for A_11t_12", status);
			throw status;
		}

		// restore A12
		try {
			status = sstepmat.make_res_array (A12, c_a21, r_a21);
			if (status != 0) throw status;
			fa12.seekg(0, std::ios::beg);
			fa12.read(reinterpret_cast<char*>(A12), r_a21*c_a21*sizeof( double ));
		}
		catch (int ex) {
			write_log("Memory for A12", ex);
			throw ex;
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception reading A12 from the file", 3);
			throw 3;
		}

		fa12.close(); remove(filename15.c_str());

		sstepmat.matr_prod (A11f, A12, A_11t_12, r_a11, r_a12, c_a11, c_a12);

		mkl_free(A11f); // we do not need A11 anymore
		mkl_free(A12);

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		printf ("\n A11_i * A12: \n");
		for (auto i=0; i<min(r_a11,20); i++) {
			for (auto j=0; j<min(c_a12,20); j++) {
				printf ("%12.5G", A_11t_12[j+i*c_a12]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT ---------------------------------------------------
		 */

		runtimelog (true, "A21 * A11(-1)_A12 ...", 0, false, false);
		runtimelog (false, "Restore A21 from binary ...", 0, false, false);

		try {
			status = sstepmat.make_res_array (A21, r_a21, c_a21);
			if (status != 0) throw status;
			fa21.seekg(0, std::ios::beg);
			fa21.read(reinterpret_cast<char*>(A21), r_a21*c_a21*sizeof( double ));
		}
		catch (int ex) {
			write_log("Memory for A21", ex);
			throw ex;
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception reading A21 from the file", 3);
			throw 3;
		}

		fa21.close(); remove(filename12.c_str());

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		std::cout<<std::endl;
		printf ("\n restored A21: \n");
		for (auto i=0; i<min(r_a21,20); i++) {
			for (auto j=0; j<min(c_a21,20); j++) {
				printf ("%12.5G", A21[j+i*c_a21]);
			}
			printf ("\n");
		}
		std::cout<<std::endl;
		//----- End DEBUG PRINT ---------------------------------------------------
		 */

		// A21_11t_12 = A21 * A_11t_12. Rows = columns = c_a12.
		double *A21_11t_12;
		status = sstepmat.make_res_array (A21_11t_12, c_a12, c_a12);
		if (status != 0) {
			write_log("Memory for A21_11t_12", status);
			throw status;
		}

		runtimelog (false, "Do matrix product ...", 0, false, false);

		sstepmat.matr_prod (A21, A_11t_12, A21_11t_12, r_a21, r_a11, c_a21, c_a12);

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		printf ("\n Top left corner of A21 * A11_i_A12: \n");
		for (auto i=0; i<min(c_a12,20); i++) {
			for (auto j=0; j<min(c_a12,20); j++) {
				printf ("%12.5G", A21_11t_12[j+i*c_a12]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT --------------------------------------------------
		 */

		mkl_free(A_11t_12); // we do not need A_11t_12 anymore
		mkl_free(A21); // we do not need A21 and A12 anymore

		runtimelog (true, "A22 - A21_A11(-1)_A12 ...", 0, false, false);
		runtimelog (false, "Restore A22 from binary ...", 0, false, false);

		try {
			status = sstepmat.make_res_array_h (A22, r_a22);
			if (status != 0) throw status;
			fa22.seekg(0, std::ios::beg);
			fa22.read(reinterpret_cast<char*>(A22), ((r_a22*r_a22 + r_a22)/2)*sizeof( double ));
			if (fa22.fail()) throw 4;
		}
		catch (int ex) {
			write_log("Memory for A22", ex);
			throw ex;
		}
		catch (std::exception const& e) {
			write_log (e.what(), 3);
			write_log ("Exception reading A22 from the file", 3);
			throw 3;
		}

		fa22.close(); remove(filename13.c_str());

		runtimelog (false, "Do matrix substitution ...", 0, false, false);

		// Get A22inv: A22 = A22 - A21_11t_12 (only 'L' part);
		for (size_t i = 1; i <= genotypedID.size(); i++) {
			for (size_t j = 1; j<= i; j++) {
				A22[i*(i-1)/2 + j-1] = A22[i*(i-1)/2 + j-1] - A21_11t_12[(i - 1)*genotypedID.size() + j - 1];
			}
		}

		mkl_free(A21_11t_12); // we do not need A21_11t_12 anymore

		/*
		//----- Start DEBUG PRINT -------------------------------------------------
		printf ("\n Top left corner of A22inv: \n");
		for (auto i=1; i<=min(genotypedID.size(),20); i++) {
			for (auto j=1; j<=i; j++) {
				printf ("%12.5G", A22[i*(i - 1)/2 + j - 1]);
			}
			printf ("\n");
		}
		//----- End DEBUG PRINT --------------------------------------------------
		 */

		/*
		{

			 for debugging: A22_inv * A22 = I

			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<"DEBUG: Restore A22 from the binary file ..."<<std::endl<<std::endl;
			//----- End DEBUG PRINT ---------------------------------------------------

			double *Agen;
			status = sstepmat.make_res_array_h (Agen, genotypedID.size());
			if (status != 0) {
				write_log("Memory for Agen", status);
				throw status;
			}

			size_t sz = (genotypedID.size()*genotypedID.size() + genotypedID.size())/2;
			try {
				float *a; a = (float *)malloc( sizeof( float ));
				if (a == NULL) throw 40;
				fA22.seekg(0, std::ios::beg);
				for (auto i = 0; i < sz; i++) {
					a[0] = 0.0f;
					fA22.read(reinterpret_cast<char*> (a), sizeof( float ));
					Agen[i] = static_cast<double>(a[0]);
				}
				free(a);
			}
			catch (int ex) {
				write_log("Memory for a", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading file of A22 (s.step part) matrix", 3);
				throw 3;
			}

			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<std::endl;
			printf ("\n DEBUG: directly generated Agen: \n");
			for (auto i=1; i<=min(genotypedID.size(),20); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", Agen[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}

			double ma_el = 0.0;
			double mi_el = 0.0;
			size_t ma_num = 0;
			size_t mi_num = 0;
			for (size_t i = 0; i<genotypedID.size(); i++){
				for (size_t j = 0; j<=i; j++){
					if (i != j){
						if (Agen[i*(i + 1)/2 + j] != 0.0){
							if(Agen[i*(i + 1)/2 + j] < mi_el){
								mi_num++;
								mi_el = Agen[i*(i + 1)/2 + j];
							}
							if(Agen[i*(i + 1)/2 + j] > ma_el){
								ma_num++;
								ma_el = Agen[i*(i + 1)/2 + j];
							}
						}
					}
				}
			}
			std::cout<<"mi_num, mi_el:"<<mi_num<<", "<<mi_el<<std::endl;
			std::cout<<"ma_num, ma_el:"<<ma_num<<", "<<ma_el<<std::endl;
			std::cout<<std::endl;

			//----- End DEBUG PRINT ---------------------------------------------------

			// make Agen full to do multiplication further on
			double *Agenf;
			status = sstepmat.make_res_array (Agenf, genotypedID.size(), genotypedID.size());
			if (status != 0) {
				write_log("Memory for Agenf", status);
				throw status;
			}

			double *A22f;
			status = sstepmat.make_res_array (A22f, genotypedID.size(), genotypedID.size());
			if (status != 0) {
				write_log("Memory for A22f", status);
				throw status;
			}

//#pragma omp parallel for
			for (size_t i = 0; i < genotypedID.size(); i++) {
				for (size_t j = 0; j <= i; j++) {
					Agenf[(i)*genotypedID.size() + j] = Agenf[(j)*genotypedID.size() + i] = Agen[i*(i + 1)/2 + j];
					A22f[(i)*genotypedID.size() + j] = A22f[(j)*genotypedID.size() + i] = A22[i*(i + 1)/2 + j];
				}
			}

			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of Agenf: \n");
			for (auto i=0; i<min(genotypedID.size(),20); i++) {
				for (auto j=0; j<min(genotypedID.size(),20); j++) {
					printf ("%12.5G", Agenf[i*genotypedID.size() + j]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT --------------------------------------------------

			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of A22f: \n");
			for (auto i=0; i<min(genotypedID.size(),20); i++) {
				for (auto j=0; j<min(genotypedID.size(),20); j++) {
					printf ("%12.5G", A22f[i*genotypedID.size() + j]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT --------------------------------------------------

			//mkl_free(Agen); // we do not need A11 anymore

			double *I;
			status = sstepmat.make_res_array (I, genotypedID.size(), genotypedID.size());
			if (status != 0) {
				write_log("Memory for I", status);
				throw status;
			}

			size_t dim = genotypedID.size();
			//cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, 1.0, Agenf, dim, A22f, dim, 0.0, I, dim);
			sstepmat.matr_prod(A22f, Agenf, I, genotypedID.size(), genotypedID.size(), genotypedID.size(), genotypedID.size());

			double big = 0.00000001;
			double small = -0.0000001;

			size_t ones = 0;
			size_t mins = 0;
			double min_e = small;
			double max_e = big;
			size_t maxes = 0;
			for (size_t i = 0; i < genotypedID.size()*genotypedID.size(); i++) {
				if (I[i] != 0.0){
					if (I[i] >= 1.0+small && I[i] <= 1.0+big) {
						ones++;
					}
					else if(I[i] < small && I[i] < min_e){
						mins++;
						min_e = I[i];
					}
					else if(I[i] > big && I[i] > max_e){
						max_e = I[i];
					}
					else if (I[i] > big) maxes++;
				}
			}

			std::cout<<"ones, genotyped: "<<ones<<", "<<genotypedID.size()<<std::endl;
			std::cout<<"mins: "<<mins<<std::endl;
			std::cout<<"maxes : "<<maxes<<std::endl;
			if (mins != 0) std::cout<<"min_e: "<<min_e<<std::endl;
			if (maxes != 0) std::cout<<"max_e: "<<max_e<<std::endl;

			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of I: \n");
			for (auto i=0; i<min(genotypedID.size(),20); i++) {
				for (auto j=0; j<=i; j++) {
					printf ("%12.5G", I[i*genotypedID.size() + j]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT --------------------------------------------------

		}
		 */

		runtimelog (true, "Write out A22(-1) to a file ...", 0, false, false);

		double *aVal;
		status = sstepmat.make_res_array(aVal, 1, 1);
		if (status != 0) {
			write_log("Memory for aVal", status);
			throw status;
		}

		/*write out inversed A22*/
		size_t wrAllElemA = 0;
		size_t wrZerElemA = 0;

		for (size_t i = 0; i < (r_a22*r_a22 + r_a22)/2; i++) {
			try {
				wrAllElemA++;

				if (std::abs(A22[i]) > tol_a22inv)
					aVal[0] = A22[i];
				else {
					aVal[0] = 0.0;
					wrZerElemA++;
				}

				fa22i.write(reinterpret_cast<char*>(aVal), sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception writing A22i to the file", 3);
				throw 3;
			}
		}

		runtimelog (false, "Elements in A22i turned to zero, % ", static_cast<double>(100*wrZerElemA/wrAllElemA), false, false);
		runtimelog (false, "Zero threshold applied to A22i ", tol_a22inv, false, false);

		mkl_free(aVal);
		mkl_free(A22); // we do not need A22 anymore

		/*
		 //=============== Single step block: direct inversion of A22 =============================================
		if (coreID.size() == 0 || (coreID.size() == genotypedID.size())) {

			runtimelog (true, "Enter the Single step block", 0, false, false);

			double *Ageni;

			status = sstepmat.make_res_array_h(Ageni, genotypedID.size()); // 'L' part
			if (status != 0) {
				write_log ("Memory for A22i", status);
				throw status;
			}

			runtimelog (true, "Restore A22 from the binary file ...", 0, false, false);

			size_t sz = (genotypedID.size()*genotypedID.size() + genotypedID.size())/2;
			try {
				float *a; a = (float *)malloc( sizeof( float ));
				if (a == NULL) throw 40;
				fA22.seekg(0, std::ios::beg);
				for (size_t i = 0; i < sz; i++) {
					a[0] = 0.0f;
					fA22.read(reinterpret_cast<char*> (a), sizeof( float ));
					Ageni[i] = static_cast<double>(a[0]);
				}
				free(a);
			}
			catch (int ex) {
				write_log("Memory for a", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading file of A22 (s.step part) matrix", 3);
				throw 3;
			}

			runtimelog (true, "Inversion of A22 ...", 0, false, false);

			status = sstepmat.matr_inv_h(Ageni, genotypedID.size());
			if (status != 0) {
				write_log ("Exception in A22 inversion (s.step part) matrix", status);
				throw status;
			}

			runtimelog (true, "Writing A22_inv to the binary file ...", 0, false, false);

			try {
				fA22i.write(reinterpret_cast<char*>(Ageni), sz*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception writing inverted A22 (s.step part) to the file", 3);
				throw 3;
			}


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of inverse of A22: \n");
			for (auto i=1; i<=min(genotypedID.size(),26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", Ageni[i*(i - 1)/2 + j - 1]);
					//printf ("%12.5G", Ageni[j+i*genotypedID.size()]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------


			mkl_free(Ageni); // free memory since we have A22inv in file

			runtimelog (true, "Exit the Single step block", 0, false, false);

		}
		// =============== End of Single step block ==============================================================
		 */

		mkl_free_buffers();

		std::vector<int> g_r;
		std::vector<int> g_c;
		std::vector<float> g_v;

		/* If the marker data file mentioned in the parameters file - produce G matrix, otherwise read existing G matrix file */
		//std::string dataFile = apyDatStr.file_marker;

		if( !marker_file.empty() ) {
			
			runtimelog (true, "Producing G matrix using marker file ...", 0, false, false);

			GMAT g(marker_file);
			g.makeG();
			//g.G.print("G matrix from APY");
			g.getValues(g_r, g_c, g_v, gdiag_val);
			g.getIDs(gmatID);
			
		}
		else {

			runtimelog (true, "Read G matrix from file ...", 0, false, false);

			tmpSatatus = true;
			theLock = 0;
			modifFlag = 0;


	#pragma omp parallel sections num_threads (2)
			{
	#pragma omp section
				{
					status = read_genmT (g_r, g_c, g_v);
					tmpSatatus = false;
					if (status != 0) throw status;
				}

	#pragma omp section
				{
					while (tmpSatatus)
						writeOut (false, "Remains to read, % ", 1.0, false, false);

				}
			}

		}

		runtimelog (false, "number of unique gmatIDs", static_cast<double>(gmatID.size()), false, false);

		runtimelog (true, "Check IDs ...", 0, false, false);

		status = check_ID();
		if (status != 0) throw status;

		runtimelog (true, "Get averages of G ...", 0, false, false);

		// Get mean values
		double g_ofd_mean = 0.0;
		double g_all_mean = 0.0;
		double g_diag_mean = 0.0;

		double alpha = 0.0;
		double betha = 1.0;

		if (if_adjust) {

#pragma omp parallel for reduction (+: g_diag_mean, g_ofd_mean)
			for(size_t i = 0; i < g_r.size(); i++) {
				if (g_r[i] == g_c[i]) g_diag_mean = g_diag_mean + static_cast<double> (g_v[i]);
				else g_ofd_mean = g_ofd_mean + static_cast<double> (g_v[i]);
			}
			g_all_mean = (g_diag_mean + 2*g_ofd_mean)/(genotypedID.size()*genotypedID.size());
			g_diag_mean = g_diag_mean/genotypedID.size();

			betha = (a_all_mean - a_diag_mean) / (g_all_mean - g_diag_mean);
			alpha = a_diag_mean - g_diag_mean * betha;

			g_ofd_mean = 2*g_ofd_mean/(genotypedID.size()*genotypedID.size()-genotypedID.size());
		}

		//----- Start DEBUG PRINT -----------------------------------------------------------
		/*
		std::cout<<std::endl;
		std::cout<<"g_diag_mean : "<<g_diag_mean<<std::endl;
		std::cout<<"g_all_mean : "<<g_all_mean<<std::endl;
		std::cout<<"g_ofd_mean : "<<2*g_ofd_mean/(genotypedID.size()*genotypedID.size()-genotypedID.size())<<std::endl;
		std::cout<<"alpha : "<<alpha<<std::endl;
		std::cout<<"betha : "<<betha<<std::endl;
		std::cout<<"adjust, w : "<<if_adjust<<", "<<wa<<std::endl;
		if (if_adjust) {
			double betha2 = (g_diag_mean - g_ofd_mean)/(a_diag_mean - a_ofd_mean);
			double alpha2 = a_diag_mean*betha2 - g_diag_mean;
			std::cout<<"alpha2 : "<<alpha2<<std::endl;
			std::cout<<"betha2 : "<<betha2<<std::endl;

		}

		std::cout<<std::endl;
		std::cout<<"G matrix : "<<std::endl;
		for (auto i = 0; i < g_v.size(); i++) {
			std::cout<<g_r[i]<<", "<<g_c[i]<<", "<<g_v[i]<<std::endl;
		}
		 */

		double betha2, alpha2;
		if (if_adjust) {
			betha2 = (g_diag_mean - g_ofd_mean)/(a_diag_mean - a_ofd_mean);
			alpha2 = a_diag_mean*betha2 - g_diag_mean;
		}

		runtimelog (false, "program's buffers uses, GiB", static_cast<double> (mkl_mem_stat(&N_AllocatedBuffers))/(1024*1024*1024), false, false);
		runtimelog (false, "number of genotyped IDs", static_cast<double>(genotypedID.size()), false, false);
		runtimelog (false, "number of values in G", static_cast<double>(g_v.size()), false, false);
		runtimelog (false, "number of rows in G", static_cast<double>(g_r.size()), false, false);
		runtimelog (false, "number of cols in G", static_cast<double>(g_c.size()), false, false);
		runtimelog (false, "The size of G in memory, GiB", static_cast<double> (sizeof (int[g_r.size()]) + sizeof (int[g_c.size()]) + sizeof (float[g_v.size()])) /(1024*1024*1024), false, false);

		runtimelog (true, "Statistics of G & A elements ...", 0, false, false);
		runtimelog (false, "G diagonals mean", static_cast<double>(g_diag_mean), false, false);
		runtimelog (false, "G off-diagonals mean", static_cast<double>(g_ofd_mean), false, false);
		runtimelog (false, "G all elements mean", static_cast<double>(g_all_mean), false, false);
		runtimelog (false, "A diagonals mean", static_cast<double>(a_diag_mean), false, false);
		runtimelog (false, "A off-diagonals mean", static_cast<double>(a_ofd_mean), false, false);
		runtimelog (false, "A all elements mean", static_cast<double>(a_all_mean), false, false);
		runtimelog (false, "alpha", static_cast<double>(alpha), false, false);
		runtimelog (false, "betha", static_cast<double>(betha), false, false);
		runtimelog (false, "alpha2", static_cast<double>(alpha2), false, false);
		runtimelog (false, "betha2", static_cast<double>(betha2), false, false);


		runtimelog (true, "Get scaling of G (in vector form): g_v[i] = g_v[i] * betha + alpha ...", 0, false, false);

#pragma omp parallel for
		for(size_t i = 0; i < g_r.size(); i++) {
			g_v[i] = static_cast<float> (g_v[i] * betha + alpha);
		}

		runtimelog (false, "program's buffers uses, GiB", static_cast<double> (mkl_mem_stat(&N_AllocatedBuffers))/(1024*1024*1024), false, false);


		/*=============== SINGLE STEP BLOCK:  =====================================================================*/

		/* If single-step, get G(-1) (through the direct inversion), then subtract A22(-1), then produce H(-1) */

		if (coreID.size() == 0 || (coreID.size() == genotypedID.size())) {

			runtimelog (true, "Enter the Single step block", 0, false, false);
			runtimelog (true, "Build the L-part of G; copy vector g_v[i] ...", 0, false, false);

			double *G;

			status = sstepmat.make_res_array_h(G, genotypedID.size());
			if (status != 0) {
				write_log("Memory for G matrix (s.step part)", status);
				throw status;
			}

			std::map<int, int> rid_map;

			status = getRecodedIdMap(rid_map, genotypedID);
			if (status != 0) throw status;

			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (g_r.size()/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for(size_t i = 0; i < g_r.size(); i++) {
				size_t r = rid_map[g_r[i]];
				size_t c = rid_map[g_c[i]];
				size_t ind = r*(r - 1)/2 + c -1;
				if (c > r) ind = c*(c - 1)/2 + r -1;
				G[ind] = static_cast<double> (g_v[i]);
			}

			g_r.clear();
			g_r.shrink_to_fit();
			g_c.clear();
			g_c.shrink_to_fit();
			g_v.clear();
			g_v.shrink_to_fit();
			rid_map.clear();

			// Complete scaling of G matrix
			// This is much slower way since we read A22 element-by-element directly from disc, which is not in parallel
			// This is just in order to save RAM

			try {
				fA22.seekg(0, std::ios::beg);
				float *a; a = (float *)malloc( 1*sizeof( float ));
				if (a == NULL) throw 40;
				size_t lda = genotypedID.size();
				for (size_t i = 0; i < (lda*lda + lda)/2; i++) {
					fA22.read(reinterpret_cast<char*> (a), 1*sizeof( float ));
					if (G[i] == 0.0) G[i] = alpha;
					G[i] = static_cast<double> ((1-wa) * G[i]  + wa * a[0]);
				}
				free(a);
			}
			catch (int ex) {
				write_log("Memory for a", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Reading A22 (s.step part)", 3);
				throw 3;
			}

			fA22.close(); remove(filename1.c_str());

			/*
			std::vector<int> codeList;
			getList_genID_recod (genotypedID, g_r, g_c, codeList);

#pragma omp parallel for
			for(auto i = 0; i < g_r.size(); i++) {
				G[codeList[i]] = static_cast<double> (g_v[i]);
			}

			g_r.clear();
			g_r.shrink_to_fit();
			g_c.clear();
			g_c.shrink_to_fit();
			g_v.clear();
			g_v.shrink_to_fit();
			codeList.clear();
			codeList.shrink_to_fit();

			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<"Restore A22 from the file:"<<std::endl<<std::endl;
			//----- End DEBUG PRINT ---------------------------------------------------

			// Restore A22 from the file
			sstepmat.make_res_array_h (Agen, genotypedID.size());

			try {
				fA22.seekg(0, std::ios::beg);
				size_t sz = (genotypedID.size()*genotypedID.size() + genotypedID.size())/2;
				fA22.read(reinterpret_cast<char*>(Agen), sz*sizeof( float ));
			}
			catch (std::exception const& e) {
				write_log (where, 10);
				std::cerr << "Exception reading file associated with fA22: "<< e.what() << std::endl;
				return 10;
			}

			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<"Complete scaling of G:"<<std::endl<<std::endl;
			//----- End DEBUG PRINT ---------------------------------------------------

			// complete scaling
#pragma omp parallel for
			for(auto i = 1; i <= genotypedID.size(); i++) {
				for(auto j = 1; j <= i; j++) {
					size_t ind = i*(i - 1)/2 + j - 1;
					double g = G[ind];
					if (g == 0.0) g = alpha;
					G[ind] = (1-wa) * g  + wa * Agen[ind];
				}
			}

			mkl_free(Agen);
			 */


			/*
			//----- Start DEBUG PRINT -----------------------------------------------------------
			printf ("\n complete G: \n");
			for (auto i=1; i<=min(genotypedID.size(),26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", G[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT -------------------------------------------------------------
			 */

			if (!save) {

				runtimelog (true, "SAVE option is sat to FALSE !..", 0, false, false);

				runtimelog (true, "Do full-stored G inverse ...", 0, false, false);

				runtimelog (true, "Make G full-stored ...", 0, false, false);

				// make G full-store to do not-memory save inversion further on
				double *Gfull;
				status = sstepmat.make_res_array (Gfull, genotypedID.size(), genotypedID.size());
				if (status != 0) {
					write_log("Memory for Gfull", status);
					throw status;
				}

#pragma omp parallel for
				for (size_t i = 1; i <= genotypedID.size(); i++) {
					for (size_t j = 1; j <= i; j++) {
						Gfull[(i - 1)*genotypedID.size() + j - 1] = Gfull[(j - 1)*genotypedID.size() + i - 1] = G[i*(i - 1)/2 + j - 1];
					}
				}

				runtimelog (true, "Inverting Gfull ...", 0, false, false);

				lapack_int *pvt;
				status = sstepmat.make_pivot_array(pvt, genotypedID.size());
				if (status != 0) {
					write_log("Memory for pvt (s.step part)", status);
					throw status;
				}

				status = sstepmat.matr_inv(Gfull, pvt, genotypedID.size(), genotypedID.size());
				if (status != 0) {
					write_log("Inversion of Gfull (s.step part)", status);
					throw status;
				}

				mkl_free(pvt);

				runtimelog (true, "Return to the halve-stored format for G ...", 0, false, false);

				/* return to the halve-stored format for G */
#pragma omp parallel for
				for (size_t i = 1; i <= genotypedID.size(); i++) {
					for (size_t j = 1; j <= i; j++) {
						G[i*(i - 1)/2 + j - 1] = Gfull[(i - 1)*genotypedID.size() + j - 1];
					}
				}

				mkl_free(Gfull);

			}
			else {

				/* memory-save inversion of G*/
				runtimelog (true, "SAVE option is sat to TRUE !..", 0, false, false);

				runtimelog (true, "Do halve-stored G inverse ...", 0, false, false);

				runtimelog (true, "Inverting G ...", 0, false, false);

				status = sstepmat.matr_inv_h(G, genotypedID.size());
				if (status != 0) {
					write_log("Inversion of G (s.step part)", status);
					throw status;
				}

			}


			/*
			//----- Start DEBUG PRINT -----------------------------------------------------------
			std::cout<<std::endl;
			printf ("\n inverse of G: \n");
			for (auto i=1; i<=min(genotypedID.size(),26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", G[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			std::cout<<std::endl;
			AllocatedBytes = mkl_mem_stat(&N_AllocatedBuffers);
			printf("\nprogram uses %li bytes in %d buffers",AllocatedBytes,N_AllocatedBuffers);
			std::cout<<std::endl<<std::endl;
			//----- End DEBUG PRINT -------------------------------------------------------------
			 */

			runtimelog (true, "Making G - A22 ...", 0, false, false);

			double *A22i; // only one cell of double is needed since we shall not store all A22(-1) in memory
			status = sstepmat.make_res_array (A22i, 1, 1);
			if (status != 0) {
				write_log("Memory for A22i (s.step part)", status);
				throw status;
			}

			/* Variables to calculate sparsity of Ginv and Hsub */

			size_t allVals = 0;
			size_t zerrG = 0;
			size_t zerrH = 0;
			double tmpG, tmpGA;

			try {
				//fA22i.seekg(0, std::ios::beg);
				fa22i.seekg(0, std::ios::beg);
				size_t sz = (genotypedID.size()*genotypedID.size()+genotypedID.size())/2; // half-store case
				for (size_t i = 0; i < sz; i++) {
					//fA22i.read(reinterpret_cast<char*>(A22i), sizeof( double )); /* direct A22 inversion */
					fa22i.read(reinterpret_cast<char*>(A22i), sizeof( double )); /* partitioned A22 inversion */

					tmpG = *(G+i);

					*(G+i) = *(G+i)-*(A22i);

					tmpGA = *(G+i);

					allVals++;
					if ( tmpG == 0.0 ) zerrG++;
					if ( tmpGA == 0.0 ) zerrH++;
				}
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Read A22i (s.step part)", 3);
				throw 3;
			}


			runtimelog (false, "Elements in G_inv which are zero, % ", static_cast<double>(100*zerrG/allVals), false, false);
			runtimelog (false, "Elements in [G - A22] which are zero, % ", static_cast<double>(100*zerrH/allVals), false, false);

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of [G(-1) - A22(-1)]: \n");
			for (auto i=1; i<=min(genotypedID.size(),26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", G[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			runtimelog (true, "Writing to a file the part G-A22 ...", 0, false, false);

			if (isHtxt) {
				status = write_toFile_h ("H_sub.dat", G, genotypedID, false); // for half-stored format; TXT file
				if (status != 0) throw status;
			}

			status = write_toFile_h ("H_sub.bin", G, genotypedID, true); // for half-stored format; BIN file
			if (status != 0) throw status;

			mkl_free(A22i);

			/* Optional part related to full H_inv*/

			if (isH) {

				runtimelog (true, "Making full H(-1). Restore full A(-1) ...", 0, false, false);

				try {
					strAinv *a; a = (strAinv *)malloc( sizeof( strAinv ));
					if (a == NULL) throw 40;
					dKey k;
					fAfi.seekg(0, std::ios::beg);
					for (size_t i = 0; i < ainv_sz; i++) {
						fAfi.read(reinterpret_cast<char*> (a), sizeof( strAinv ));
						k.day =  a->id_1;
						k.id = a->id_2;
						ainv[k] = a->val;
					}
					free(a);
				}
				catch (int ex) {
					write_log("Memory for a", ex);
					throw ex;
				}
				catch (std::exception const& e) {
					write_log (e.what(), 3);
					write_log ("Read ainv (s.step part)", 3);
					throw 3;
				}

				runtimelog (true, "Build H(-1) = A(-1) + [G-A22] ...", 0, false, false);

				// this is for half-store format
				for (size_t i = 1; i <= genotypedID.size(); i++) {
					dKey k;
					k.day = genotypedID[i-1];
					for (size_t j = 1; j <=i; j++) {
						k.id = genotypedID[j-1];
						ainv[k] = ainv[k] + G[i*(i-1)/2 + j-1];
					}
				}

				mkl_free(G);

				runtimelog (true, "Writing to a file full H(-1) ...", 0, false, false);

				status = write_toFile2 ("H_mat", ainv);
				if (status != 0) throw status;


				/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Gi-A22i: \n");
			for (auto i=1; i<=min(genotypedID.size(),26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", G[i*(i-1)/2 + j-1]);
				}
				printf ("\n");
			}
			std::cout<<std::endl;
			std::cout<<"full A inv:"<<std::endl;
			auto debug_count = 0;
			for(auto const& elem : ainv) {
				if (debug_count < 200) { std::cout << "key: "<<elem.first.day <<" "<<elem.first.id<< ",  values: "<<elem.second<<"\n"; debug_count++; }
			}
			std::cout<<std::endl;
			//----- End DEBUG PRINT ---------------------------------------------------
				 */

				ainv.clear();

			}
			else
				mkl_free(G);

			//fA22.close(); remove(filename1.c_str());
			fA22i.close(); remove(filename2.c_str());
			fAfi.close(); remove(filename3.c_str());
			//fAri.close(); remove(filename4.c_str());
			fGcci.close(); remove(filename5.c_str());
			fGnc.close(); remove(filename6.c_str());
			fGnni.close(); remove(filename7.c_str());
			fGcn.close(); remove(filename8.c_str());
			fG12.close(); remove(filename9.c_str());
			fGi.close(); remove(filename10.c_str());
			//fa11.close(); remove(filename11.c_str());
			//fa21.close(); remove(filename12.c_str());
			//fa22.close(); remove(filename13.c_str());
			fa22i.close(); remove(filename14.c_str());
			//fa12.close(); remove(filename15.c_str());

			runtimelog (false, "Execution status", 0, false, true);

			return 0;

			/*=============== END OF SINGLE STEP BLOCK ==============================================================*/

		}
		else {

			/*=============== APY BLOCK ========================================================================*/

			runtimelog (true, "Make temporal L-part of G matrix of single precision, will be used to build the G sub-matrices", 0, false, false);

			/* Make temporal G matrix of single precision, will be used to build the G sub-matrices */

			float *Gf;
			status = sstepmat.make_res_array_h(Gf, genotypedID.size());
			if (status != 0) {
				write_log("Memory for G", status);
				throw status;
			}

			std::map<int, int> rid_map;
			status = getRecodedIdMap(rid_map, genotypedID);
			if (status != 0) throw status;

			//unsigned int n_threads = std::thread::hardware_concurrency();
			unsigned int block_size = static_cast<unsigned int> (g_r.size()/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for(size_t i = 0; i < g_r.size(); i++) {
				size_t r = rid_map[g_r[i]];
				size_t c = rid_map[g_c[i]];
				size_t ind = r*(r - 1)/2 + c -1;
				if (c > r) ind = c*(c - 1)/2 + r -1;
				Gf[ind] = g_v[i];
			}

			g_r.clear();
			g_r.shrink_to_fit();
			g_c.clear();
			g_c.shrink_to_fit();
			g_v.clear();
			g_v.shrink_to_fit();
			rid_map.clear();

			// Here APY starts;
			// we start with G part since at this moment there is g_container in the memory.

			runtimelog (true, "Making complete (scaled) G ...", 0, false, false);

			// Complete scaling of G matrix
			// this is the fastest way, and it is commented out since it uses A22 all in RAM

			//----- Start DEBUG PRINT -------------------------------------------------

			/*
			// Restore A22 from the file
			status = sstepmat.make_res_array_h (Agen, genotypedID.size());
			if (status != 0) {
				write_log("Memory for A22", status);
				throw status;
			}

			try {
				fA22.seekg(0, std::ios::beg);
				size_t sz = (genotypedID.size()*genotypedID.size() + genotypedID.size())/2;
				fA22.read(reinterpret_cast<char*>(Agen), sz*sizeof( float ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Read A22 from the file", 3);
				throw 3;
			}

			//----- Start DEBUG PRINT -------------------------------------------------

			printf ("\n Agen: \n");
			for (auto i=1; i<=min(genotypedID.size(),26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", Agen[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}



			printf ("\n Gf: \n");
			for (auto i=1; i<=min(genotypedID.size(),26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", Gf[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			 */


			//----- End DEBUG PRINT ---------------------------------------------------

			/*

			// complete scaling
#pragma omp parallel for
			for(auto i = 1; i <= genotypedID.size(); i++) {
				for(auto j = 1; j <= i; j++) {
					size_t ind = i*(i - 1)/2 + j - 1;
					double g = Gf[ind];
					if (g == 0.0) g = alpha;
					Gf[ind] = static_cast<float> ((1-wa) * g  + wa * Agen[ind]);
				}
			}


			mkl_free(Agen);
			 */

			//----- End DEBUG PRINT ---------------------------------------------------


			/* Complete scaling of G matrix
			This is much slower way since we read A22 element-by-element directly from disc, which is not in parallel
			This is just in order to save RAM */

			try {
				fA22.seekg(0, std::ios::beg);
				float *a; a = (float *)malloc( 1*sizeof( float ));
				if (a == NULL) throw 40;
				size_t lda = genotypedID.size();
				for (size_t i = 0; i < (lda*lda + lda)/2; i++) {
					fA22.read(reinterpret_cast<char*> (a), 1*sizeof( float ));
					if (Gf[i] == 0.0) Gf[i] = alpha;
					Gf[i] = static_cast<float> ((1-wa) * Gf[i]  + wa * a[0]);
				}
				free(a);
			}
			catch (int ex) {
				write_log("Memory for a", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Read A22 from the file", 3);
				throw 3;
			}

			fA22.close(); remove(filename1.c_str());

			/*
			//----- Start DEBUG PRINT -----------------------------------------------------------
			printf ("\n complete G: \n");
			for (auto i=1; i<=min(genotypedID.size(),26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", Gf[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT -------------------------------------------------------------
			 */

			runtimelog (true, "Making full store Gcc ...", 0, false, false);

			// make Gcc
			double *Gcc;

			size_t r_gcc, c_gcc;
			r_gcc = c_gcc = coreID.size();

			status = sstepmat.make_res_array (Gcc, r_gcc, c_gcc);
			if (status != 0) {
				write_log("Memory for Gcc", status);
				throw status;
			}

			// make the list of positions of coreIDs in genotypedIDs
			std::map<int, int> corePositions;
			status = findRecodedIdMap(corePositions, genotypedID, coreID);
			if (status != 0) throw status;

			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (coreID.size()/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for ( size_t i = 0; i < coreID.size(); i++ ) {
				//auto r = find_invect2(genotypedID, coreID[i])+1;
				size_t r = corePositions[ coreID[i] ];
				for (size_t j = 0; j <= i; j++) {
					//auto c = find_invect2(genotypedID, coreID[j])+1;
					size_t c = corePositions[ coreID[j] ];
					//Gcc[i * coreID.size() + j] = Gf[(r - 1) * genotypedID.size() + c - 1]; // for ful-store Gf
					Gcc[i * coreID.size() + j] = Gcc[j * coreID.size() + i] = static_cast<double> (Gf[r*(r - 1)/2 + c - 1]); // for half-store Gf
				}
			}

			corePositions.clear();

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of Gcc: \n");
			for (auto i=0; i<min(r_gcc,26); i++) {
				for (auto j=0; j<min(c_gcc,26); j++) {
					printf ("%12.5G", Gcc[j+i*c_gcc]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			runtimelog (true, "Inverting Gcc ...", 0, false, false);

			lapack_int *gpvt;

			status = sstepmat.make_pivot_array(gpvt, r_gcc);
			if (status != 0) {
				write_log("Memory for gpvt", status);
				throw status;
			}

			status = sstepmat.matr_inv(Gcc, gpvt, r_gcc, c_gcc);
			if (status != 0) {
				write_log("Inversion of Gcc", status);
				throw status;
			}

			mkl_free(gpvt);

			/*

			// For 'L' part of Gcc:

			double *Gcc;

			int r_gcc, c_gcc;
			r_gcc = c_gcc = coreID.size();

			sstepmat.make_res_array_h (Gcc, r_gcc);

#pragma omp parallel for
			for ( auto i = 1; i <= coreID.size(); i++ ) {
				auto r = find_invect2(genotypedID, coreID[i-1])+1;
				for (auto j = 1; j <= i; j++) {
					auto c = find_invect2(genotypedID, coreID[j-1])+1;
					Gcc[i*(i - 1)/2 + j - 1] = static_cast<double> (Gf[r*(r - 1)/2 + c - 1]); // for half-store Gf
				}
			}

			finish = std::chrono::high_resolution_clock::now();
			elapsed = finish - start;
			std::cout << "Elapsed time of making making full store Gcc:" << elapsed.count() << " s\n\n";
			std::cout<<std::endl<<std::endl;
			std::cout << "The size of Gcc in mmemory, GiB:"<< static_cast<double> (sizeof (double[r_gcc*c_gcc])) /(1024*1024*1024) << "\n\n";

			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of 'L' part of Gcc: \n");
			for (auto i=1; i<=min(r_gcc,26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", Gcc[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------

			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<"invert 'L' part of Gcc:"<<std::endl<<std::endl;
			//----- End DEBUG PRINT ---------------------------------------------------

			start = std::chrono::high_resolution_clock::now();

			// invert Gcc
			sstepmat.matr_inv_h(Gcc, r_gcc);

			finish = std::chrono::high_resolution_clock::now();
			elapsed = finish - start;
			std::cout << "Elapsed time of 'L' part Gcc inverse:" << elapsed.count() << " s\n\n";

			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of 'L' part of inverted Gcc: \n");
			for (auto i=1; i<=min(r_gcc,26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", Gcc[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of 'L' part of inverted Gcc: \n");
			for (auto i=1; i<=min(r_gcc,26); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", Gcc[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			runtimelog (true, "Write out Gcc to a file ...", 0, false, false);

			try {
				fGcci.write(reinterpret_cast<char*>(Gcc), r_gcc*c_gcc*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Writing Gcc inverse to the file", 3);
				throw 3;
			}

			mkl_free(Gcc);

			runtimelog (true, "Make Gnc ...", 0, false, false);

			double *Gnn, *Gnc;

			size_t r_gnn, c_gnn;
			size_t r_gnc, c_gnc;
			r_gnn = c_gnn = genotypedID.size() - coreID.size();
			r_gnc = r_gnn;
			c_gnc = c_gcc;

			runtimelog (true, "Make a vector of non-core IDs ...", 0, false, false);

			std::vector<int> noncoreID;

			for (auto i = 0; i < genotypedID.size(); i++) {
				int id = genotypedID[i];
				if (!find_invect(coreID, id)) noncoreID.push_back(id);
			}

			//----- Start DEBUG PRINT -------------------------------------------------
			/*
			for (auto const& i : noncoreID)
				std::cout<<"noncore ID : "<<i<<std::endl;
			 */
			//----- End DEBUG PRINT ---------------------------------------------------

			std::sort(noncoreID.begin(), noncoreID.end());

			status = sstepmat.make_res_array (Gnc, r_gnc, c_gnc);
			if (status != 0) {
				write_log("Memory for Gnc", status);
				throw status;
			}

			runtimelog (true, "Build Gnc from G matrix ...", 0, false, false);

			// make combined vector coreNoncoreIDs
			std::vector<int> coreNonCoreIDs;

			coreNonCoreIDs.insert(coreNonCoreIDs.end(), coreID.begin(), coreID.end());

			coreNonCoreIDs.insert(coreNonCoreIDs.end(), noncoreID.begin(), noncoreID.end());

			status = findRecodedIdMap(corePositions, genotypedID, coreNonCoreIDs);
			if (status != 0) throw status;

			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (noncoreID.size()/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for ( size_t i = 0; i < noncoreID.size(); i++ ) {
				//auto r = find_invect2(genotypedID, noncoreID[i])+1;
				size_t r = corePositions[ noncoreID[i] ];
				for (size_t j = 0; j < coreID.size(); j++) {
					//auto c = find_invect2(genotypedID, coreID[j])+1;
					size_t c = corePositions[ coreID[j] ];
					size_t ind;
					if (noncoreID[i] >= coreID[j]) ind = r*(r - 1)/2 + c -1;
					else ind = c*(c - 1)/2 + r -1;
					Gnc[i * coreID.size() + j] = static_cast<double> (Gf[ind]); // for half-stored Gf
				}
			}

			coreNonCoreIDs.clear();

			runtimelog (true, "Write out Gnc to a file ...", 0, false, false);

			try {
				fGnc.write(reinterpret_cast<char*>(Gnc), r_gnc*c_gnc*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Writing Gnc to the file", 3);
				throw 3;
			}

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of Gnc: \n");
			for (auto i=0; i<min(r_gnc,26); i++) {
				for (auto j=0; j<min(c_gnc,26); j++) {
					printf ("%12.5G", Gnc[j+i*c_gnc]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */


			/*========================== Fast but memory inefficient way to calculate Gnn =============================================*/

			runtimelog (true, "Make Gcn (transpose Gnc to get Gcn) and write it to a file ...", 0, false, false);

			double *Gcn;

			status = sstepmat.make_res_array (Gcn, c_gnc, r_gnc);
			if (status != 0) {
				write_log("Memory for Gcn", status);
				throw status;
			}

			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (r_gnc/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i = 0; i < r_gnc; i++) {
				for (size_t j = 0; j < c_gnc; j++) {
					Gcn[j*r_gnc+i] = Gnc[i*c_gnc+j];
				}
			}

			try {
				fGcn.write(reinterpret_cast<char*>(Gcn), r_gnc*c_gnc*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Writing Gcn to the file", 3);
				return 3;
			}

			mkl_free(Gcn);

			runtimelog (false, "Calculate Gnn and then - invert it ...", 0, false, false);

			double *GncGcc;

			status = sstepmat.make_res_array (GncGcc, r_gnc, c_gnc);
			if (status != 0) {
				write_log("Memory for GncGcc", status);
				throw status;
			}

			// restore Gcc from file:
			try {
				status = sstepmat.make_res_array (Gcc, r_gcc, c_gcc);
				if (status != 0) throw status;
				fGcci.seekg(0, std::ios::beg);
				fGcci.read(reinterpret_cast<char*>(Gcc), r_gcc*c_gcc*sizeof( double ));
			}
			catch (int ex) {
				write_log ("Memory for Gcc", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading Gcci file", 3);
				throw 3;
			}

			sstepmat.matr_prod(Gnc, Gcc, GncGcc, r_gnc, r_gcc, r_gcc, c_gcc);

			mkl_free(Gcc);
			mkl_free(Gnc);

			// restore Gcn from file:
			try {
				status = sstepmat.make_res_array (Gcn, c_gnc, r_gnc);
				if (status != 0) throw status;
				fGcn.seekg(0, std::ios::beg);
				fGcn.read(reinterpret_cast<char*>(Gcn), r_gnc*c_gnc*sizeof( double ));
			}
			catch (int ex) {
				write_log ("Memory for Gcn", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading Gcn file", 3);
				throw 3;
			}

			double *Gnnfull;

			status = sstepmat.make_res_array (Gnnfull, r_gnc, r_gnc);
			if (status != 0) {
				write_log("Memory for Gnnfull", status);
				throw status;
			}

			sstepmat.matr_prod(GncGcc, Gcn, Gnnfull, r_gnc, c_gnc, r_gcc, r_gnc);

			mkl_free(Gcn);
			mkl_free(GncGcc);

			status = sstepmat.make_res_array (Gnn, r_gnn, 1); // it has the number only on 'D', so make it one-dimentional representation
			if (status != 0) {
				write_log("Memory for Gnn", status);
				throw status;
			}

			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (noncoreID.size()/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i = 1; i <= noncoreID.size(); i++) {
				size_t r = corePositions[ noncoreID[i-1] ];
				Gnn[i - 1] = 1.0 / (static_cast<double> (Gf[r*(r - 1)/2 + r - 1]) - Gnnfull[(i - 1) * r_gnc + (i - 1)]); // half-store case
			}

			mkl_free(Gnnfull);
			mkl_free(Gf);

			try {
				fGnni.write(reinterpret_cast<char*>(Gnn), r_gnn*1*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception writing Gnni to the file", 3);
				throw 3;
			}

			/*========================== End of fast but memory inefficient way to calculate Gnn =============================================*/

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n complete Gnn_i (new): \n");
			for (auto i=0; i<min(r_gnn,26); i++) {
				for (auto j=0; j<=i; j++) {
					if (i == j) printf ("%12.5G", Gnn[i]);
					else printf ("%12.5G", 0.0);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			/*
			//====================== slow but memory efficient way to produce Gnn_i ========================================================

			// manage threads:
		    n_threads = std::thread::hardware_concurrency();
		    unsigned int loop_threads = 1;
		    unsigned int matr_threads = 1;
		    if (n_threads > 2) {
		    	loop_threads = 2;
		    	matr_threads = n_threads / loop_threads;
		    	while (matr_threads > loop_threads) {
		    		loop_threads += 2;
		    		matr_threads = n_threads / loop_threads;
		    	}
		    }

		    block_size = static_cast<unsigned int> (noncoreID.size()/(64*loop_threads));

		    std::cout << n_threads << " concurrent threads are supported.\n";
		    std::cout << loop_threads << " loop threads are supported.\n";
		    std::cout << matr_threads << " matr threads are supported.\n";

			// get Gnn and invert it
			omp_set_nested(1);
#pragma omp parallel for schedule(dynamic, block_size) num_threads(loop_threads)
			for (auto i = 1; i <= noncoreID.size(); i++) {

				auto r = corePositions[ noncoreID[i-1] ];
				double *c1; // c1 = Gcc_i * Gic_t
				double *gic; //gic = Gnc(i, 1:c_gnc);
				double *c2; // c2 = Gic * c1
				sstepmat.make_res_array (c1, r_gcc, 1);
				sstepmat.make_res_array (c2, 1, 1);
				sstepmat.make_res_array (gic, c_gnc, 1);

				for (auto k = 1; k <= c_gnc; k++) {
					gic[k-1] = Gnc[(i - 1) * c_gnc + (k - 1)];
				}

#pragma omp parallel for schedule(dynamic) num_threads(matr_threads)
				for (auto l = 0; l < r_gcc; l++) {
					c1[l] = 0.0;
					for (auto j = 0; j < c_gnc; j++) {
						c1[l] += Gcc[j+l*c_gcc]*gic[j];
					}
				}

				c2[0] = 0.0;
				for (auto l = 0; l < c_gnc; l++) {
					c2[0] += gic[l] * c1[l];
				}

				//Gnn[(i - 1)*(r_gnn + 1)] = 1.0 / (Gf[(r - 1) * genotypedID.size() + r - 1] - c2[0]); // full-store case
				Gnn[i - 1] = 1.0 / (static_cast<double> (Gf[r*(r - 1)/2 + r - 1]) - c2[0]); // half-store case
			}

			omp_set_nested(0);
			corePositions.clear();

			finish = std::chrono::high_resolution_clock::now();
			elapsed = finish - start;
			std::cout << "Elapsed time of making Gnn (inverted):" << elapsed.count() << " s\n\n";

			try {
				//fGnni.write(reinterpret_cast<char*>(Gnn), r_gnn*c_gnn*sizeof( double ));
				fGnni.write(reinterpret_cast<char*>(Gnn), r_gnn*1*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (where, 10);
				std::cerr << "Exception writing file: "<< e.what() << std::endl;
				return 10;
			}


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n complete Gnn_i: \n");
			for (auto i=0; i<min(r_gnn,26); i++) {
				for (auto j=0; j<=i; j++) {
					if (i == j) printf ("%12.5G", Gnn[i]);
					else printf ("%12.5G", 0.0);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------


			mkl_free(Gf); // we will not use G any more
			//=========================================================================================================
			 */

			runtimelog (true, "Restore Gcn and Gcc from the files ...", 0, false, false);

			try {
				status = sstepmat.make_res_array (Gcn, c_gnc, r_gnc);
				if (status != 0) throw status;
				fGcn.seekg(0, std::ios::beg);
				fGcn.read(reinterpret_cast<char*>(Gcn), r_gnc*c_gnc*sizeof( double ));
			}
			catch (int ex) {
				write_log("Memory for Gcn", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading Gcn from the file", 3);
				throw 3;
			}

			fGcn.close(); remove(filename8.c_str());

			try {
				status = sstepmat.make_res_array (Gcc, r_gcc, c_gcc);
				if (status != 0) throw status;
				fGcci.seekg(0, std::ios::beg);
				fGcci.read(reinterpret_cast<char*>(Gcc), r_gcc*c_gcc*sizeof( double ));
			}
			catch (int ex) {
				write_log("Memory for Gcci", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading Gcci from the file", 3);
				throw 3;
			}

			// make G12: G12 = - Gcc_i * Gcn * Gnn_i.
			// 1) G12 = Gcc_i * Gcn;
			// 2) G12 = G12 * Gnn_i;
			// 3) G12 = -1 * G12;

			runtimelog (true, "Matrix product, produce Gcci*Gcn ...", 0, false, false);

			double *G12;
			status = sstepmat.make_res_array (G12, r_gcc, r_gnc);// !!!
			if (status != 0) {
				write_log("Memory for G12", status);
				throw status;
			}

			sstepmat.matr_prod(Gcc, Gcn, G12, r_gcc, c_gnc, c_gcc, r_gnc);

			mkl_free(Gcc);
			mkl_free(Gcn);

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of Gcc_cn: \n");
			for (auto i=0; i<min(r_gcc,26); i++) {
				for (auto j=0; j<min(r_gnc,26); j++) {
					printf ("%12.5G", G12[j+i*r_gnc]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			//double *G12;
			size_t r_g12 = r_gcc;
			size_t c_g12 = c_gnn;
			//sstepmat.make_res_array (G12, r_g12, c_g12);

			runtimelog (true, "Matrix product, produce Gcc_cn*Gnn (G12 matrix) ...", 0, false, false);

			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (r_gcc/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i=0; i<r_gcc; i++) {
				for (size_t j=0; j<r_gnc; j++) {
					G12[j+i*r_gnc] = G12[j+i*r_gnc]*Gnn[j];
				}
			}

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of Gcc_cn * Gnn: \n");
			for (auto i=0; i<min(r_gcc,26); i++) {
				for (auto j=0; j<min(r_gnc,26); j++) {
					printf ("%12.5G", G12[j+i*r_gnc]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			mkl_free(Gnn);

			/*
			make G11: G11 = (I + G12 * Gnc) * Gcc_i.
			1) G12_nc = G12 * Gnc.
			2) G12_nc = G12_nc + I.
										not exists anymore 2) G12_nc_cc = G12_nc * Gcc_i.
			3) G11 = G12_nc * Gcc_i.
			 							not exists anymore 3) G11 = Gcc_i + G12_nc_cc.
			 */

			// restore Gnc
			try {
				status = sstepmat.make_res_array (Gnc, r_gnc, c_gnc);
				if (status != 0) throw status;
				fGnc.seekg(0, std::ios::beg);
				fGnc.read(reinterpret_cast<char*>(Gnc), r_gnc*c_gnc*sizeof( double ));
			}
			catch (int ex) {
				write_log ("Memory for Gnc", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading Gnc from the file", 3);
				throw 3;
			}

			fGnc.close(); remove(filename6.c_str());

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of Gnc: \n");
			for (auto i=0; i<min(r_gnc,26); i++) {
				for (auto j=0; j<min(c_gnc,26); j++) {
					printf ("%12.5G", Gnc[j+i*c_gnc]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			runtimelog (true, "Produce G11 ...", 0, false, false);

			double *G12_nc;

			status = sstepmat.make_res_array (G12_nc, r_g12, c_gnc);
			if (status != 0) {
				write_log("Memory for G12_nc", status);
				throw status;
			}

			sstepmat.matr_prod(G12, Gnc, G12_nc, r_g12, r_gnc, c_g12, c_gnc);

			mkl_free(Gnc);

			// complete making G12 by mult. by -1.
			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> ((r_g12 * c_g12)/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i = 0; i < (r_g12 * c_g12); i++)
				G12[i] = G12[i] * (-1.0);


			try {
				fG12.write(reinterpret_cast<char*>(G12), r_g12*c_g12*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception writing G12 to the file", 3);
				throw 3;
			}

			mkl_free(G12);

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of G12_nc: \n");
			for (auto i=0; i<min(r_g12,26); i++) {
				for (auto j=0; j<min(c_gnc,26); j++) {
					printf ("%12.5G", G12_nc[j+i*c_gnc]);
				}
				printf ("\n");
			}
			//----- Start DEBUG PRINT -------------------------------------------------
			 */


			// adding identity matrix to G12_nc
			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (r_g12/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i=0; i<r_g12; i++) {
				G12_nc[i+i*c_gnc] = G12_nc[i+i*c_gnc] + 1.0;
			}

			// restore Gcc
			try {
				status = sstepmat.make_res_array (Gcc, r_gcc, c_gcc);
				if (status != 0) throw status;
				fGcci.seekg(0, std::ios::beg);
				fGcci.read(reinterpret_cast<char*>(Gcc), r_gcc*c_gcc*sizeof( double ));
			}
			catch (int ex) {
				write_log("Memory for Gcci", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading Gcci from the file", 3);
				throw 3;
			}

			fGcci.close(); remove(filename5.c_str());

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of Gcc: \n");
			for (auto i=0; i<min(r_gcc,26); i++) {
				for (auto j=0; j<min(c_gcc,26); j++) {
					printf ("%12.5G", Gcc[j+i*c_gcc]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */


			/*
			double *G12_nc_cc;
			sstepmat.make_res_array (G12_nc_cc, r_gcc, c_gnc);
			sstepmat.matr_prod(G12_nc, Gcc, G12_nc_cc, r_gcc, r_gcc, c_gnc, c_gnc);

			mkl_free(G12_nc);
			 */

			double *G11;
			status = sstepmat.make_res_array (G11, r_gcc, c_gcc);
			if (status != 0) {
				write_log("Memory for G11", status);
				throw status;
			}

			sstepmat.matr_prod(G12_nc, Gcc, G11, r_gcc, r_gcc, c_gnc, c_gnc);
			//sstepmat.matr_add (Gcc, G12_nc_cc, G11, false, false, r_gcc, c_gcc);

			mkl_free(Gcc);
			mkl_free(G12_nc);

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of G11: \n");
			for (auto i=0; i<min(r_gcc,26); i++) {
				for (auto j=0; j<min(c_gcc,26); j++) {
					printf ("%12.5G", G11[j+i*c_gcc]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			runtimelog (true, "Produce G(-1) matrix directly to a binary file ...", 0, false, false);

			// restore G12
			try {
				status = sstepmat.make_res_array (G12, r_g12, c_g12);
				if (status != 0) throw status;
				fG12.seekg(0, std::ios::beg);
				fG12.read(reinterpret_cast<char*>(G12), r_g12*c_g12*sizeof( double ));
			}
			catch (int ex) {
				write_log("Memory for G12", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading G12 from the file", 3);
				throw 3;
			}

			fG12.close(); remove(filename9.c_str());

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of G12: \n");
			for (auto i=0; i<min(r_g12,26); i++) {
				for (auto j=0; j<min(c_g12,26); j++) {
					printf ("%12.5G", G12[j+i*c_g12]);
				}
				printf ("\n");
			}
			std::cout<<std::endl;
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			// Produce Ginv matrix directly to a binary file:

			double *gVal;
			gVal = (double *)malloc( (coreID.size()+noncoreID.size())*sizeof( double ));
			if (gVal == NULL) throw 40;
			//std::cout<<"size of gVal = "<<gVal.size()<<std::endl;


			// write out to a file core x noncore half of G_inv matrix
			// writing into the file row by row
			for (size_t i = 0; i < coreID.size(); i++) {
#pragma omp parallel sections
				{
#pragma omp section
					{
						for (size_t j = 0; j < coreID.size(); j++) {
							gVal[j] = G11[i*coreID.size()+j];
						}
					}
#pragma omp section
					{
						for (size_t k = 0; k < noncoreID.size(); k++) {
							gVal[coreID.size()+k] = G12[i*noncoreID.size()+k];
						}
					}
				}

				try {
					fGi.write(reinterpret_cast<char*> (gVal), (coreID.size()+noncoreID.size())*sizeof( double ));
				}
				catch (std::exception const& e) {
					write_log (e.what(), 4);
					write_log ("Exception writing Gci to the file", 4);
					throw 4;
				}

			}

			mkl_free(G11);


			// making G21: G21 = transpose(G12).

			double *G21;

			status = sstepmat.make_res_array (G21, c_g12, r_g12);
			if (status != 0) {
				write_log("Memory for G21", status);
				throw status;
			}

			//n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (r_g12/(1*n_threads));

#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i = 0; i < r_g12; i++) {
				for (size_t j = 0; j < c_g12; j++) {
					G21[j*r_g12+i] = G12[i*c_g12+j];
				}
			}

			mkl_free(G12);

			/*
			// making G21: G21 = transpose(G12).
			double *G21;
			sstepmat.make_res_array (G21, r_g12, c_g12);
			for (auto i = 0; i < (r_g12 * c_g12); i++)
				G21[i] = G12[i];

			mkl_free(G12);

			start = std::chrono::high_resolution_clock::now();
			// transpose G12 to get G21
			mkl_dimatcopy ('r', 't', r_g12, c_g12, 1.0, G21, c_g12, r_g12);
			finish = std::chrono::high_resolution_clock::now();
			elapsed = finish - start;
			std::cout << "Elapsed time of transpose G12 to get G21:" << elapsed.count() << " s\n\n";
			 */


			// restore Gnn
			try {
				status = sstepmat.make_res_array (Gnn, r_gnn, 1);
				if (status != 0) throw status;
				fGnni.seekg(0, std::ios::beg);
				fGnni.read(reinterpret_cast<char*>(Gnn), r_gnn*1*sizeof( double ));
			}
			catch (int ex) {
				write_log("Memory for Gnni", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading Gnni from the file", 3);
				throw 3;
			}

			fGnni.close(); remove(filename7.c_str());

			// write out to a file noncore(core) x noncore half of G_inv matrix
			for (size_t i = 0; i < noncoreID.size(); i++) {
#pragma omp parallel sections
				{
#pragma omp section
					{
						for (size_t j = 0; j < coreID.size(); j++) {
							gVal[j] = G21[i*coreID.size()+j];
						}
					}
#pragma omp section
					{
						for (size_t k = 0; k < noncoreID.size(); k++) {
							//gVal[coreID.size()+k] = Gnn[i*noncoreID.size()+k];
							if (i == k) gVal[coreID.size()+k] = Gnn[i];
							else gVal[coreID.size()+k] = 0.0;
						}
					}
				}
				//for (auto i = 0; i < (coreID.size()+noncoreID.size()); i++) printf ("%12.5G", gVal[i]);//std::cout<<gVal[i]<<" ";
				//std::cout<<std::endl;

				try {
					fGi.write(reinterpret_cast<char*> (gVal), (coreID.size()+noncoreID.size())*sizeof( double ));
				}
				catch (std::exception const& e) {
					write_log (e.what(), 3);
					write_log ("Exception writing Gi to the file", 3);
					throw 3;
				}

			}


			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n complete Gnn_i: \n");
			for (auto i=0; i<min(r_gnn,26); i++) {
				for (auto j=0; j<=i; j++) {
					if (i == j) printf ("%12.5G", Gnn[i]);
					else printf ("%12.5G", 0.0);
				}
				printf ("\n");
			}

			printf ("\n Top left corner of G21: \n");
			for (auto i=0; i<min(c_g12,26); i++) {
				for (auto j=0; j<min(r_g12,26); j++) {
					printf ("%12.5G", G21[j+i*r_g12]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			free(gVal);
			mkl_free(G21);
			mkl_free(Gnn);
			mkl_free_buffers();

			// make Ginv just for debugging purpose since we already have Ginv in the file
			/*
			{
				double *Ginv;
				int r_ginv = genotypedID.size();
				int c_ginv = genotypedID.size();
				sstepmat.make_res_array (Ginv, r_ginv, c_ginv);

				gVal = (double *)malloc( (coreID.size()+noncoreID.size())*sizeof( double ));

				fGi.seekg(0, std::ios::beg);

				// read from a file core x noncore half of G_inv matrix
				for (auto i = 0; i < coreID.size(); i++) {
					try {
						fGi.read(reinterpret_cast<char*> (gVal), (coreID.size()+noncoreID.size())*sizeof( double ));
					}
					catch (std::exception const& e) {
						write_log (where, 10);
						std::cerr << "Exception writing to file: "<< e.what() << std::endl;
						return 10;
					}
					for (auto m = 0; m < (coreID.size()+noncoreID.size()); m++)
						Ginv[i*genotypedID.size()+m] = gVal[m];
				}
				// read from a file noncore(core) x noncore half of G_inv matrix
				for (auto i = 0; i < noncoreID.size(); i++) {
					try {
						fGi.read(reinterpret_cast<char*> (gVal), (coreID.size()+noncoreID.size())*sizeof( double ));
					}
					catch (std::exception const& e) {
						write_log (where, 10);
						std::cerr << "Exception writing to file: "<< e.what() << std::endl;
						return 10;
					}
					for (auto m = 0; m < (coreID.size()+noncoreID.size()); m++)
						Ginv[(coreID.size()+i)*genotypedID.size()+m] = gVal[m];

				}

				//----- Start DEBUG PRINT -------------------------------------------------
				printf ("\n Top left corner of restored Ginv: \n");
				for (auto i=0; i<min(r_ginv,20); i++) {
					for (auto j=0; j<min(c_ginv,20); j++) {
						printf ("%12.5G", Ginv[j+i*c_ginv]);
					}
					printf ("\n");
				}
				//----- End DEBUG PRINT ---------------------------------------------------

				mkl_free(Ginv);
			}
			 */

			/*
			// ==========  working with reduced A_inv :

			runtimelog (true, "Making A22(-1). Prepare vectors and restore r_ainv from the file ...", 0, false, false);

			size_t r_ginv = genotypedID.size();
			size_t c_ginv = genotypedID.size();

			// Restore r_ainv
			try {
				dKey k;
				fAri.seekg(0, std::ios::beg);
				strAinv *a; a = (strAinv *)malloc( sizeof( strAinv ));
				if (a == NULL) throw 40;
				for (auto i = 0; i < r_ainv_size; i++) {
					fAri.read(reinterpret_cast<char*> (a), sizeof( strAinv ));
					k.day = a->id_1;
					k.id = a->id_2;
					r_ainv[k] = a->val;
				}
				free(a);
			}
			catch (int ex) {
				write_log("Memory for a", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading r_ainv from the file", 3);
				throw 3;
			}


			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<"r_ainv:"<<std::endl;
			for (auto const& elem: r_ainv) {
				std::cout<<elem.first.day<<" "<<elem.first.id<<" "<<elem.second<<std::endl;
			}
			//----- End DEBUG PRINT ---------------------------------------------------


			runtimelog (true, "Build A11 and write it to a file ...", 0, false, false);

			double *A11; // non-genotyped IDs matrix

			size_t r_a11 = nongenotID.size();
			size_t c_a11 = nongenotID.size();

			// give the range of IDs to define the submatrix (USE the ORIGINAL IDs from the pedigree)

			//sstepmat.make_res_array (A11, r_a11, c_a11); // full-store matrix
			status = sstepmat.make_res_array_h (A11, r_a11); //half-store matrix
			if (status != 0) {
				write_log("Memory for A11", status);
				throw status;
			}

			n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (nongenotID.size()/(1*n_threads));

			//size_t a11_el = 0;
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i = 1; i <= nongenotID.size(); i++ ) {
				dKey k;
				k.day = nongenotID[i-1];
				for (size_t j = 1; j <= i; j++ ) {
					k.id = nongenotID[j-1];
					if (r_ainv.count(k)) {
						A11[i*(i - 1)/2 + j - 1] = r_ainv[k];
						//a11_el++;
					}
				}
			}
			//std::cout << "a11_el = " <<a11_el << "\n";

			try {
				fa11.write(reinterpret_cast<char*>(A11), ((r_a11*r_a11 + r_a11)/2)*sizeof( double )); //half-store matrix
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception writing A11 to the file", 3);
				throw 3;
			}

			runtimelog (false, "The size of A11 in mmemory, GiB", static_cast<double> ( sizeof (double[(r_a11*r_a11 + r_a11)/2]) ) /(1024*1024*1024), false, false);


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of matrix A11: \n");
			for (auto i=1; i<=min(nongenotID.size(),20); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", A11[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------


			mkl_free(A11);

			runtimelog (true, "Build A22 and write it to a file ...", 0, false, false);

			double *A22; // genotyped

			size_t r_a22 = genotypedID.size();
			size_t c_a22 = genotypedID.size();

			status = sstepmat.make_res_array_h (A22, r_a22); // half-store matrix
			if (status != 0) {
				write_log("Memory for A22", status);
				throw status;
			}

			n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (genotypedID.size()/(1*n_threads));

			//size_t a22_el = 0;
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i = 1; i <= genotypedID.size(); i++ ) {
				dKey k;
				k.day = genotypedID[i-1];
				for (size_t j = 1; j <= i; j++ ) {
					k.id = genotypedID[j-1];
					//A22[(i - 1)*genotypedID.size() + j - 1] = A22[(j - 1)*genotypedID.size() + i - 1] = r_ainv[k]; // full-store matrix
					if (r_ainv.count(k)) {
						A22[i*(i - 1)/2 + j - 1] = r_ainv[k];
						//a22_el++;
					}
				}
			}
			//std::cout << "a22_el = " <<a22_el << "\n";

			try {
				fa22.write(reinterpret_cast<char*>(A22), ((r_a22*r_a22 + r_a22)/2)*sizeof( double ));
				if (fa22.fail()) throw 4;
			}
			catch (int ex) {
				write_log ("Writing A22 to binary file", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception writing A22 to the file", 3);
				throw 3;
			}

			runtimelog (false, "The size of A22 in mmemory, GiB", static_cast<double> ( sizeof (double[(r_a22*r_a22 + r_a22)/2]) ) /(1024*1024*1024), false, false);


			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<std::endl;
			printf ("\n Top left corner of matrix A22: \n");
			for (auto i=1; i<=min(genotypedID.size(),29); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", A22[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------


			mkl_free(A22);

			runtimelog (true, "Build A21 and write it to a file ...", 0, false, false);

			double *A21; // genotyped - non-genotyped
			size_t c_a21 = nongenotID.size();
			size_t r_a21 = genotypedID.size();

			status = sstepmat.make_res_array (A21, r_a21, c_a21);
			if (status != 0) {
				write_log("Memory for A21", status);
				throw status;
			}

			n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (genotypedID.size()/(1*n_threads));

			//size_t elm_inA21 = 0;
#pragma omp parallel for schedule(static, block_size) num_threads(n_threads)
			for (size_t i = 1; i <= genotypedID.size(); i++ ) {
				dKey k;
				k.day = genotypedID[i-1];
				bool changed_key = false;
				for (size_t j = 1; j <= nongenotID.size(); j++ ) {
					k.id = nongenotID[j-1];
					 if k.day < k.id -> the key does not exists!
			 * Hence we will have zero element in a full matrix.
					if (k.day < k.id) {
						size_t t_id = k.id;
						k.id = k.day;
						k.day = t_id;
						changed_key = true;
					}
					if (r_ainv.count(k)) {
						A21[(i - 1)*nongenotID.size() + j - 1] = r_ainv[k];
						//elm_inA21++;
					}
					if (changed_key){
						k.day = genotypedID[i-1];
						changed_key = false;
					}
				}
			}

			//std::cout<<"elements in A21 = "<< elm_inA21<<std::endl;

			r_ainv.clear();

			try {
				fa21.write(reinterpret_cast<char*>(A21), r_a21*c_a21*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception writing A21 to the file", 3);
				throw 3;
			}

			runtimelog (false, "The size of A21 in mmemory, GiB", static_cast<double> ( sizeof (double[(r_a21*c_a21)]) ) /(1024*1024*1024), false, false);

			runtimelog (true, "Transpose A21 to get A12 ...", 0, false, false);

			// transpose A21 to get A12
			double *A12;

			status = sstepmat.make_res_array (A12, c_a21, r_a21);
			if (status != 0) {
				write_log("Memory for A12", status);
				throw status;
			}

			n_threads = std::thread::hardware_concurrency();
			block_size = static_cast<unsigned int> (r_a21/(64*n_threads));

#pragma omp parallel for schedule(dynamic, block_size) num_threads(n_threads)
			for (size_t i = 0; i < r_a21; i++) {
				for (size_t j = 0; j < c_a21; j++) {
					A12[j*r_a21+i] = A21[i*c_a21+j];
				}
			}

			try {
				fa12.write(reinterpret_cast<char*>(A12), r_a21*c_a21*sizeof( double ));
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception writing A12 to the file", 3);
				throw 3;
			}

			size_t r_a12 = c_a21;
			size_t c_a12 = r_a21;


			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<std::endl;
			printf ("\n A21: \n");
			for (auto i=0; i<min(r_a21,19); i++) {
				for (auto j=0; j<min(c_a21,19); j++) {
					printf ("%12.5G", A21[j+i*c_a21]);
				}
				printf ("\n");
			}
			std::cout<<std::endl;
			printf ("\n A12: \n");
			for (auto i=0; i<min(c_a21,26); i++) {
				for (auto j=0; j<min(r_a21,20); j++) {
					printf ("%12.5G", A12[j+i*r_a21]);
				}
				printf ("\n");
			}
			std::cout<<std::endl;

			double ma_el = 0.0;
			double mi_el = 0.0;
			size_t ma_num = 0;
			size_t mi_num = 0;
			for (size_t i = 0; i<r_a21*c_a21; i++){
				if (A12[i] != 0.0){
					if(A12[i] < mi_el){
						mi_num++;
						mi_el = A12[i];
					}
					if(A12[i] > ma_el){
						ma_num++;
						ma_el = A12[i];
					}
				}
			}
			std::cout<<"mi_num, mi_el:"<<mi_num<<", "<<mi_el<<std::endl;
			std::cout<<"ma_num, ma_el:"<<ma_num<<", "<<ma_el<<std::endl;
			std::cout<<std::endl;

			//----- End DEBUG PRINT ---------------------------------------------------


			mkl_free(A21);
			mkl_free(A12);

			runtimelog (true, "Restore A11 and inverting it ...", 0, false, false);

			try {
				status = sstepmat.make_res_array_h (A11, r_a11);
				if (status != 0) throw status;
				fa11.seekg(0, std::ios::beg);
				fa11.read(reinterpret_cast<char*>(A11), ((r_a11*r_a11 + r_a11)/2)*sizeof( double )); // half-store
			}
			catch (int ex)
			{
				write_log("Memory for Gcci", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading A11 from the file", 3);
				throw 3;
			}


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n A11: \n");
			for (auto i=1; i<=min(nongenotID.size(),29); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", A11[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------



			// for full-store matrix
			lapack_int *pvt;
			sstepmat.make_pivot_array(pvt, r_a11);
			sstepmat.matr_inv(A11, pvt, r_a11, c_a11);
			mkl_free(pvt);

			// for half-store matrix
			status = sstepmat.matr_inv_h(A11, r_a11);
			if (status != 0) {
				write_log("Inversion of A11", status);
				throw status;
			}


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n inverse of matrix A11: \n");
			for (auto i=1; i<=min(nongenotID.size(),20); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", A11[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------



			// A_11t_12 = A11_inv * A12. Rows = r_a11; columns = c_a12.

			runtimelog (true, "Make A11_i full-stored to do multiplication further on ...", 0, false, false);

			// make A11_inv full to do multiplication further on
			double *A11f;
			status = sstepmat.make_res_array (A11f, r_a11, c_a11);
			if (status != 0) {
				write_log("Memory for A11f", status);
				throw status;
			}

#pragma omp parallel for
			for (size_t i = 1; i <= nongenotID.size(); i++) {
				for (size_t j = 1; j <= i; j++) {
					A11f[(i - 1)*nongenotID.size() + j - 1] = A11f[(j - 1)*nongenotID.size() + i - 1] = A11[i*(i - 1)/2 + j - 1];
				}
			}

			mkl_free(A11); // we do not need A11 anymore

			runtimelog (false, "The size of A11_full in mmemory, GiB", static_cast<double> ( sizeof (double[(r_a11*c_a11)]) ) /(1024*1024*1024), false, false);


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of A11f: \n");
			for (auto i=0; i<min(r_a11,19); i++) {
				for (auto j=0; j<min(c_a11,19); j++) {
					printf ("%12.5G", A11f[j+i*c_a11]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------


			runtimelog (true, "A11(-1) * A12 ...", 0, false, false);

			double *A_11t_12;

			status = sstepmat.make_res_array (A_11t_12, r_a11, c_a12);
			if (status != 0) {
				write_log("Memory for A_11t_12", status);
				throw status;
			}

			// restore A12
			try {
				status = sstepmat.make_res_array (A12, c_a21, r_a21);
				if (status != 0) throw status;
				fa12.seekg(0, std::ios::beg);
				fa12.read(reinterpret_cast<char*>(A12), r_a21*c_a21*sizeof( double ));
			}
			catch (int ex) {
				write_log("Memory for A12", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading A12 from the file", 3);
				throw 3;
			}

			sstepmat.matr_prod (A11f, A12, A_11t_12, r_a11, r_a12, c_a11, c_a12);

			mkl_free(A11f); // we do not need A11 anymore
			mkl_free(A12);


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n A11_i * A12: \n");
			for (auto i=0; i<min(r_a11,20); i++) {
				for (auto j=0; j<min(c_a12,20); j++) {
					printf ("%12.5G", A_11t_12[j+i*c_a12]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------


			runtimelog (true, "A21 * A11(-1)_A12 ...", 0, false, false);
			runtimelog (false, "Restore A21 from binary ...", 0, false, false);

			try {
				status = sstepmat.make_res_array (A21, r_a21, c_a21);
				if (status != 0) throw status;
				fa21.seekg(0, std::ios::beg);
				fa21.read(reinterpret_cast<char*>(A21), r_a21*c_a21*sizeof( double ));
			}
			catch (int ex) {
				write_log("Memory for A21", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading A21 from the file", 3);
				throw 3;
			}


			//----- Start DEBUG PRINT -------------------------------------------------
			std::cout<<std::endl;
			printf ("\n restored A21: \n");
			for (auto i=0; i<min(r_a21,20); i++) {
				for (auto j=0; j<min(c_a21,20); j++) {
					printf ("%12.5G", A21[j+i*c_a21]);
				}
				printf ("\n");
			}
			std::cout<<std::endl;
			//----- End DEBUG PRINT ---------------------------------------------------


			// A21_11t_12 = A21 * A_11t_12. Rows = columns = c_a12.
			double *A21_11t_12;
			status = sstepmat.make_res_array (A21_11t_12, c_a12, c_a12);
			if (status != 0) {
				write_log("Memory for A21_11t_12", status);
				throw status;
			}

			runtimelog (false, "Do matrix product ...", 0, false, false);

			sstepmat.matr_prod (A21, A_11t_12, A21_11t_12, r_a21, r_a11, c_a21, c_a12);


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of A21 * A11_i_A12: \n");
			for (auto i=0; i<min(c_a12,20); i++) {
				for (auto j=0; j<min(c_a12,20); j++) {
					printf ("%12.5G", A21_11t_12[j+i*c_a12]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT --------------------------------------------------


			mkl_free(A_11t_12); // we do not need A_11t_12 anymore
			mkl_free(A21); // we do not need A21 and A12 anymore

			runtimelog (true, "A22 - A21_A11(-1)_A12 ...", 0, false, false);
			runtimelog (false, "Restore A22 from binary ...", 0, false, false);

			try {
				status = sstepmat.make_res_array_h (A22, r_a22);
				if (status != 0) throw status;
				fa22.seekg(0, std::ios::beg);
				fa22.read(reinterpret_cast<char*>(A22), ((r_a22*r_a22 + r_a22)/2)*sizeof( double ));
				if (fa22.fail()) throw 4;
			}
			catch (int ex) {
				write_log("Memory for A22", ex);
				throw ex;
			}
			catch (std::exception const& e) {
				write_log (e.what(), 3);
				write_log ("Exception reading A22 from the file", 3);
				throw 3;
			}

			runtimelog (false, "Do matrix substitution ...", 0, false, false);

			// Get A22inv: A22 = A22 - A21_11t_12 (only 'L' part);
			for (size_t i = 1; i <= genotypedID.size(); i++) {
				for (size_t j = 1; j<= i; j++) {
					A22[i*(i-1)/2 + j-1] = A22[i*(i-1)/2 + j-1] - A21_11t_12[(i - 1)*genotypedID.size() + j - 1];
				}
			}

			mkl_free(A21_11t_12); // we do not need A21_11t_12 anymore


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of A22inv: \n");
			for (auto i=1; i<=min(genotypedID.size(),20); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", A22[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT --------------------------------------------------



			{

				 for debugging: A22_inv * A22 = I

				//----- Start DEBUG PRINT -------------------------------------------------
				std::cout<<"DEBUG: Restore A22 from the binary file ..."<<std::endl<<std::endl;
				//----- End DEBUG PRINT ---------------------------------------------------

				double *Agen;
				status = sstepmat.make_res_array_h (Agen, genotypedID.size());
				if (status != 0) {
					write_log("Memory for Agen", status);
					throw status;
				}

				size_t sz = (genotypedID.size()*genotypedID.size() + genotypedID.size())/2;
				try {
					float *a; a = (float *)malloc( sizeof( float ));
					if (a == NULL) throw 40;
					fA22.seekg(0, std::ios::beg);
					for (auto i = 0; i < sz; i++) {
						a[0] = 0.0f;
						fA22.read(reinterpret_cast<char*> (a), sizeof( float ));
						Agen[i] = static_cast<double>(a[0]);
					}
					free(a);
				}
				catch (int ex) {
					write_log("Memory for a", ex);
					throw ex;
				}
				catch (std::exception const& e) {
					write_log (e.what(), 3);
					write_log ("Exception reading file of A22 (s.step part) matrix", 3);
					throw 3;
				}

				//----- Start DEBUG PRINT -------------------------------------------------
				std::cout<<std::endl;
				printf ("\n DEBUG: directly generated Agen: \n");
				for (auto i=1; i<=min(genotypedID.size(),20); i++) {
					for (auto j=1; j<=i; j++) {
						printf ("%12.5G", Agen[i*(i - 1)/2 + j - 1]);
					}
					printf ("\n");
				}

				double ma_el = 0.0;
				double mi_el = 0.0;
				size_t ma_num = 0;
				size_t mi_num = 0;
				for (size_t i = 0; i<genotypedID.size(); i++){
					for (size_t j = 0; j<=i; j++){
						if (i != j){
							if (Agen[i*(i + 1)/2 + j] != 0.0){
								if(Agen[i*(i + 1)/2 + j] < mi_el){
									mi_num++;
									mi_el = Agen[i*(i + 1)/2 + j];
								}
								if(Agen[i*(i + 1)/2 + j] > ma_el){
									ma_num++;
									ma_el = Agen[i*(i + 1)/2 + j];
								}
							}
						}
					}
				}
				std::cout<<"mi_num, mi_el:"<<mi_num<<", "<<mi_el<<std::endl;
				std::cout<<"ma_num, ma_el:"<<ma_num<<", "<<ma_el<<std::endl;
				std::cout<<std::endl;

				//----- End DEBUG PRINT ---------------------------------------------------

				// make Agen full to do multiplication further on
				double *Agenf;
				status = sstepmat.make_res_array (Agenf, genotypedID.size(), genotypedID.size());
				if (status != 0) {
					write_log("Memory for Agenf", status);
					throw status;
				}

				double *A22f;
				status = sstepmat.make_res_array (A22f, genotypedID.size(), genotypedID.size());
				if (status != 0) {
					write_log("Memory for A22f", status);
					throw status;
				}

//#pragma omp parallel for
				for (size_t i = 0; i < genotypedID.size(); i++) {
					for (size_t j = 0; j <= i; j++) {
						Agenf[(i)*genotypedID.size() + j] = Agenf[(j)*genotypedID.size() + i] = Agen[i*(i + 1)/2 + j];
						A22f[(i)*genotypedID.size() + j] = A22f[(j)*genotypedID.size() + i] = A22[i*(i + 1)/2 + j];
					}
				}

				//----- Start DEBUG PRINT -------------------------------------------------
				printf ("\n Top left corner of Agenf: \n");
				for (auto i=0; i<min(genotypedID.size(),20); i++) {
					for (auto j=0; j<min(genotypedID.size(),20); j++) {
						printf ("%12.5G", Agenf[i*genotypedID.size() + j]);
					}
					printf ("\n");
				}
				//----- End DEBUG PRINT --------------------------------------------------

				//----- Start DEBUG PRINT -------------------------------------------------
				printf ("\n Top left corner of A22f: \n");
				for (auto i=0; i<min(genotypedID.size(),20); i++) {
					for (auto j=0; j<min(genotypedID.size(),20); j++) {
						printf ("%12.5G", A22f[i*genotypedID.size() + j]);
					}
					printf ("\n");
				}
				//----- End DEBUG PRINT --------------------------------------------------

				//mkl_free(Agen); // we do not need A11 anymore

				double *I;
				status = sstepmat.make_res_array (I, genotypedID.size(), genotypedID.size());
				if (status != 0) {
					write_log("Memory for I", status);
					throw status;
				}

				size_t dim = genotypedID.size();
				//cblas_dgemm (CblasRowMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, 1.0, Agenf, dim, A22f, dim, 0.0, I, dim);
				sstepmat.matr_prod(A22f, Agenf, I, genotypedID.size(), genotypedID.size(), genotypedID.size(), genotypedID.size());

				double big = 0.00000001;
				double small = -0.0000001;

				size_t ones = 0;
				size_t mins = 0;
				double min_e = small;
				double max_e = big;
				size_t maxes = 0;
				for (size_t i = 0; i < genotypedID.size()*genotypedID.size(); i++) {
					if (I[i] != 0.0){
						if (I[i] >= 1.0+small && I[i] <= 1.0+big) {
							ones++;
						}
						else if(I[i] < small && I[i] < min_e){
							mins++;
							min_e = I[i];
						}
						else if(I[i] > big && I[i] > max_e){
							max_e = I[i];
						}
						else if (I[i] > big) maxes++;
					}
				}

				std::cout<<"ones, genotyped: "<<ones<<", "<<genotypedID.size()<<std::endl;
				std::cout<<"mins: "<<mins<<std::endl;
				std::cout<<"maxes : "<<maxes<<std::endl;
				if (mins != 0) std::cout<<"min_e: "<<min_e<<std::endl;
				if (maxes != 0) std::cout<<"max_e: "<<max_e<<std::endl;

				//----- Start DEBUG PRINT -------------------------------------------------
				printf ("\n Top left corner of I: \n");
				for (auto i=0; i<min(genotypedID.size(),20); i++) {
					for (auto j=0; j<=i; j++) {
						printf ("%12.5G", I[i*genotypedID.size() + j]);
					}
					printf ("\n");
				}
				//----- End DEBUG PRINT --------------------------------------------------

			}


			runtimelog (true, "Write out A22(-1) to a file ...", 0, false, false);

			double *aVal;
			status = sstepmat.make_res_array(aVal, 1, 1);
			if (status != 0) {
				write_log("Memory for aVal", status);
				throw status;
			}

			// write out inversed A22
			for (size_t i = 0; i < (r_a22*r_a22 + r_a22)/2; i++) {
				try {
					aVal[0]  = A22[i];
					fa22i.write(reinterpret_cast<char*>(aVal), sizeof( double ));
				}
				catch (std::exception const& e) {
					write_log (e.what(), 3);
					write_log ("Exception writing A22i to the file", 3);
					throw 3;
				}
			}

			mkl_free(aVal);
			mkl_free(A22); // we do not need A22 anymore
			 */

			//double *A22inv; // temporal, it is not in use anymore; now A22 == A22inv !

			/*
			// for full store matrix
			// Get A22inv: A22inv = A22 - A21_11t_12;
			double *A22inv;
			sstepmat.make_res_array (A22inv, r_a22, c_a22);
			sstepmat.matr_sub(A22, A21_11t_12, A22inv, false, false, r_a22, c_a22);

			mkl_free(A22); // we do not need A22 anymore
			mkl_free(A21_11t_12); // we do not need A21_11t_12 anymore

			printf ("\n Top left corner of A22inv: \n");
			for (auto i=0; i<min(r_a22,19); i++) {
				for (auto j=0; j<min(c_a22,19); j++) {
					printf ("%12.5G", A22inv[j+i*c_a22]);
				}
				printf ("\n");
			}
			 */


			// Rearrange A22, which is [genotypedID(sorted)], to A22 [coreID(sorted), noncoreID(sorted)]

			// A22 = P*A22*P', where P is a permutation matrix
			// 1) make P
			// 2) A22_1 = P*A22
			// 3) make P'
			// 4) A22_2 = A22_1*P'

			/*
			double *P, *Pt;
			sstepmat.make_res_array(P, r_ginv, c_ginv);
			sstepmat.make_res_array(Pt, r_ginv, c_ginv);
			 */

			runtimelog (true, "Make permutation map ...", 0, false, false);

			std::vector<int> a22list;

			a22list.reserve( coreID.size() + noncoreID.size() ); // preallocate memory
			a22list.insert( a22list.end(), coreID.begin(), coreID.end() );
			a22list.insert( a22list.end(), noncoreID.begin(), noncoreID.end() );

			std::map<int, int> p_map;

			//std::vector<int> p_vect; // permutation vector
			auto i_p = 1;
			for (auto const& elem : a22list ) {
				p_map[elem] = i_p;
				//p_vect.push_back(find_invect2(genotypedID, elem));
				i_p++;
			}

			/*
			i_p = 1;
			for (auto const& elem : p_vect ) {
				size_t ind = (i_p - 1) * r_ginv + elem;
				P[ind] = 1.0;
				i_p++;
			}
			 */

			/* Read file with A22_inv and do permutation to make write A22_inv in memory */

			runtimelog (true, "Read file with A22(-1) and do permutation to make write A22_i in memory ...", 0, false, false);

			status = sstepmat.make_res_array_h (A22, r_a22);
			if (status != 0) {
				write_log("Memory for A22", status);
				throw status;
			}

			status = sstepmat.make_res_array(aVal, 1, 1);
			if (status != 0) {
				write_log("Memory for aVal", status);
				throw status;
			}

			fa22i.seekg(0, std::ios::beg);

			for (size_t i = 0; i < genotypedID.size(); i++) { //(r_a22*r_a22 + r_a22)/2
				for (size_t j = 0; j <= i; j++) {
					size_t ind;
					size_t row = p_map[genotypedID[i]];
					size_t col = p_map[genotypedID[j]];
					ind = row*(row-1)/2 + col - 1;
					if (col > row) ind = col*(col-1)/2 + row - 1;
					try {
						fa22i.read(reinterpret_cast<char*>(aVal), sizeof( double ));
						A22[ind] = aVal[0];
						//std::cout<<"   "<<ind<<",  "<<aVal[0]<<", "<<A22[ind]<<", row, col "<<row<<", "<<col<<std::endl;
					}
					catch (std::exception const& e) {
						write_log (e.what(), 3);
						write_log ("... when reading A22i", 3);
						throw 3;
					}
				}
			}

			fa22i.close(); remove(filename14.c_str());

			mkl_free(aVal);

			/* ----------------- DEBUG -------------- */
			/*
			status = sstepmat.make_res_array(aVal, 1, 1);
			if (status != 0) {
				write_log("Memory for aVal", status);
				throw status;
			}
			for (size_t i = 0; i < (r_a22*r_a22 + r_a22)/2; i++) {
				try {
					aVal[0] = A22[i];
					fa22i_prm.write(reinterpret_cast<char*>(aVal), sizeof( double ));
				}
				catch (std::exception const& e) {
					write_log (e.what(), 3);
					write_log ("Exception writing A22i_prm to the file", 3);
					throw 3;
				}
			}
			mkl_free(aVal);
			 */
			/* --------- END DEBUG -------------- */


			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n permuted A22_inv: \n");
			for (auto i=1; i<=min(genotypedID.size(),29); i++) {
				for (auto j=1; j<=i; j++) {
					printf ("%12.5G", A22[i*(i - 1)/2 + j - 1]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */

			/*
			double * A22_1;
			sstepmat.make_res_array(A22_1, r_a22, c_a22);
			sstepmat.matr_prod(P, A22inv, A22_1, r_a22, r_a22, c_a22, c_a22);


			mkl_free(A22inv);


			mkl_dimatcopy ('r', 't', r_ginv, r_ginv, 1.0, P, c_ginv, r_ginv); // transpose P

			double * A22_2;
			sstepmat.make_res_array(A22_2, r_a22, c_a22);
			sstepmat.matr_prod(A22_1, P, A22_2, r_a22, r_a22, c_a22, c_a22);


			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n Top left corner of A22_1: \n");
			for (auto i=0; i<min(r_ginv,26); i++) {
				for (auto j=0; j<min(c_ginv,26); j++) {
					printf ("%12.5G", A22_1[j+i*c_ginv]);
				}
				printf ("\n");
			}
			printf ("\n Top left corner of A22_2: \n");
			for (auto i=0; i<min(r_ginv,26); i++) {
				for (auto j=0; j<min(c_ginv,26); j++) {
					printf ("%12.5G", A22_2[j+i*c_ginv]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------


			mkl_free(P);
			mkl_free(A22_1);
			 */

			runtimelog (true, "Make H_sub ...", 0, false, false);

			// make Hsub = Ginv - Ainv.
			// Here we modify A22, so Hsub == A22.

			/* Variables to calculate sparsity of Ginv and Hsub */

			size_t allVals = 0;
			size_t zerrG = 0;
			size_t zerrH = 0;

			gVal = (double *)malloc( (coreID.size()+noncoreID.size())*sizeof( double ));
			if (gVal == NULL) throw 40;

			fGi.seekg(0, std::ios::beg);
			/* Read from a file core x non-core half of G_inv matrix */
			for (size_t i = 0; i < (coreID.size()+noncoreID.size()); i++) {
				try {
					fGi.read(reinterpret_cast<char*> (gVal), (coreID.size()+noncoreID.size())*sizeof( double ));
				}
				catch (std::exception const& e) {
					write_log (e.what(), 3);
					write_log ("... when reading Gi", 3);
					throw 3;
				}
				for (size_t m = 0; m <= i /*(coreID.size()+noncoreID.size())*/; m++) {
					//A22_2[i*genotypedID.size()+m] = gVal[m] - A22_2[i*genotypedID.size()+m];

					A22[i*(i+1)/2 + m] = gVal[m] - A22[i*(i+1)/2 + m];

					allVals++;
					if ( gVal[m] == 0.0 ) zerrG++;
					if ( A22[i*(i+1)/2 + m] == 0.0 ) zerrH++;
				}
			}

			free(gVal);

			fGi.close(); remove(filename10.c_str());

			runtimelog (false, "Elements in G_inv which are zero, % ", static_cast<double>(100*zerrG/allVals), false, false);
			runtimelog (false, "Elements in [G - A22] which are zero, % ", static_cast<double>(100*zerrH/allVals), false, false);

			runtimelog (true, "Write H_sub to a file ...", 0, false, false);


			if (isHtxt) {
				status = write_toFile_h ("H_sub.dat", A22, a22list, false);
				if (status != 0) throw status;
			}

			status = write_toFile_h ("H_sub.bin", A22, a22list, true);
			if (status != 0) throw status;

			/*
			//----- Start DEBUG PRINT -------------------------------------------------
			printf ("\n G-A22_inv: \n");
			for (auto i=0; i<min(r_ginv,26); i++) {
				for (auto j=0; j<=i; j++) {
					printf ("%12.5G", A22[j+i*(i+1)/2]);
				}
				printf ("\n");
			}
			//----- End DEBUG PRINT ---------------------------------------------------
			 */


			/* Optional part related to full H_inv*/

			if (isH) {

				runtimelog (true, "Make H ...", 0, false, false);

				// Restore ainv
				try {
					dKey k;
					fAfi.seekg(0, std::ios::beg);
					strAinv *a; a = (strAinv *)malloc( sizeof( strAinv ));
					if (a == NULL) throw 40;
					for (size_t i = 0; i < ainv_sz; i++) {
						fAfi.read(reinterpret_cast<char*> (a), sizeof( strAinv ));
						k.day = a->id_1;
						k.id = a->id_2;
						ainv[k] = a->val;
					}
					free(a);
				}
				catch (int ex) {
					write_log("Memory for a", ex);
					throw ex;
				}
				catch (std::exception const& e) {
					write_log (e.what(), 3);
					write_log ("... when reading ainv", 3);
					throw 3;
				}

				// Make full H_inv
				for (size_t i = 1; i <= a22list.size(); i++) {
					dKey k;
					k.day = a22list[i-1];
					int switchFlag = 0;
					for (size_t j = 1; j <= i; j++) {
						if (switchFlag) k.day = a22list[i-1];
						k.id = a22list[j-1];
						if (k.day < k.id) {
							k.day = a22list[j-1];
							k.id = a22list[i-1];
							switchFlag = 1;
						}
						//std::cout<<"row col: "<<k.day<<", "<<k.id<<", ainv: "<<ainv[k]<<", + A22 "<<A22[i*(i-1)/2 + j-1]<<", ind = "<<i*(i-1)/2 + j-1<<std::endl;
						ainv[k] = ainv[k] + A22[i*(i-1)/2 + j-1];
					}
				}

				mkl_free(A22);

				runtimelog (true, "Write H to a file ...", 0, false, false);

				status = write_toFile2 ("H_mat", ainv);
				if (status != 0) throw status;

				ainv.clear();

				runtimelog (true, "Completed.", 0, false, false);
			}
			else {
				mkl_free(A22);

				runtimelog (true, "Completed.", 0, false, false);
			}

			// clean also all vectors, or put this job to a class destructor !

		} // =============== End of APY block ==============================================================

		//fA22.close(); remove(filename1.c_str());
		fA22i.close(); remove(filename2.c_str());
		fAfi.close(); remove(filename3.c_str());
		//fAri.close(); remove(filename4.c_str());
		//fGcci.close(); remove(filename5.c_str());
		//fGnc.close(); remove(filename6.c_str());
		//fGnni.close(); remove(filename7.c_str());
		//fGcn.close(); remove(filename8.c_str());
		//fG12.close(); remove(filename9.c_str());
		//fGi.close(); remove(filename10.c_str());
		//fa11.close(); remove(filename11.c_str());
		//fa21.close(); remove(filename12.c_str());
		//fa22.close(); remove(filename13.c_str());
		//fa22i.close(); remove(filename14.c_str());
		//fa12.close(); remove(filename15.c_str());

		runtimelog (false, "Execution status", 0, false, true);

	}
	catch (int e_int) {
		write_log (where, e_int);
		return e_int;
	}
	catch (std::exception const& e) {
		write_log (e.what(), 1);
		write_log ("... handling in the final catch()", -1);
		return -1;
	}
	catch(...) {
		write_log ("Other exception handled in the final catch()", -2);
		return -2;
	}

	return 0;
}

//===============================================================================================================

int APY::write_toFile_h (std::string where, double *what, std::vector<int> id_list, bool isBinary) {

	std::string where_err("writing to a file");

	size_t lda = id_list.size();
	size_t THRESHOLD = lda*(lda+1)*4;
	std::string buffer;
	buffer.reserve(THRESHOLD);

	int num_writes = 0;

	std::ofstream fHbin;
	std::ofstream fHtxt;

	size_t wrAllElem = 0;
	size_t wrZerElem = 0;

	try {
		if (isBinary) {
			fHbin.open (where, fHbin.binary | fHbin.trunc | fHbin.out);
			fHbin.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
			if (!fHbin.is_open()) throw 3;

			strAinv *a; a = (strAinv *)malloc( sizeof( strAinv ));
			if (a == NULL) throw 40;

			for (size_t i = 1; i <= id_list.size(); i++) {
				for (size_t j = 1; j <= i; j++) {
					a->id_1 = id_list[i-1];
					a->id_2 = id_list[j-1];
					a->val = what[i*(i - 1) / 2 + j - 1];
					wrAllElem++;
					if ( std::abs(a->val) > tol_Hsub )
						fHbin.write(reinterpret_cast<char*> (a), sizeof( strAinv ));
					else
						wrZerElem++;
				}
			}
			free(a);
			fHbin.close();

		}
		else {
			fHtxt.open (where, fHtxt.out);
			fHtxt.exceptions ( std::ifstream::failbit | std::ifstream::badbit );
			if (!fHtxt.is_open()) throw 3;

			for (size_t i = 1; i <= id_list.size(); i++) {
				for (size_t j = 1; j <= i; j++) {
					size_t inx = i*(i - 1) / 2 + j - 1;
					//==================================================================
					wrAllElem++;
					if ( std::abs(what[inx]) > tol_Hsub ) {
						std::string s = std::to_string(id_list[i-1]);
						s.append(" ");
						s.append(std::to_string(id_list[j-1]));
						s.append(" ");
						s.append(std::to_string(what[inx]));

						if (buffer.length() + s.length() + 1 >= THRESHOLD)
						{
							fHtxt << buffer;
							buffer.resize(0);
							num_writes++;
						}
						buffer.append(s);
						buffer.append(1, '\n');
					}
					else
						wrZerElem++;
					//==================================================================
				}
			}
			fHtxt << buffer;
			fHtxt.close();
		}
	}
	catch (std::exception const& e) {
		if (fHtxt.eof()) {
			return 0;
		}
		else {
			write_log (where_err, 10);
			//std::cerr << "Exception opening/reading file: "<< e.what() << std::endl;
			return 10;
		}
	}
	catch (...)
	{
		write_log (where, 1);
		return 1;
	}

	runtimelog (false, "Elements in [G - A22] turned to zero, % ", static_cast<double>(100*wrZerElem/wrAllElem), false, false);
	runtimelog (false, "Zero threshold applied to H_sub ", tol_Hsub, false, false);

	return 0;
}

//===============================================================================================================

int APY::write_toFile2 (std::string where, std::map <dKey, double> &what) {

	std::string where_err("writing H to a file");

	std::ofstream ped;
	ped.exceptions ( std::ifstream::failbit | std::ifstream::badbit );

	size_t THRESHOLD = what.size()*8;
	std::string buffer;
	buffer.reserve(THRESHOLD);

	int num_writes = 0;

	try {
		size_t row, col;
		double val;
		ped.open(where, std::ofstream::out);
		for (auto const &i : what) {
			row = i.first.day;
			col = i.first.id;
			val = i.second;
			//==================================================================
			if ( std::abs(val) > tol_Hsub ) {
				std::string s = std::to_string(i.first.day);
				s.append(" ");
				s.append(std::to_string(i.first.id));
				s.append(" ");
				s.append(std::to_string(i.second));
				if (buffer.length() + s.length() + 1 >= THRESHOLD)
				{
					ped << buffer;
					buffer.resize(0);
					num_writes++;
				}
				buffer.append(s);
				buffer.append(1, '\n');
			}

			//==================================================================
			//ped << row << " " << col << " " << val << std::endl;
		}
		ped << buffer;
		ped.close();
	}
	catch (std::exception const& e) {
		if (ped.eof()) {
			return 0;
		}
		else {
			write_log (where_err, 10);
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

int APY::check_ID () {

	std::string where("check_ID");
	int status1 = 0;
	int status2 = 0;

	try
	{

#pragma omp parallel sections default(shared) num_threads(2)
		{
#pragma omp section
			{
				status1 = check_gmatID ();
			}
#pragma omp section
			{
				status2 = check_gtypedID ();
			}
		}
		gmatID.clear(); // we do not deed it anymore

	}
	catch (...)
	{
		write_log (where, 1);
		return 1;
	}

	return std::max(status1, status2);
}

//===============================================================================================================
