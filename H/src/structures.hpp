/*
 * structures.hpp
 *
 * Definition of structures and overloaded operator < used across multiple classes
 *
 *  Created on: Aug 10, 2017
 *      Author: vimi
 */
#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

#include <vector>

//#pragma once

/* structure to hold pair(sire, dam) */
	struct pedPair {
		int id_1; // sire
		int id_2; // dame
	};

	struct dKey {
		int day;
		int id;
	};

	struct idlist {
		std::vector<int> r_list;
		std::vector<int> c_list;
		std::vector<int> r_list_rec;
		std::vector<int> c_list_rec;
	};

	struct strAinv {
		int id_1;
		int id_2;
		double val;
	};

	bool operator<(const dKey& lhs, const dKey& rhs);

#endif