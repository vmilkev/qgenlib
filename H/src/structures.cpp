/*
 * structures.cpp
 *
 *  Created on: Sep 7, 2017
 *      Author: vimi
 */


#include "structures.hpp"

bool operator<(const dKey& lhs, const dKey& rhs)
{
	if (lhs.day != rhs.day)
		return lhs.day < rhs.day;
	else
		return lhs.id < rhs.id;
}
