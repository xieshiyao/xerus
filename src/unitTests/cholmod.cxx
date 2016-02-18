
// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2016 Benjamin Huber and Sebastian Wolf. 
// 
// Xerus is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
// 
// Xerus is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Affero General Public License for more details.
// 
// You should have received a copy of the GNU Affero General Public License
// along with Xerus. If not, see <http://www.gnu.org/licenses/>.
//
// For further information on Xerus visit https://libXerus.org 
// or contact us at contact@libXerus.org.


#include <xerus.h>

#include "../../include/xerus/misc/test.h"
using namespace xerus;

void error_handler(int status, const char *file, int line, const char *message) {
	if (status < 0) {
		LOG(fatal, "CHOLMOD had a fatal error in " << file << ":" << line << " (status: " << status << ") msg: " << message);
	} else {
		LOG(cholmod_warning, "CHOLMOD warns in " << file << ":" << line << " (status: " << status << ") msg: " << message);
	}
}

UNIT_TEST2(Cholmod, basics) {
	auto c = new cholmod_common();
	cholmod_start(c);
	c->itype = CHOLMOD_LONG;
	c->dtype = CHOLMOD_DOUBLE;
	c->error_handler = &error_handler;
	auto x = cholmod_allocate_sparse(10, 10, 10, 1, 1, 0, CHOLMOD_REAL, c);
	
	REQUIRE(x, "Failed");
}});
