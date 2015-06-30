// Xerus - A General Purpose Tensor Library
// Copyright (C) 2014-2015 Benjamin Huber and Sebastian Wolf. 
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

#include "include/xerus.h"
#include <stdio.h>
#include <fstream>

void storeVeloData(xerus::FullTensor &_v, std::string _fname) {
	std::vector<size_t> oldDim = _v.dimensions;
	_v.reinterpret_dimensions(std::vector<size_t>({25,27,25,3}));
	std::ofstream outX(_fname+"_vx.dat");
	std::ofstream outY(_fname+"_vy.dat");
	std::ofstream outZ(_fname+"_vz.dat");
	for (size_t x=0; x<25; x+=1) {
		for (size_t y=0; y<27; y+=1) {
			outX << _v[{x,y,10,0}] << " ";
			outY << _v[{x,y,10,1}] << " ";
			outZ << _v[{x,y,10,2}] << " ";
		}
		outX << std::endl;
		outY << std::endl;
		outZ << std::endl;
	}
	outX.close();
	outY.close();
	outZ.close();
	_v.reinterpret_dimensions(oldDim);
}

void swap_endianness(size_t *n) {
	size_t res=0;
	while (*n != 0) {
		res <<= 8;
		res += (*n) | 0xFF;
		*n >>= 8;
	}
	*n = res;
}

//g++ -std=c++11 main.cxx -I /usr/include/suitesparse build/lib/libxerus.a ../lib/lapack/liblapacke.a ../lib/lapack/liblapack.a /usr/lib64/atlas/libf77blas.a /usr/lib64/atlas/libcblas.a /usr/lib64/atlas/libatlas.a -lcxsparse -lbfd -liberty -lz -ldl -lgfortran
int main() {
	std::ifstream in("data/fgv_013.bin", std::ios::binary);
	xerus::FullTensor velocity({3,600,242,600});
	xerus::Index i1,i2,i3,i4;
	in.read(reinterpret_cast<char*>(velocity.data.get()), 25*27*25*3*sizeof(double));
	in.close();
	velocity(i1,i2,i3,i4) = velocity(i4,i3,i2,i1);

// 	velocity(i1,i2,i3) = velocity(i1,i2,i3, 0);
	velocity.reinterpret_dimensions(std::vector<size_t>({3,2,5,2,5,2, 11,2,11, 3,2,5,2,5,2, 3}));
	
	storeVeloData(velocity, "channels/channel_full");
	
	xerus::TTTensor ttv(velocity, 1e-10, 1000);
	
	xerus::value_t velo_norm = xerus::frob_norm(velocity);
	
	std::ofstream out("channels/channel_ttapprox.dat");
	
	for (size_t r = 1000; r > 0; r -= 10) {
		std::cout << r << '\n' << std::flush;
		ttv.round(r);
		std::cout << r << " ALS" << '\r' << std::flush;
		
// 		xerus::decomposition_als(ttv, velocity, 1e-5, 2);
		/*
		std::vector<double> perf;
		xerus::ProjectionALSVariant pALS(xerus::ProjectionALS);
		pALS.printProgress = true; pALS.preserveCorePosition = false;
		pALS(ttv, ttvOpt, 1e-4, &perf);*/
		xerus::FullTensor approx(ttv);
        std::cout << "Current residual: " << xerus::frob_norm(approx-velocity)/velo_norm << std::endl;
		out << r << " " 
		    << xerus::frob_norm(approx-velocity)/velo_norm << " " 
 			<< ttv.datasize() << std::endl;
		storeVeloData(approx, "channels/channel_r"+std::to_string(r));
		xerus::FullTensor diff = approx-velocity;
		for(size_t i = 0; i < diff.size; ++i) {
			diff[i] = std::abs(diff[i]); ///approx[i]
		}
		storeVeloData(diff, "channels/channel_diff_r"+std::to_string(r));
	}
	out.close();
	
    return 0;
}
