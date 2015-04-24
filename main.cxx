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

#include <TimeMeasure.h>

#include <suiteSparse/cs.h>

#include "xerus.h"
#include "sparseMatrix.h"
#include <stdio.h>
#include <fstream>

std::mt19937_64 rnd;
std::normal_distribution<double> normalDist(0,1000);
std::uniform_real_distribution<double> avgOneDist(0.5, 2);
std::uniform_real_distribution<double> distZeroOne(0, 1);

std::map<std::string, std::map<std::vector<size_t>, size_t>> results;

void add_call(const std::string& _callName, const std::vector<size_t>& _parameters, size_t _time) {
        results[_callName][_parameters] += _time;
}

void storeVeloData(xerus::FullTensor &_v, std::string _fname) {
	std::vector<size_t> oldDim = _v.dimensions;
	_v.reinterpret_dimensions(std::vector<size_t>({2048,512,3}));
	std::ofstream outX(_fname+"_vx.dat");
	std::ofstream outY(_fname+"_vy.dat");
	std::ofstream outZ(_fname+"_vz.dat");
	for (size_t x=0; x<2048; x+=2) {
		for (size_t y=0; y<512; y+=2) {
			outX << _v[{x,y,0}] << " ";
			outY << _v[{x,y,1}] << " ";
			outZ << _v[{x,y,2}] << " ";
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

int main() {
	std::ifstream in("channel_velocity_data_xy.dat", std::ios::binary);
	xerus::FullTensor velocity({2048,512,3});
	/*
	 * 1 2048
	 * 2 1024
	 * 4  512
	 * 8  256
	 * 16 128
	 * 32  64
	 * 64  32
	 * 
	 * */
	float tmp;
	for (size_t i=0; i<2048*512*3; ++i) {
		in.read(reinterpret_cast<char*>(&tmp), 4);
		velocity.data.get()[i] = (double)tmp;
	}
	in.close();
	
	velocity.reinterpret_dimensions(std::vector<size_t>({32,64,16,32,3}));
	
	storeVeloData(velocity, "channel_full");
	
	xerus::TTTensor ttv(velocity);
	xerus::TTOperator I(xerus::TTOperator::construct_identity({32,64,16,32,3,32,64,16,32,3}));
	size_t r = std::max(ttv.ranks())-1;
	xerus::value_t velo_norm = xerus::frob_norm(velocity);
	std::ofstream out("channel_ttapprox.dat");
	xerus::TTTensor ttvOpt(ttv);
// 	while (2051*r + 512*r*r > 2048*512*3) r-=1;
	while (r*r*(64+16+32)+r*(32+3) > 2048*512*3) r-=1;
	for (r=10; r>0; --r) {
		std::cout << r << '\r' << std::flush;
		ttv.round(r);
		std::vector<double> perf;
		xerus::ALS(I, ttv, ttvOpt, 1e-4, &perf);
		xerus::FullTensor approx(ttv);
		out << r << " " 
		    << xerus::frob_norm(approx-velocity)/velo_norm << " " 
 			<< ttv.datasize() << std::endl;
		storeVeloData(approx, "channel_r"+to_string(r));
	}
	out.close();
	
    return 0;
}
