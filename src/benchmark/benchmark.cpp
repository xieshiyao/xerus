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

#include <fstream>
#include <string>
#include <iostream>
#include <stdio.h>

#include "xerus.h"
#include <TimeMeasure.h>


std::string exec(const std::string cmd) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "ERROR";
    char buffer[128];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 128, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;
}

const size_t maxN = 2500; // Maximal single dimension
const size_t createN = 2*maxN*maxN; /// We want to have two times NxN random data
const size_t maxOffset = maxN*maxN-1; /// One NxN is for offset

std::mt19937_64 rnd;
std::normal_distribution<double> normalDist(0,1000);
std::uniform_int_distribution<size_t> offsetDist(0, maxOffset);

double* offset(double* const _p) {
    return _p+offsetDist(rnd);
}

std::map<std::string, std::map<std::vector<size_t>, size_t>> results;

void add_call(const std::string& _callName, const std::vector<size_t>& _parameters, const size_t _time) {
        results[_callName][_parameters] += _time;
}


int main() {
    LOG(benchmark, "Initializing...");
    
    // Find out what system we are on
    std::string cpuModel = exec( "cat /proc/cpuinfo | grep -m 1 'model name'" );
    cpuModel = explode(cpuModel, ':')[1];
    replace(cpuModel, "Intel(R)", "");
    replace(cpuModel, "AMD", "");
    replace(cpuModel, "Core(TM)", "");
    replace(cpuModel, "(tm)", "");
    replace(cpuModel, "CPU", "");
    replace(cpuModel, "@", "");
    reduce(cpuModel, " \t", "_");
    LOG(benchmark, "CPU model deduced as " << cpuModel);
    
    std::map<size_t, std::string> blasImplementations;
    blasImplementations.insert(std::pair<size_t, std::string>(1, "Reference_Blas"));
    blasImplementations.insert(std::pair<size_t, std::string>(2, "Atlas_serial"));
    blasImplementations.insert(std::pair<size_t, std::string>(3, "Atlas_parallel"));
    blasImplementations.insert(std::pair<size_t, std::string>(4, "Openblas_serial"));
    blasImplementations.insert(std::pair<size_t, std::string>(5, "Openblas_parallel"));
    blasImplementations.insert(std::pair<size_t, std::string>(6, "Other"));
    
    LOG(benchmark, "Please select the BLAS implementation used:");
    for(const std::pair<size_t, std::string>& impl : blasImplementations) {
        LOG(benchmark, "("<<impl.first<<") " << impl.second);
    }
    
    // Init random data while waiting for user input
    double* A = new double[createN];
    double* B = new double[createN];
    
    for(size_t i = 0; i < createN; ++i) {
        A[i] = normalDist(rnd);
        B[i] = normalDist(rnd);
    }
    
    double* const C = new double[createN];
    double* const D = new double[createN];
    double* const E = new double[createN];
    
    size_t blasImplementation;
    std::cin >> blasImplementation;
    
    LOG(benchmark, "Selected " << blasImplementations[blasImplementation]);
    
    std::string directory = std::string("benchmark/")+cpuModel+"/"+blasImplementations[blasImplementation];
    
    exec("mkdir -p "+directory);
    
    TimeMeasure timer;
    
    LOG(benchmark, "Finished Initializing. Starting first Round...");
    
    for(size_t i = 1; i <= 100; ++i) {
        for(size_t n = 10; n < maxN; n = ((n*19)/18)+1) {
            LOG(benchmark, "Starting n="<<n);
            const size_t nSqr = n*n;
            
            //-------------------------------- MEMSET -----------------------------------
            timer.step();
            array_set_zero(offset(A), nSqr);
            add_call("Set zero", {n}, timer.get());
            
            //-------------------------------- COPY -----------------------------------
            //Memcpy
            timer.step();
            array_copy(offset(A), offset(C), nSqr);
            add_call("Copy", {n}, timer.get());
            
            //-------------------------------- SCALE -----------------------------------
            timer.step();
            array_scale(offset(A), normalDist(rnd), nSqr);
            add_call("Scale", {n}, timer.get());
            
            //-------------------------------- ADD -----------------------------------
            timer.step();
            array_add(offset(A), normalDist(rnd), offset(B), nSqr);
            add_call("Add", {n}, timer.get());
            
            timer.step();
            array_scale_add(normalDist(rnd), offset(A), normalDist(rnd), offset(B), nSqr);
            add_call("Add Scaled", {n}, timer.get());
            
            //-------------------------------- 2-Norm PRODUCT -----------------------------------
            timer.step();
            blasWrapper::two_norm(offset(A), nSqr);
            add_call("Two Norm", {n}, timer.get());
            
            //-------------------------------- DOT PRODUCT ----------------------------------- 
            timer.step();
            blasWrapper::dot_product(offset(A), nSqr, offset(B));
            add_call("Dot Product", {n}, timer.get());
            
            //-------------------------------- matrix_vector_product -----------------------------------
            timer.step();
            blasWrapper::matrix_vector_product(offset(C), n, normalDist(rnd), offset(A), n, false, offset(B));
            add_call("Matrix vector Product", {n}, timer.get());
            
            //-------------------------------- dyadic_vector_product -----------------------------------
                        
            timer.step();
            blasWrapper::dyadic_vector_product(offset(C), n, n, normalDist(rnd), offset(A), offset(B));
            add_call("Dyadic Vector Product", {n}, timer.get());
            
            
            
            //-------------------------------- Matrix Matrix product -----------------------------------
            timer.step();
            blasWrapper::matrix_matrix_product(offset(C), n, n, normalDist(rnd), offset(A), false, n, offset(B), false);
            add_call("Matrix Matrix Product", {n}, timer.get());
            
            //-------------------------------- QR -----------------------------------
            timer.step();
            blasWrapper::qr(offset(C), offset(D), offset(A), n, n);
            add_call("QR", {n}, timer.get());
            
            //-------------------------------- Solve -----------------------------------
            timer.step();
            blasWrapper::solve(offset(C), offset(A), n, offset(B));
            add_call("Solve", {n}, timer.get());
             
            //-------------------------------- Solve LS -----------------------------------
            timer.step();
            blasWrapper::solve_least_squares(offset(C), offset(A), n, n, offset(B));
            add_call("Solve Least Squares", {n}, timer.get());
            
            //-------------------------------- SVD -----------------------------------
            timer.step();
            blasWrapper::svd(offset(C), offset(D), offset(E), offset(A), n, n);
            add_call("SVD", {n}, timer.get());
        }
        
        LOG(benchmark, "Finished one Round, writing to files...");
        for(const std::pair<std::string, std::map<std::vector<size_t>, size_t>>& call : results) {
            std::fstream file(directory+"/"+call.first, std::fstream::out);
            for(const std::pair<std::vector<size_t>, size_t>& callEntry : call.second) {
                for(size_t par : callEntry.first) {
                    file << par << " ";
                }
                file << callEntry.second/i << std::endl;
            }
            file.close();
        }
        LOG(benchmark, "Finished writing to files. Starting next round...");
    }
        
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;
    delete[] E;
    return 0;
}
