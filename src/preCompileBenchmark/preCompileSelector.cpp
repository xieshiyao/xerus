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

#include "../xerus.h"

#include "timeMeasure.h"

using namespace xerus; 

#ifdef FULL_SELECTION_
    const size_t createN = (1024*1024*1024/8)/2;
    const size_t maxN = createN-128*1024*1024/8; // At least 128 MB of random space 
#else
    const size_t createN = (1024*1024*1024/8)/4;
    const size_t maxN = createN-16*1024*1024/8; // At least 16 MB of random space 
#endif

  
const size_t maxNSquare = (size_t) sqrt((double) maxN);

std::mt19937_64 rnd;

std::map<std::string, std::string> functionDeclarations; // One per group
std::map<std::string, std::string> functionDefinitions; // One per candidate
std::map<std::string, std::map<size_t, std::map<std::string, size_t>>> results;
//        ^^Group^^           ^Param^          ^^candidate^^  ^Time^


double* offset(double* const _p, const size_t _maxSize) {
    std::uniform_int_distribution<size_t> offsetDist(0, createN-_maxSize);
    return _p+offsetDist(rnd);
}

_inline_ void add_call(const std::string& _groupName, const size_t _param, const std::string& _callName, const size_t _time) {
        results[_groupName][_param][_callName] += _time;
}

void determine_winner(const std::string& _groupName, std::fstream& _headerFile) {
    std::map<std::string, double> totalScores;
    
    for(const std::pair<size_t, std::map<std::string, size_t>>& param : results[_groupName]) {
        size_t totalTime = 0;
        for(const std::pair<std::string, size_t>& functionScore : param.second) {
            totalTime += functionScore.second;
        }
        
        totalTime = std::max(totalTime, 1ul); //Ensure we have at least 1
        
        for(const std::pair<std::string, size_t>& functionScore : param.second) {
            totalScores[functionScore.first] += double(functionScore.second)/double(totalTime);
        }
    }
    
    double bestScore = 1e100;
    std::string winner;
    LOG(Selection, "Results for group " << _groupName << ":");
    for(const std::pair<std::string, double>& candidate : totalScores) {
        LOG(Selection, "Candidate " <<candidate.first << " has a total score of " << candidate.second);
        if(candidate.second < bestScore) {
            bestScore = candidate.second;
            winner = candidate.first;
        }
    }
    LOG(Selection, "The winner in the group " << _groupName << " is the function: " << winner <<std::endl);
    
    _headerFile << functionDeclarations[_groupName] << " {" << std::endl;
    _headerFile << "\t" << functionDefinitions[winner] << std::endl;
    _headerFile << "}" << std::endl << std::endl;
}

 


_inline_ void self_set(double* const __restrict _x, const size_t _n) {
    for(size_t i=0; i<_n; ++i) { _x[i] = 0.0; }
}
   
_inline_ void self_copy(double* const __restrict _x, const double* const __restrict _y, const size_t _n) {
    for(size_t i=0; i<_n; ++i) { _x[i] = _y[i]; }
}

_inline_ void self_scaled_copy(double* const __restrict _x, const double _alpha, const double* const __restrict _y, const size_t _n) {
    for(size_t i=0; i<_n; ++i) { _x[i] = _alpha*_y[i]; }
}

_inline_ void self_scale(double* const __restrict _x, const size_t _n, const double _alpha) {
    for(size_t i = 0; i < _n; i++ ) {  _x[i] *= _alpha; }
}

_inline_ void self_add(double* const __restrict _x, const size_t _n, const double _alpha, const double* const __restrict _y) {
    for(size_t i = 0; i < _n; i++ ) { _x[i] += _alpha*_y[i]; }
}

_inline_ void self_scale_add(double* const __restrict _x, const size_t _n, const double _alpha, const double _beta, const double* const __restrict _y) {
    for(size_t i = 0; i < _n; i++ ) { _x[i] = _alpha*_x[i]+_beta*_y[i]; }
}


 

int main() {
    REQUIRE(createN > 128*1024*1024, " Arrays are too small");
    
    LOG(Selection, "Initializing...");
    TimeMeasure timer;
    size_t timeControl;
    #ifdef FULL_SELECTION_
        const size_t maxTime = 5*1000*1000;
    #else
        const size_t maxTime = 1*1000*1000;
    #endif
    
    
    std::fstream headerFile("misc/selectedFunctions.h", std::fstream::out);
    headerFile << "#pragma once" << std::endl << std::endl;
    headerFile << "#include <cstring>" << std::endl;
    headerFile << "#include \"blasLapackWrapper.h\"" << std::endl << std::endl << std::endl;
    
    headerFile << "START_MISC_NAMESPACE" << std::endl << std::endl;
    
    // Standard template functions for ordinary types
    headerFile << "template <typename T>" << std::endl;
    headerFile << "_inline_ void array_set_zero(T* const __restrict _x, const size_t _n) {" << std::endl;
    headerFile << "    memset(_x, 0, _n*sizeof(T));" << std::endl;
    headerFile << "}" << std::endl << std::endl;
        
    headerFile << "template <typename T>" << std::endl;
    headerFile << "_inline_ void array_copy(T* const __restrict _x, const T* const _y, const size_t _n) {" << std::endl;
    headerFile << "    memcpy(_x, _y, _n*sizeof(T));" << std::endl;
    headerFile << "}" << std::endl << std::endl;
    

    

    
    std::normal_distribution<double> normalDist(0, 1000);
    double* A = new double[createN];
    double* B = new double[createN];
    for(size_t i = 0; i < createN; ++i) {
        A[i] = normalDist(rnd);
        B[i] = normalDist(rnd);
    }
    
    LOG(benchmark, "Finished Initializing.");
    
     
    LOG(Selection, "Try to determinie best function to set memory equals zero.");
    
    functionDeclarations["set_zero"] = " _inline_ void array_set_zero(double* const __restrict _x, const size_t _n)";
    functionDefinitions["memset"] = "memset(_x, 0.0, _n*sizeof(double));";
    functionDefinitions["self_set"] = "for(size_t i=0; i<_n; ++i) { _x[i] = 0; }";
    functionDefinitions["atlas_set"] = "REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), \"Dimension to large for blas/lapack\"); \n\tcatlas_dset((int) _n, 0.0, _x, 1);";
        
    timeControl = uTime()+maxTime;
    while(uTime() < timeControl) {
        for(size_t m = 100; m < maxN; m = ((m+1)*27)/26) {
            timer.step();
            memset(offset(A, m), 0, m*sizeof(double));
            add_call("set_zero", m, "memset", timer.get());
            
            timer.step();
            self_set(offset(A, m), m);
            add_call("set_zero", m, "self_set", timer.get());
            
            #ifdef _ATLAS
                timer.step();
                REQUIRE(m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
                catlas_dset(m, 0.0, offset(A, m), 1);
                add_call("set_zero", m, "atlas_set", timer.get());
            #endif
        }
    }
    determine_winner("set_zero", headerFile);
    
    
    
    LOG(Selection, "Try to determinie best function to copy a vector to another.");
    
    functionDeclarations["copy"] = "_inline_ void array_copy(double* const __restrict _x, const double* const _y, const size_t _n)";
    functionDefinitions["memCpy"] = "memcpy(_x, _y, _n*sizeof(double));";
    functionDefinitions["selfCpy"] = "for(size_t i=0; i<_n; ++i) { _x[i] = _y[i]; }";
    functionDefinitions["blasCpy"] = "REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), \"Dimension to large for blas/lapack\"); \n\tcblas_dcopy((int) _n, _y, 1, _x, 1);";

    timeControl = uTime()+maxTime;
    while(uTime() < timeControl) {
        for(size_t m = 100; m < maxN; m = ((m+1)*27)/26) {
            timer.step();
            memcpy(offset(B, m), offset(A, m), m*sizeof(double));
            add_call("copy", m, "memCpy", timer.get());
                
            timer.step();
            self_copy(offset(B, m), offset(A, m), m);
            add_call("copy", m, "selfCpy", timer.get());
            
            timer.step();
            REQUIRE(m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
            cblas_dcopy((int) m, offset(A, m), 1, offset(B, m), 1);
            add_call("copy", m, "blasCpy", timer.get());
        }
    }
    determine_winner("copy", headerFile);
    
    
    
    
    LOG(Selection, "Try to determinie best function to scaled copy a vector to another.");
    
    functionDeclarations["scaledCopy"] = "_inline_ void array_scaled_copy(double* const __restrict _x, const double _alpha, const double* const _y, const size_t _n)";
    functionDefinitions["selfScaledCpy"] = "for(size_t i=0; i<_n; ++i) { _x[i] = _alpha*_y[i]; }";
    functionDefinitions["blasScaledCpy"] = "REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), \"Dimension to large for blas/lapack\"); \n\tcblas_dcopy((int) _n, _y, 1, _x, 1);\n\tcblas_dscal((int) _n, _alpha, _x, 1);";
    functionDefinitions["atlasScaledCpy"] = "REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), \"Dimension to large for blas/lapack\"); \n\tcatlas_daxpby((int) _n, _alpha, _y, 1, 0, _x, 1);";
    
    timeControl = uTime()+maxTime;
    while(uTime() < timeControl) {
        for(size_t m = 10; m < maxN; m = ((m+1)*27)/26) {
            timer.step();
            self_scaled_copy(offset(B, m), normalDist(rnd), offset(A, m), m);
            add_call("scaledCopy", m, "selfScaledCpy", timer.get());
            
            double* const x = offset(A, m);
            timer.step();
            REQUIRE(m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
            cblas_dcopy((int) m, x, 1, offset(B, m), 1);
            cblas_dscal((int) m, normalDist(rnd), x, 1);
            add_call("scaledCopy", m, "blasScaledCpy", timer.get());
            
            #ifdef _ATLAS
                timer.step();
                REQUIRE(m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
                catlas_daxpby((int) m, normalDist(rnd), offset(B, m), 1, 0, offset(A, m), 1);
                add_call("scaledCopy", m, "atlasScaledCpy", timer.get());
            #endif
        }
    }
    determine_winner("scaledCopy", headerFile);
    
    
    
    
    
    LOG(Selection, "Try to determinie best function to scale a vector by a constant factor.");
     
    functionDeclarations["scale"] = "_inline_ void array_scale(double* const __restrict _x, const double _alpha, const size_t _n)";
    functionDefinitions["selfScale"] = "for(size_t i = 0; i < _n; i++ ) {  _x[i] *= _alpha; }";
    functionDefinitions["blasScale"] = "REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), \"Dimension to large for blas/lapack\");\n\tcblas_dscal((int) _n, _alpha, _x, 1);";
    
    timeControl = uTime()+maxTime;
    while(uTime() < timeControl) {
        for(size_t m = 100; m < maxN; m = ((m+1)*27)/26) {
            timer.step();
            self_scale(offset(A, m), m, normalDist(rnd));
            add_call("scale", m, "selfScale", timer.get());
            
            timer.step();
            REQUIRE(m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
            cblas_dscal((int) m, normalDist(rnd), offset(A, m), 1);
            add_call("scale", m, "blasScale", timer.get());
        }
    }
    determine_winner("scale", headerFile);
    
    
    
     
    
    
    LOG(Selection, "Try to determinie best function to add a vector to another.");
    
    functionDeclarations["add"] = "_inline_ void array_add(double* const __restrict _x, const double _alpha, const double* const _y, const size_t _n)";
    functionDefinitions["selfAdd"] = "for(size_t i = 0; i < _n; i++ ) { _x[i] += _alpha*_y[i]; }";
    functionDefinitions["blasAdd"] = "REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), \"Dimension to large for blas/lapack\");\n\tcblas_daxpy((int) _n, _alpha, _y, 1, _x, 1);";
    functionDefinitions["AtlasAdd"] = "REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), \"Dimension to large for blas/lapack\");\n\tcatlas_daxpby((int) _n, _alpha, _y, 1, 1.0, _x, 1);";
    
    timeControl = uTime()+maxTime;
    while(uTime() < timeControl) {
        for(size_t m = 100; m < maxN; m = ((m+1)*27)/26) {
            timer.step();
            self_add(offset(A, m), m, normalDist(rnd), offset(B, m));
            add_call("add", m, "selfAdd", timer.get());
            
            timer.step();
            REQUIRE(m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
            cblas_daxpy((int) m, normalDist(rnd), offset(B, m), 1, offset(A, m), 1);
            add_call("add", m, "blasAdd", timer.get());
            
            #ifdef _ATLAS
                timer.step();
                REQUIRE(m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
                catlas_daxpby((int) m, normalDist(rnd), offset(B, m), 1, normalDist(rnd), offset(A, m), 1);
                add_call("add", m, "AtlasAdd", timer.get());
            #endif
        }
    }
    determine_winner("add", headerFile);
    
    
    
    
    
    LOG(Selection, "Try to determinie best function to  add a vector to a rescaled other.");
    

    functionDeclarations["scaleAdd"] = "_inline_ void array_scale_add(const double _alpha, double* const __restrict _x, const double _beta, const double* const _y, const size_t _n)";
    functionDefinitions["selfScaleAdd"] = "for(size_t i = 0; i < _n; i++ ) { _x[i] = _alpha*_x[i]+_beta*_y[i]; }";
    functionDefinitions["AtlasScaleAdd"] = "REQUIRE(_n <= (size_t) std::numeric_limits<int>::max(), \"Dimension to large for blas/lapack\");\n\tcatlas_daxpby((int) _n, _beta, _y, 1, _alpha, _x, 1);";
    
    timeControl = uTime()+maxTime;
    while(uTime() < timeControl) {
        for(size_t m = 100; m < maxN; m = ((m+1)*27)/26) {
            timer.step();
            self_scale_add(offset(A, m), m, normalDist(rnd), normalDist(rnd), offset(B, m));
            add_call("scaleAdd", m, "selfScaleAdd", timer.get());
            
            #ifdef _ATLAS
                timer.step();
                REQUIRE(m <= (size_t) std::numeric_limits<int>::max(), "Dimension to large for blas/lapack");
                catlas_daxpby((int) m, normalDist(rnd), offset(B, m), 1, normalDist(rnd), offset(A, m), 1);
                add_call("scaleAdd", m, "AtlasScaleAdd", timer.get());
            #endif
        }
    }
    determine_winner("scaleAdd", headerFile);
    
    
    
    headerFile << "END_MISC_NAMESPACE" << std::endl << std::endl;
    headerFile.close();
        
    delete[] A;
    delete[] B;
    return 0;
}
