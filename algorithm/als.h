#pragma once

#include "../ttTensor.h"

namespace xerus {

    class ALSVariant {
    public:
        uint sites;
        value_t minimumLocalResidual;
        bool printProgress;
        // TODO printProgressToFile;
        // TODO bool useEnergyFunctional;
        // TODO std::function endCriterion
        //TODO all arguments should be tensorNetworks to allow efficient CG in DMRG case
        std::function<void(const TensorNetwork &, Tensor &, const Tensor &)> localSolver;
        
        static void lapack_solver(const TensorNetwork &_A, Tensor &_x, const Tensor &_b);
        
        ALSVariant(uint _sites, value_t _minimumLocalResidual, std::function<void(const TensorNetwork &, Tensor &, const Tensor &)> _localSolver) 
                : sites(_sites), minimumLocalResidual(_minimumLocalResidual), printProgress(false), localSolver(_localSolver) {
            REQUIRE(_sites>0, "");
            REQUIRE(_minimumLocalResidual>=0, "");
        }
        
        double operator()(const TTOperator &_A, TTTensor &_x, const TTTensor &_b, value_t _convergenceEpsilon = 1e-6,  std::vector<value_t> *_perfData = nullptr) const;
    };

    const ALSVariant ALS(1, 1e-13, ALSVariant::lapack_solver);
    const ALSVariant DMRG(2, 1e-13, ALSVariant::lapack_solver);

}

