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

/**
 * @file
 * @brief Implementation of the ADF variants. 
 */

#include <xerus/algorithms/uqAdf.h>
 
#include <xerus/indexedTensorMoveable.h>

#include <xerus/misc/basicArraySupport.h>
#include <xerus/misc/simpleNumerics.h>
#include <xerus/misc/internal.h>


namespace xerus {
    
    class InternalSolver {
        const Tensor one = Tensor::ones({1});
        
        const size_t N;
        const size_t d;
        
        const std::vector<std::vector<Tensor>> positions;
        const std::vector<Tensor>& solutions;
        
        TTTensor& x;
        
        
		std::vector<std::vector<Tensor>> leftStack; // From corePosition 1 to d-2
		std::vector<std::vector<Tensor>> rightStack;  // From corePosition 1 to d-1
		
		std::vector<std::vector<Tensor>> leftIsStack;
		std::vector<std::vector<Tensor>> leftOughtStack;
		
        
        
    public:
        
        static Tensor randVar_to_position(const double _v, const size_t _polyDegree) {
            const std::vector<xerus::misc::Polynomial> stochasticBasis = xerus::misc::Polynomial::build_orthogonal_base(_polyDegree, [](const double){return 1.0;}, -1., 1.);
            
            Tensor p({stochasticBasis.size()});
            for (size_t i = 0; i < stochasticBasis.size(); ++i) {
                p[i] = stochasticBasis[i](_v);
            }
            
            return p;
        }
        
        static std::vector<std::vector<Tensor>> create_positions(const TTTensor& _x, const std::vector<std::vector<double>>& _randomVariables) {
            std::vector<std::vector<Tensor>> positions(_x.degree());
            
            for(size_t corePosition = 1; corePosition < _x.degree(); ++corePosition) {
                positions[corePosition].reserve(_randomVariables.size());
                for(size_t j = 0; j < _randomVariables.size(); ++j) {
                    positions[corePosition].push_back(randVar_to_position(_randomVariables[j][corePosition-1], _x.dimensions[corePosition]));
                }
            }
            
            return positions;
        }
        
        
        InternalSolver(TTTensor& _x, const std::vector<std::vector<double>>& _randomVariables, const std::vector<Tensor>& _solutions) : 
            N(_randomVariables.size()), 
            d(_x.degree()),
            positions(create_positions(_x, _randomVariables)),
            solutions(_solutions),
            x(_x),
            leftStack(d, std::vector<Tensor>(N)), 
            rightStack(d, std::vector<Tensor>(N)),
            leftIsStack(d, std::vector<Tensor>(N)), 
            leftOughtStack(d, std::vector<Tensor>(N))
            {
                REQUIRE(_randomVariables.size() == _solutions.size(), "ERROR");
                LOG(bug, "N=" << N << ", d=" << d);
        }
        
        
        void calc_left_stack(const size_t _corePosition) {
//             LOG(bug, "Calc left stack @ " << _corePosition);
            REQUIRE(_corePosition+1 < d, "Invalid corePosition");
            
			if(_corePosition == 0) {
				Tensor shuffledX = x.get_component(0);
				shuffledX.reinterpret_dimensions({x.dimensions[0], x.rank(0)});
				
				for(size_t j = 0; j < N; ++j) {
					contract(leftIsStack[_corePosition][j], shuffledX, true, shuffledX, false, 1); // TODO simple Identity!
					contract(leftOughtStack[_corePosition][j], solutions[j], shuffledX, 1);
				}
				
			} else { // _corePosition > 0
				const Tensor shuffledX = reshuffle(x.get_component(_corePosition), {1, 0, 2});
				
				Tensor measCmp, tmp;
				for(size_t j = 0; j < N; ++j) {
					contract(measCmp, positions[_corePosition][j], shuffledX, 1);
					
					
					contract(tmp, measCmp, true, leftIsStack[_corePosition-1][j], false,  1);
					contract(leftIsStack[_corePosition][j], tmp, false, measCmp, false, 1);
					contract(leftOughtStack[_corePosition][j], leftOughtStack[_corePosition-1][j], measCmp, 1);
					
					
					
					if(_corePosition > 1) {
						contract(leftStack[_corePosition][j], leftStack[_corePosition-1][j], measCmp, 1);
					} else { // _corePosition == 1
						leftStack[_corePosition][j] = measCmp;
					}
				}
			}
		}
		
        
        void calc_right_stack(const size_t _corePosition) {
//             LOG(bug, "Calc right stack @ " << _corePosition);
            REQUIRE(_corePosition > 0 && _corePosition < d, "Invalid corePosition");
            const Tensor shuffledX = reshuffle(x.get_component(_corePosition), {1, 0, 2});
            
            Tensor tmp;
            for(size_t j = 0; j < N; ++j) {
                contract(tmp, positions[_corePosition][j], shuffledX, 1);
                
                if(_corePosition < d-1) {
                    contract(rightStack[_corePosition][j], tmp, rightStack[_corePosition+1][j], 1);
                } else {
                    contract(rightStack[_corePosition][j], tmp, one, 1);
                }
            }
        }
        
        std::vector<Tensor> calc_residual(const size_t _corePosition) {
//             LOG(bug, "Calculate residual @ " << _corePosition);
            std::vector<Tensor> residual;
            
            Tensor tmp;
            if(_corePosition > 0) {
                const Tensor shuffledX = reshuffle(x.get_component(_corePosition), {1, 0, 2});
                
                for(size_t j = 0; j < N; ++j) {
                    // Current node
                    contract(tmp, positions[_corePosition][j], shuffledX, 1);
                    
                    // Right part
                    if(_corePosition < d-1) {
                        contract(tmp, tmp, rightStack[_corePosition+1][j], 1);
                    } else {
                        contract(tmp, tmp, one, 1);
                    }
                    
                    // Left Part
                    if(_corePosition > 1) {
                        contract(tmp, leftStack[_corePosition-1][j], tmp, 1);
                    }
                    
                    // Spacial Part
                    contract(tmp, x.get_component(0), tmp, 1);
                    contract(tmp, one, tmp, 1); // TODO weg
                    
                    residual.push_back(tmp - solutions[j]);
                }
            } else { // _corePosition == 0
                
                for(size_t j = 0; j < N; ++j) {
                    // Spacial Part
                    contract(tmp, x.get_component(0), rightStack[1][j], 1);
                    contract(tmp, one, tmp, 1); // TODO weg
                    
                    residual.push_back(tmp - solutions[j]);
                }
            }
            
            return residual;
        }
        
        
        Tensor calculate_delta(const std::vector<Tensor>& _residual, const size_t _corePosition) {
//             LOG(bug, "Calculate delta @ " << _corePosition);
            //calculate_projected_gradient
			Tensor delta(x.get_component(_corePosition).dimensions);
			Tensor delta2(x.get_component(_corePosition).dimensions);
            Tensor dyadComp;
            
// 			if(_corePosition > 0) {
// 				for(size_t j = 0; j < N; ++j) {
// 					// Calculate "direction" 
// 					if(_corePosition < d-1) {
// 						contract(dyadComp, positions[_corePosition][j], rightStack[_corePosition+1][j], 0);
// 					} else {
// 						dyadComp = positions[_corePosition][j];
// 						dyadComp.reinterpret_dimensions({dyadComp.dimensions[0], 1});
// 					}
// 					
// 					if(_corePosition > 1) {
// 						contract(dyadComp, leftStack[_corePosition-1][j], dyadComp, 0);
// 						contract(dyadComp, x.get_component(0), dyadComp, 1);
// 					} else {
// 						contract(dyadComp, x.get_component(0), dyadComp, 0);
// 					}
// 					
// 					contract(dyadComp, one, dyadComp, 1); // TODO weg
// 					
// 					// Scale with residual
// 					contract(dyadComp, _residual[j], dyadComp, 1);
// 					
// 					delta += dyadComp;
// 				}
// 			} else { // _corePosition == 0
// 				for(size_t j = 0; j < N; ++j) {
// 					contract(dyadComp, _residual[j], rightStack[_corePosition+1][j], 0);
// 					
// 					dyadComp.reinterpret_dimensions(delta.dimensions);
// 					delta += dyadComp;
// 				}
// 			}
            
//             LOG(bug, "Calculate alternative delta @ " << _corePosition);
			
			Tensor tmp;
			if(_corePosition > 0) {
				const Tensor shuffledX = reshuffle(x.get_component(_corePosition), {1, 0, 2});
				for(size_t j = 0; j < N; ++j) {
					// Calculate "is"
					contract(tmp, positions[_corePosition][j], shuffledX, 1);
					contract(tmp, tmp, true, leftIsStack[_corePosition-1][j], false, 1);
					if(_corePosition < d-1) {
						contract(tmp, rightStack[_corePosition+1][j], true, tmp, false, 1);
					} else {
						contract(tmp, one, true, tmp, false, 1);
					}
					contract(tmp, tmp, positions[_corePosition][j], 0);
					
					if(_corePosition < d-1) {
						contract(dyadComp, tmp, rightStack[_corePosition+1][j], 0);
					} else {
						dyadComp = tmp;
						dyadComp.reinterpret_dimensions({dyadComp.dimensions[0], dyadComp.dimensions[1], 1});
					}
					
					delta2 += dyadComp;
					
					// Calculate "Ought"
					contract(tmp, leftOughtStack[_corePosition-1][j], positions[_corePosition][j], 0);
					
					if(_corePosition < d-1) {
						contract(dyadComp, tmp, rightStack[_corePosition+1][j], 0);
					} else {
						dyadComp = tmp;
						dyadComp.reinterpret_dimensions({dyadComp.dimensions[0], dyadComp.dimensions[1], 1});
					}
					
					delta2 -= dyadComp;
					
				}
			} else { // _corePosition == 0
				for(size_t j = 0; j < N; ++j) {
					contract(dyadComp, x.get_component(0), rightStack[_corePosition+1][j], 1);
					dyadComp.reinterpret_dimensions({x.dimensions[0]});
					dyadComp -= solutions[j];
					
					dyadComp.reinterpret_dimensions({1, x.dimensions[0]});
					contract(dyadComp, dyadComp, rightStack[_corePosition+1][j], 0);
					delta2 += dyadComp;
				}
			}
            
//             LOG(bug, "Difference: " << frob_norm(delta - delta2));
            
//             LOG(norm, "Norm @ " << _corePosition << ": " << frob_norm(delta2));
            return delta2;
        }
        
        
        double calculate_norm_A_projGrad(const Tensor& _delta, const size_t _corePosition) {
            double norm = 0.0;
            
            if(_corePosition == 0) {
                Tensor tmp;
                
                for(size_t j = 0; j < N; ++j) {
                    
                    // Spacial Part
                    contract(tmp, _delta, rightStack[1][j], 1);
                    
                    norm += misc::sqr(frob_norm(tmp));
                }
            } else { // _corePosition > 0
                const Tensor shuffledX = reshuffle(_delta, {1, 0, 2});
                Tensor tmp;
                
                for(size_t j = 0; j < N; ++j) {
                    
                    // Current node
                    contract(tmp, positions[_corePosition][j], shuffledX, 1);
                    
                    // Right part
                    if(_corePosition < d-1) {
                        contract(tmp, tmp, rightStack[_corePosition+1][j], 1);
                    } else {
                        contract(tmp, tmp, one, 1);
                    }
                    
                    // Left Part
                    if(_corePosition > 1) {
                        contract(tmp, leftStack[_corePosition-1][j], tmp, 1);
                    }
                    
                    // Spacial Part
                    contract(tmp, x.get_component(0), tmp, 1);
                    
                    norm += misc::sqr(frob_norm(tmp));
                }
            }
            
            return std::sqrt(norm);
        }
        
        void update_x(const Tensor& _delta, const double _normAProjGrad, const size_t _corePosition) {
            const value_t PyR = misc::sqr(frob_norm(_delta));
            
            x.component(_corePosition) -= (PyR/misc::sqr(_normAProjGrad))*_delta;
        }
        
        
        double residual_norm(const std::vector<Tensor> _resiudal) {
            double norm = 0;
            for(const auto& r : _resiudal) {
                norm += misc::sqr(frob_norm(r));
            }
            
            return std::sqrt(norm);
        }
        
        double calc_solutions_norm() {
            double norm = 0;
            for(const auto& s : solutions) {
                norm += misc::sqr(frob_norm(s));
            }
            
            return std::sqrt(norm);
        }
        
        void solve() {
			
			const size_t maxIterations = 10;
			
			for(size_t iteration = 0; maxIterations == 0 || iteration < maxIterations; ++iteration) {
				x.require_correct_format();
				LOG(start, x.dimensions << " and " << x.ranks());
				for(size_t corePosition = 0; corePosition < x.degree(); ++corePosition) {
					REQUIRE(x.get_component(corePosition).all_entries_valid(), "Invalid entries at " << corePosition);
				}
				x.move_core(0, true);
				
				// Rebuild right stack
				for(size_t corePosition = d-1; corePosition > 0; --corePosition) {
					calc_right_stack(corePosition);
				}
				
				
				for(size_t corePosition = 0; corePosition < x.degree(); ++corePosition) {
					
					const auto residual = calc_residual(corePosition);
					
					if(corePosition == 0) {
						LOG(bla, "Residual norm: " << residual_norm(residual)/calc_solutions_norm());
					}
					
					const auto delta = calculate_delta(residual, corePosition);
					
					const auto normAProjGrad = calculate_norm_A_projGrad(delta, corePosition);
					
					update_x(delta, normAProjGrad, corePosition);
					
					// If we have not yet reached the end of the sweep we need to take care of the core and update our stacks
					if(corePosition+1 < d) {
						x.move_core(corePosition+1, true);
						calc_left_stack(corePosition);
					}
				}
			}
        }
    };
    
    
    
    void uq_adf(TTTensor& _x, const std::vector<std::vector<double>>& _randomVariables, const std::vector<Tensor>& _solutions) {
        InternalSolver solver(_x, _randomVariables, _solutions);
        return solver.solve();
    }
	
} // namespace xerus
