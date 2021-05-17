#ifndef __SOLVER_HPP
#define __SOLVER_HPP

#include <limits>

#include "BoundaryConditions.hpp"
#include "Problem.hpp"
#include "Solution.hpp"

template<typename T>
struct SolverStats{

    int itercount;
    T residual;
    T l2_norm;
    T linf_norm;

    SolverStats(){
        reset();
    }

    void reset(){
        itercount = 0;
        residual = std::numeric_limits<T>::max();
        l2_norm = std::numeric_limits<T>::max();
        linf_norm = std::numeric_limits<T>::max();
    }
};

template<typename T>
class Solver{

    private:

        int _max_iters;

        Problem<T> _problem;
        SolverStats<T> _stats;
        Solution<T> _solution;

    public:
        Solver(Problem<T>& problem){

            _problem = problem;
            _solution = Solution<real_t>(_problem);

        }

        void sor_iteration(){
            auto P = &_solution.pressure;
        }

        void solve_pressure(){

            _stats.reset();

            do{
                sor_iteration();
                //compute_norm();

            }while (_stats.itercount < _max_iters);
        }

        /**
         * 
        */
        void step(){

            // Apply boundary conditions

            // Compute momenta

            // Construct the right-hand side of the pressure equation

        }

        /**
         * 
        */
        void solve(){

            auto tmax = _problem._tstepper.max_time();

            do{
                // Advance the timestep
                auto dim = _problem.geometry().dimension();
                auto cell_sizes = _problem.geometry().cell_sizes();
                auto Re = _problem.flow_parameters().Re();
                auto max_vels = _solution.max_velocity_components(dim);

                _problem._tstepper.advance(cell_sizes, Re, max_vels);

                step();

                // Save Timestep Solution
                
            } while (_problem._tstepper.current_time() <= tmax);
            
        }

};

#endif