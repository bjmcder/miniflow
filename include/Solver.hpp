#ifndef __SOLVER_HPP
#define __SOLVER_HPP

#include <limits>

#include "Solution.hpp"
#include "Problem.hpp"

template<typename T>
struct Solution;

template<typename T>
struct Problem;

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

        /**
         * 
        */
        Solver(){}

        /**
         * 
        */
        Solver(Problem<T>& problem){

            _problem = problem;
            std::cout << "\tConstructing solution\n";
            _solution = Solution<real_t>(_problem);
        }

        /**
         * 
        */
        void compute_momenta(){

        }

        /**
         * 
        */
        void sor_iteration(){
            auto P = &_solution.pressure;
        }

        /**
         * 
        */
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
            std::cout << "Applying boundary conditions:\n";
            _problem.boundaries().apply_conditions(_solution);

            _solution.U.print_field2d(0);
            _solution.V.print_field2d(0);

            // Compute momenta

            // Construct the right-hand side of the pressure equation

        }

        /**
         * 
        */
        void solve(){

            auto tmax = _problem.timestepper().max_time();

            do{

                std::cout << "Taking step t = "
                          << _problem.timestepper().current_time()
                          << "\n";
                step();

                // Advance the timestep
                auto dim = _problem.geometry().dimension();
                auto cell_sizes = _problem.geometry().cell_sizes();
                auto Re = _problem.flow_parameters().Re();
                std::cout << "doing max velocity\n";
                auto max_vels = _solution.max_velocity_components(dim);

                std::cout << "advancing timestep\n";
                _problem.timestepper().advance(cell_sizes, Re, max_vels);

                // Save Timestep Solution
                throw;
                
            } while (_problem.timestepper().current_time() <= tmax);
        }

};

#endif