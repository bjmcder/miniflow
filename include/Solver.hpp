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
struct SolverSettings{
    T upwind_factor;
    T relax_factor;
    T tolerance;
    int max_iters;

    /**
     * 
    */
    SolverSettings(){}

    /**
     * 
    */
    SolverSettings(const T& gamma,
                   const T& omega,
                   const T& tol,
                   const int& max_it){
        
        upwind_factor = gamma;
        relax_factor = omega;
        tolerance = tol;
        max_iters = max_it;
    }
};

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
        SolverSettings<T> _settings;

    public:

        /**
         * 
        */
        Solver(){}

        /**
         * 
        */
        Solver(Problem<T>& problem, SolverSettings<T>& settings){

            _problem = problem;
            _settings = settings;

            std::cout << "\tConstructing solution\n";
            _solution = Solution<T>(_problem);
        }

        T avg(T item1, T item2){
            return 0.5*(item1+item2);
        }

        /**
         * 
        */
        T derivative(int axis1, int axis2=0, int i=0, int j=0, int k=0){

            auto dh = _problem.geometry().cell_sizes();
            auto gamma = _settings.upwind_factor;

            T central = 0.0;
            T upwind = 0.0;

            auto& S1 = _solution.tentative_momentum(axis1);
            auto& S2 = _solution.tentative_momentum(axis2);

            auto base_idx = std::vector<size_t>({(size_t)i,(size_t)j,(size_t)k});

            auto fwd_idx1 = base_idx;
            auto fwd_idx2 = base_idx;

            fwd_idx1[axis2] += 1;
            fwd_idx2[axis1] += 1;

            T central_fwd = avg(S1(base_idx), S1(fwd_idx1)); 
            central_fwd *= avg(S2(base_idx), S2(fwd_idx2));

            auto bwd_idx1 = base_idx;
            auto bwdfwd_idx1 = base_idx;
            auto bwd_idx2 = base_idx;

            bwd_idx1[axis1] -= 1;
            bwdfwd_idx1[axis1] -= 1;
            bwdfwd_idx1[axis2] += 1;
            bwd_idx2[axis1] -= 1;

            T central_bwd = avg(S1(bwd_idx1), S1(bwdfwd_idx1));
            central_bwd *= avg(S2(bwd_idx2), S2(base_idx));
        
            central = (1/dh[axis2])*(central_fwd - central_bwd);
            
            T upwind_fwd = fabs(avg(S1(base_idx), S1(fwd_idx1)));
            upwind_fwd *= avg(S2(base_idx), -S2(fwd_idx2));

            T upwind_bwd = fabs(avg(S1(bwd_idx1), S1(bwdfwd_idx1)));
            upwind_bwd *= avg(S2(bwd_idx2), -S2(base_idx));

            upwind = (1/dh[axis2])*(upwind_fwd - upwind_bwd);

            return central + gamma*upwind;

        }

        /**
         * 
        */
        T advection(int axis, int i, int j=0, int k=0){

            T adv = 0.0;
            T du2dx(0), duvdy(0), duwdz(0);
            T duvdx(0), dv2dy(0), dvwdz(0);
            T duwdx(0), dvwdy(0), dw2dz(0);
            switch(axis){
                case 0:
                    du2dx = derivative(axis, axis, i, j, k);
                    duvdy = derivative(axis, 1, i, j, k);
                    duwdz = derivative(axis, 2, i, j, k);
                    adv = du2dx - duvdy - duwdz;
                    break;
                case 1:
                    duvdx = derivative(axis, 0, i, j, k);
                    dv2dy = derivative(axis, axis, i, j, k);
                    dvwdz = derivative(axis, 2, i, j, k);
                    adv = duvdx - dv2dy - dvwdz;
                    break;
                case 2:
                    duwdx = derivative(axis, 0, i, j, k); 
                    dvwdy = derivative(axis, 1, i, j, k);
                    dw2dz = derivative(axis, axis, i,j,k);
                    adv = duwdx - dvwdy - dw2dz;
                    break;
                default:
                    throw std::invalid_argument("Invalid dimension.");

            }

            return adv;
        }

        /**
         * 
        */
        T diffusion(int axis, int i, int j=0, int k=0){

            auto dim = _problem.geometry().dimension();
            auto dh = _problem.geometry().cell_sizes();

            auto& UVW = _solution.velocity_component(axis);

            T ans = 0.0;            
            T num = 0.0;

            for(int d=0; d<dim; d++){
                switch(d){
                    case 0:
                        num = (UVW(i+1,j,k) + 2*UVW(i,j,k) + UVW(i-1,j,k));
                        break;
                    case 1:
                        num = (UVW(i,j+1,k) + 2*UVW(i,j,k) + UVW(i,j-1,k));
                        break;
                    case 2:
                        num = (UVW(i,j,k+1) + 2*UVW(i,j,k) + UVW(i,j,k-1));
                        break;
                    default:
                        break;
                }
                ans += num/(dh[d]*dh[d]);
            }
            return ans;
       }

        /**
         * 
        */
        void compute_momenta(){

            auto Re = _problem.flow_parameters().Re();
            auto dt = _problem.timestepper().dt();

            auto dim = _problem.geometry().dimension();

            auto imax = _problem.geometry().ncells()[0]-1;
            auto jmax = (dim == 2) ? _problem.geometry().ncells()[1]-1 : 0;
            auto kmax = (dim == 3) ? _problem.geometry().ncells()[2]-1 : 0;

            for(int d=0; d<dim; d++){
                auto& UVW = _solution.velocity_component(d);
                auto& FGH = _solution.tentative_momentum(d);

                auto& body_forces = _problem.flow_parameters().body_forces();

                for (int k=1; k<kmax; k++){
                    for(int j=1; j<jmax; j++){
                        for(int i=1; i<imax; i++){
                            FGH(i,j,k) += UVW(i,j,k);
                            auto diff = diffusion(d,i,j,k);
                            auto advect = advection(d,i,j,k);
                            FGH(i,j,k) += dt*(Re*diff - advect + body_forces[d]);
                        }
                    }
                }
            }
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
            std::cout << "Applying boundary motions\n";
            _problem.boundaries().apply_motions(_solution);

            // Compute momenta
            compute_momenta();
            _solution.F.print_field2d(0);
            _solution.G.print_field2d(0);

            // Construct the right-hand side of the pressure equation

        }

        /**
         * 
        */
        void solve(){

            auto tmax = _problem.timestepper().max_time();
            int step_count = 0;
            do{
                step_count++;
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
                //if(step_count % write_every == 0){
                //  _solution.to_vtk(fname);
                //}
                throw;
                
            } while (_problem.timestepper().current_time() <= tmax);
        }

};

#endif