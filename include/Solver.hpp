#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <limits>
#include <omp.h>

#include "Solution.hpp"
#include "Problem.hpp"
#include "OutputSettings.hpp"
#include "VTKOutput.hpp"

// Forward declarations
template<typename T>
struct Solution;

template<typename T>
struct Problem;

/**
 * Solver settings class. This is a container for all the settings that define
 * the behavior of the Navier-Stokes solver.
*/
template<typename T>
struct SolverSettings{
    T upwind_factor;
    T relax_factor;
    T tolerance;
    int max_iters;

    /**
     * Default constructor. Initialize the solver with preset values for
     * upwinding, relaxation, convergence tolerance and maximum iterations.
    */
    SolverSettings(): upwind_factor(0.5),
                      relax_factor(1.5),
                      tolerance(1e-3),
                      max_iters(100){}

    /**
     * Preferred constructor. Initialize the solver with user-defined settings.
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

/**
 * Solver statistics class. This is a container for tracking the state of the
 * Navier-Stokes solver with regards to iteration count, residuals and norms.
*/
template<typename T>
struct SolverStats{

    int itercount;
    T residual_sumsquares;
    T residual_max;
    T l2_norm;
    T linf_norm;

    /**************************************************************************
     * Default Constructor. Initialize all settings to zero.
    **************************************************************************/
    SolverStats(): itercount(0),
                   residual_sumsquares(0),
                   residual_max(std::numeric_limits<T>::max()),
                   l2_norm(std::numeric_limits<T>::max()),
                   linf_norm(std::numeric_limits<T>::max()){
        reset();
    }

    /**************************************************************************
     * Reset all solver statistics.
    **************************************************************************/
    void reset(){
        itercount = 0;
        reset_norm();
    }

    /**************************************************************************
     * Reset all solver norms.
    **************************************************************************/
    void reset_norm(){
        residual_sumsquares = 1.0;
        residual_max = 1.0;
        l2_norm = 1.0;
        linf_norm = 1.0;
    }
};

/******************************************************************************
 * Solver class. A collection of functions that operate on a problem
 * definition to solve the Navier-Stokes equations for vector and
 * scalar fields. Based on Ch. 3 of "Numerical Simulation in Fluid Dynamics"
 * by Griebel, et al.OUTPUT_HPP
******************************************************************************/
template<typename T>
class Solver{

    private:

        Problem<T> _problem;
        SolverStats<T> _stats;
        Solution<T> _solution;
        SolverSettings<T> _settings;
        OutputSettings _output_settings;

    public:

        /**********************************************************************
         * Default constructor (does nothing).
        **********************************************************************/
        Solver(){}

        /**********************************************************************
         * Preferred constructor. Uses problem definition and solver settings
         * based on user input.
        **********************************************************************/
        Solver(Problem<T>& problem,
               SolverSettings<T>& settings,
               OutputSettings& outsettings){

            _problem = problem;
            _settings = settings;
            _solution = Solution<T>(_problem);
            _output_settings = outsettings;
        }

        /**********************************************************************
         * Compute the average value of two numeric items.
        **********************************************************************/
        T avg(const T& item1, const T& item2){
            return 0.5*(item1+item2);
        }

        /**********************************************************************
         *
        **********************************************************************/
        T upwind_difference(int dir1, int dir2, int i, int j, int k){

            const auto& dh = _problem.geometry().cell_sizes();
            const auto& gamma = _settings.upwind_factor;

            T central = 0.0;
            T upwind = 0.0;

            auto& UVW = _solution.velocity;

            auto base_idx = std::array<size_t,3>({(size_t)i,
                                                  (size_t)j,
                                                  (size_t)k});
            auto fwd_idx1 = base_idx;
            auto fwd_idx2 = base_idx;

            fwd_idx1[dir2] += 1;
            fwd_idx2[dir1] += 1;

            T central_fwd = avg(UVW(base_idx)[dir2], UVW(fwd_idx1)[dir2]);
            central_fwd *= avg(UVW(base_idx)[dir1], UVW(fwd_idx2)[dir1]);

            auto bwd_idx1 = base_idx;
            auto bwdfwd_idx1 = base_idx;
            auto bwd_idx2 = base_idx;

            bwd_idx1[dir1] -= 1;
            bwdfwd_idx1[dir1] -= 1;
            bwdfwd_idx1[dir2] += 1;
            bwd_idx2[dir1] -= 1;

            T central_bwd = avg(UVW(bwd_idx2)[dir2], UVW(bwdfwd_idx1)[dir2]);
            central_bwd *= avg(UVW(bwd_idx2)[dir1], UVW(base_idx)[dir1]);

            central = (1/dh[dir2])*(central_fwd - central_bwd);

            T upwind_fwd = fabs(avg(UVW(base_idx)[dir2], UVW(fwd_idx1)[dir2]));
            upwind_fwd *= avg(UVW(base_idx)[dir1], -UVW(fwd_idx2)[dir1]);

            T upwind_bwd = fabs(avg(UVW(bwd_idx2)[dir2], UVW(bwdfwd_idx1)[dir2]));
            upwind_bwd *= avg(UVW(bwd_idx2)[dir1], -UVW(base_idx)[dir1]);

            upwind = (1/dh[dir2])*(upwind_fwd - upwind_bwd);

            return central + gamma*upwind;
        }

        /**********************************************************************
         * Apply the discrete advection operator
        **********************************************************************/
        inline vector3d_t advection(int i, int j, int k){

            const auto& dim = _problem.geometry().dimension();

            // For each velocity component, we evaluate the spatial
            // difference in each Cartesian direction.
            vector3d_t advec = {0.0, 0.0, 0.0};

            for(int d1=0; d1<dim; d1++){
                for(int d2=0; d2<dim; d2++){
                    advec(d1) -= upwind_difference(d1, d2, i, j, k);
                }
            }
            return advec;
        }

        /**********************************************************************
         * Apply the discrete diffusion (∇²) operator to the solution
         * velocity field.
        **********************************************************************/
        inline vector3d_t diffusion(int i, int j, int k){

            // We need the velocity field and cell sizes to compute
            // the Laplacian.
            auto& UVW = _solution.velocity;

            const auto& h = _problem.geometry().cell_sizes();
            const auto& dim = _problem.geometry().dimension();

            // To simplify the looping, we make index arrays to hold the
            // indices we need for central differencing.
            auto base_idx = std::array<size_t,3>({(size_t)i,
                                                  (size_t)j,
                                                  (size_t)k});

            auto fwd_idx = base_idx;
            auto bwd_idx = base_idx;

            // For each velocity component, accumulate the Laplacian terms
            // in each Cartesian direction.
            vector3d_t diff(0.0, 0.0, 0.0);
            vector3d_t h_sq(h[0]*h[0], h[1]*h[1], h[2]*h[2]);

            for(int d=0; d<dim; d++){

                fwd_idx[d] += 1;
                bwd_idx[d] -= 1;

                // Compute ∂²/∂h² for the component, where h = x,y,or z
                diff += (UVW(fwd_idx) - 2.0*UVW(base_idx) + UVW(bwd_idx))/h_sq;

                // Reset the indices for the next iteration
                fwd_idx = base_idx;
                bwd_idx = base_idx;
            }
            return diff;
       }

        /**********************************************************************
         * Compute the intermediate velocity terms needed for subsequent
         * pressure iterations. This is done by taking the current velocity
         * and applying the advection-diffusion operator.
        **********************************************************************/
        void compute_intermediate_velocity(){

            auto& UVW = _solution.velocity;
            auto& FGH = _solution.intermediate_velocity;

            // Precompute and cache constants needed in later steps
            const auto& Re = _problem.flow_parameters().Re();
            const auto& inv_Re = 1/Re;
            const auto& dt = _problem.timestepper().dt();
            const auto& bforce = _problem.flow_parameters().body_forces();

            // Precompute the max indices for each dimension of the problem
            const auto& imax = _problem.geometry().nx()-1;
            const auto& jmax = _problem.geometry().ny()-1;
            const auto& kmax = _problem.geometry().nz()-1;

            // Apply the discrete advection-diffusion operator to each
            // non-boundary element in the problem.
            #pragma omp parallel for
            for(int k=1; k<kmax; k++){
                for(int j=1; j<jmax; j++){
                    for(int i=1; i<imax; i++){

                        auto diff = diffusion(i,j,k);
                        auto advec = advection(i,j,k);
                        auto advec_diff = dt*(inv_Re*(diff) + advec + bforce);

                        FGH(i,j,k) = UVW(i,j,k) + advec_diff;
                    }
                }
            }

            // For the boundary cells, copy in the velocities of the adjacent
            // interior cells...

            // ...West, East
            for(int k=0; k<=kmax; k++){
                for(int j=0; j<=jmax; j++){

                    FGH(0,j,k) = UVW(0,j,k);
                    FGH(imax,j,k) = UVW(imax,j,k);
                }
            }
            // ...South, North
            for(int k=0; k<=kmax; k++){
                for(int i=0; i<=imax; i++){

                    FGH(i,0,k) = UVW(i,0,k);
                    FGH(i,jmax,k) = UVW(i,jmax,k);
                }
            }
            // ...Down, Up
            for(int j=0; j<=jmax; j++){
                for(int i=0; i<=imax; i++){

                    FGH(i,j,0) = UVW(i,j,0);
                    FGH(i,j,kmax) = UVW(i,j,kmax);
                }
            }
        }

        /**********************************************************************
         * Pre-compute the right-hand side of the pressure equation.
        **********************************************************************/
        void compute_rhs(){

            const auto& dt = _problem.timestepper().dt();
            const auto& dh = _problem.geometry().cell_sizes();

            const auto& imax = _problem.geometry().nx()-1;
            const auto& jmax = _problem.geometry().ny()-1;
            const auto& kmax = _problem.geometry().nz()-1;

            auto& FGH = _solution.intermediate_velocity;
            auto& RHS = _solution.rhs;

            auto inv_dt = 1/dt;
            vector3d_t inv_h = {1/dh[0], 1/dh[1], 1/dh[2]};
            vector3d_t vsum = {0.0, 0.0, 0.0};

            #pragma omp parallel for
            for(int k=1; k<kmax; k++){
                for(int j=1; j<jmax; j++){
                    for(int i=1; i<imax; i++){

                        vsum(0) = FGH(i,j,k).u() - FGH(i-1,j,k).u();
                        vsum(1) = FGH(i,j,k).v() - FGH(i,j-1,k).v();
                        vsum(2) = FGH(i,j,k).w() - FGH(i,j,k-1).w();

                        vsum *= inv_h;

                        RHS(i,j,k) = inv_dt*vsum.sum();
                    }
                }
            }
        }

        /**********************************************************************
         * Perform a single SOR iteration on the pressure field.
        **********************************************************************/
        void sor_iteration(){

            _stats.reset_norm();

            const auto& dt = _problem.timestepper().dt();
            const auto& dh = _problem.geometry().cell_sizes();
            const auto& dim = _problem.geometry().dimension();
            const auto& omega = _settings.upwind_factor;

            const auto& imax = _problem.geometry().nx()-1;
            const auto& jmax = _problem.geometry().ny()-1;
            const auto& kmax = _problem.geometry().nz()-1;

            auto& RHS = _solution.rhs;
            auto& P = _solution.pressure;

            auto psum = std::inner_product(P.data().begin(),
                                           P.data().end(),
                                           P.data().begin(),
                                           0.0);

            // Copy the adjacent pressures into the boundary cells...

            // ...West, East
            for(int k=0; k<=kmax; k++){
                for(int j=0; j<=jmax; j++){

                    P(0,j,k) = P(1,j,k);
                    P(imax,j,k) = P(imax-1,j,k);
                }
            }
            // ...South, North
            for(int k=0; k<=kmax; k++){
                for(int i=0; i<=imax; i++){

                    P(i,0,k) = P(i,1,k);
                    P(i,jmax,k) = P(i,jmax-1,k);
                }
            }
            // ...Down, Up
            for(int j=0; j<=jmax; j++){
                for(int i=0; i<=imax; i++){

                    P(i,j,0) = P(i,j,1);
                    P(i,j,kmax) = P(i,j,kmax-1);
                }
            }

            // Pre-compute the squares of the mesh sizes
            std::array<T,3> dh2;
            for(int d=0; d<dim; d++){
                dh2[d] = dh[d]*dh[d];
            }

            // Pre-compute overrelaxation coefficients
            // b1 = (1-omega)
            // b2 = omega/(2/dx^2 + 2/dy^2 +...)
            T D1 = 0.0;
            for(const auto& val: dh2){
                D1 += 2.0/val;
            }
            T b1 = 1-omega;
            T b2 = omega/D1;

            // Solve for the pressures in the internal cells
            T del_p(0.0);
            #pragma omp parallel for   
            for(int k=1; k<kmax; k++){
                for(int j=1; j<jmax; j++){
                    for(int i=1; i<imax; i++){

                        del_p = 0.0;

                        del_p += (P(i+1,j,k) + P(i-1,j,k))/dh2[0];
                        del_p += (P(i,j+1,k) + P(i,j-1,k))/dh2[1];
                        del_p += (P(i,j,k+1) + P(i,j,k-1))/dh2[2];

                        P(i,j,k) = b1*P(i,j,k) + b2*(del_p - RHS(i,j,k));
                    }
                }
            }

            compute_norms();

            _stats.l2_norm /= psum;
            _stats.linf_norm /= psum;
        }

        /**********************************************************************
         * Compute the L2 and L-inf norms.
        **********************************************************************/
        void compute_norms(){

            auto nelem = _problem.geometry().total_cells();

            const auto& res_sq = _stats.residual_sumsquares;
            const auto& res_max = _stats.residual_max;

            auto& l2 = _stats.l2_norm;
            l2 = sqrt(res_sq/nelem);

            auto& l_inf = _stats.linf_norm;
            l_inf = _stats.residual_max;
        }

        /**********************************************************************
         * Solve for the pressure field using the successive over-relaxation
         * (SOR) method. Iterate until the L2 norm reaches a tolerance
         * threshold, or until the solver reaches a maximum number of
         * iterations.
        **********************************************************************/
        void solve_pressure(){

            _stats.reset();
            auto& iter = _stats.itercount;
            const auto& max_iters = _settings.max_iters;

            do{

                sor_iteration();

                iter++;

            }while ((iter < max_iters) &&
                    (_stats.l2_norm > _settings.tolerance));
            std::cout << "\t\tL_2 = "
                      << _stats.l2_norm
                      << "\n\t\tL_inf = "
                      << _stats.linf_norm
                      << "\n\t\titers = "
                      << iter
                      << "\n";
        }

        /**********************************************************************
         * Update the velocity field based on the computed pressure field.
        **********************************************************************/
        void update_velocity(){

            const auto& dt = _problem.timestepper().dt();
            const auto& dh = _problem.geometry().cell_sizes();
            const auto& dim = _problem.geometry().dimension();

            auto& P = _solution.pressure;

            const auto& imax = _problem.geometry().nx();
            const auto& jmax = _problem.geometry().ny();
            const auto& kmax = _problem.geometry().nz();

            auto& UVW = _solution.velocity;
            auto& FGH = _solution.intermediate_velocity;

            T dtdx(dt/dh[0]), dtdy(dt/dh[1]), dtdz(dt/dh[2]);

            #pragma omp parallel for   
            for(int k=1; k<kmax-1; k++){
                for(int j=1; j<jmax-1; j++){
                    for(int i=1;i<imax-1; i++){

                        UVW(i,j,k).u() = \
                            FGH(i,j,k).u() - dtdx*(P(i+1,j,k) - P(i,j,k));
                        UVW(i,j,k).v() = \
                            FGH(i,j,k).v() - dtdy*(P(i,j+1,k) - P(i,j,k));
                        UVW(i,j,k).w() = \
                            FGH(i,j,k).w() - dtdz*(P(i,j,k+1) - P(i,j,k));
                    }
                }
            }
        }

        /**********************************************************************
         * Perform a single timestep. This consists of applying boundary
         * conditions and motions, computing the intermediate velocity field
         * based on an advection-diffusion operator, converging the pressure
         * field using the SOR method, then finally updating the velocity
         * field.
        **********************************************************************/
        void step(){

            // Apply boundary conditions
            _problem.boundaries().apply_conditions(_solution);
            _problem.boundaries().apply_motions(_solution);

            // Compute the intermediate velocities
            compute_intermediate_velocity();

            // Construct the right-hand side of the pressure equation
            compute_rhs();

            // Solve for the pressure
            solve_pressure();

            // Update the velocity components
            update_velocity();

        }

        /**********************************************************************
         * Solve the time-dependent Navier-Stokes equations. We use an explicit
         * Euler method for time stepping, and adaptively advance the time step
         * based on sability criteria. At each step, we converge the pressure
         * field using the successive over-relaxation (SOR) method then project
         * this solution to determine the velocity field.
        ***********************************************************************/
        void solve(){

            auto tmax = _problem.timestepper().max_time();

            int step_count = 0;
            if(step_count == 0){

                auto bname = _output_settings.base_name;

                std::string oname = \
                    bname + std::to_string(step_count) + ".vti";

                auto output = VTKFile<T>(_output_settings, oname);

                output.set_geometry(_problem);
                output.store_solution(_solution);
                output.save_file();
            }

            do{
                step_count++;

                std::cout << "\n---------------------------------------------\n";

                std::cout << "\tTimestep " << step_count << " (t = "
                          << _problem.timestepper().current_time()
                          << ", "
                          << "\u0394t = "
                          << _problem.timestepper().dt()
                          << ")\n\n";

                // Solve the current timestep
                step();

                // Adaptively advance the timestepper
                auto cell_sizes = _problem.geometry().cell_sizes();
                auto Re = _problem.flow_parameters().Re();
                auto max_vels = _solution.max_velocity_components();

                _problem.timestepper().advance(cell_sizes, Re, max_vels);

                // Save Timestep Solution
                int write_every = _output_settings.write_every;
                if(step_count % write_every == 0){

                    auto bname = _output_settings.base_name;

                    std::string oname = \
                        bname + std::to_string(step_count) + ".vti";

                    auto output = VTKFile<T>(_output_settings, oname);

                    output.set_geometry(_problem);
                    output.store_solution(_solution);
                    output.save_file();
                }

            } while (_problem.timestepper().current_time() <= tmax);

            std::cout << "---------------------------------------------\n";
            std::cout << "\n*** Solution complete! (t = "
                      << _problem.timestepper().current_time() << ")***\n";
        }
};

#endif