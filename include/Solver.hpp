#ifndef SOLVER_HPP
#define SOLVER_HPP

#include <limits>

#include "Solution.hpp"
#include "Problem.hpp"

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

    /**
     * Default Constructor. Initialize all settings to zero.
    */
    SolverStats(): itercount(0),
                   residual_sumsquares(0),
                   residual_max(0),
                   l2_norm(0),
                   linf_norm(0){
        reset();
    }

    /**
     * Reset all solver statistics.
    */
    void reset(){
        itercount = 0;
        reset_norm();
    }

    /**
     * Reset all solver norms.
    */
    void reset_norm(){
        residual_sumsquares = 0;
        residual_max = 0;
        l2_norm = 0;
        linf_norm = 0;
    }
};

/**
 * Solver class. A collection of functions that operate on a problem
 * definition to solve the Navier-Stokes equations for vector and 
 * scalar fields. Based on Ch. 3 of "Numerical Simulation in Fluid Dynamics"
 * by Griebel, et al.
*/
template<typename T>
class Solver{

    private:

        Problem<T> _problem;
        SolverStats<T> _stats;
        Solution<T> _solution;
        SolverSettings<T> _settings;

    public:

        /**
         * Default constructor (does nothing).
        */
        Solver(){}

        /**
         * Preferred constructor. Uses problem definition and solver settings
         * based on user input.
        */
        Solver(Problem<T>& problem, SolverSettings<T>& settings){

            _problem = problem;
            _settings = settings;
            _solution = Solution<T>(_problem);
        }

        /**
         * Compute the average value of two numeric items.
        */
        T avg(T item1, T item2){
            return 0.5*(item1+item2);
        }

        /**
         * 
        */
        std::string idx_to_str(std::vector<size_t>& idxs){
            std::string str = "(";

            for(const auto& idx: idxs){
                str += std::to_string(idx);
                str += ", ";
            }
            str += ")";
            return str;
        }

        /**
         * 
        */
        T upwind_difference(int axis1, int axis2=0, int i=0, int j=0, int k=0){

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

            T central_fwd = avg(S2(base_idx), S2(fwd_idx1));
            central_fwd *= avg(S1(base_idx), S1(fwd_idx2));

            auto bwd_idx1 = base_idx;
            auto bwdfwd_idx1 = base_idx;
            auto bwd_idx2 = base_idx;

            bwd_idx1[axis1] -= 1;
            bwdfwd_idx1[axis1] -= 1;
            bwdfwd_idx1[axis2] += 1;
            bwd_idx2[axis1] -= 1;

            T central_bwd = avg(S2(bwd_idx2), S2(bwdfwd_idx1));
            central_bwd *= avg(S1(bwd_idx2), S1(base_idx));
        
            central = (1/dh[axis2])*(central_fwd - central_bwd);
            
            T upwind_fwd = fabs(avg(S2(base_idx), S2(fwd_idx1)));
            upwind_fwd *= avg(S1(base_idx), -S1(fwd_idx2));

            T upwind_bwd = fabs(avg(S2(bwd_idx2), S2(bwdfwd_idx1)));
            upwind_bwd *= avg(S1(bwd_idx2), -S1(base_idx));

            upwind = (1/dh[axis2])*(upwind_fwd - upwind_bwd);

            return central + gamma*upwind;
        }

        /**
         * 
        */
        T advection(int axis, int i, int j=0, int k=0){
            
            constexpr bool DEBUG = false;
            const std::string partial = "\u2202";

            if constexpr(DEBUG){
                std::cout << "*** Computing advection term for ";
                if(axis == 0) std::cout << "U (";
                if(axis == 1) std::cout << "V (";
                if(axis == 2) std::cout << "W (";

                std::cout << i << ", " << j << ", " << k << ") ***\n";

            }

            int dim = _problem.geometry().dimension();

            T adv = 0.0;
            T du2dx(0), duvdy(0), duwdz(0);
            T duvdx(0), dv2dy(0), dvwdz(0);
            T duwdx(0), dvwdy(0), dw2dz(0);

            switch(axis){
                case 0:
                    du2dx = upwind_difference(axis, axis, i, j, k);
                    if constexpr (DEBUG){
                        std::cout << partial + "(u\u00b2)/" + partial + "x = ";
                        std::cout << du2dx << "\n";
                    }
                    if(dim > 1) duvdy = upwind_difference(axis, 1, i, j, k);
                    if constexpr (DEBUG){
                        std::cout << partial + "(uv)/" + partial + "y = ";
                        std::cout << duvdy << "\n";
                    }
                    if(dim > 2) duwdz = upwind_difference(axis, 2, i, j, k);
                    if constexpr (DEBUG){
                        std::cout << partial + "(uw)/" + partial + "z = ";
                        std::cout << duwdz << "\n";
                    }

                    adv = -du2dx - duvdy - duwdz;
                    break;

                case 1:

                    duvdx = upwind_difference(axis, 0, i, j, k);
                    if constexpr (DEBUG){
                        std::cout << partial + "(uv)" + partial + "x = ";
                        std::cout << duvdx << "\n";
                    }

                    dv2dy = upwind_difference(axis, axis, i, j, k);
                    if constexpr (DEBUG){
                        std::cout << partial + "(v\u00b2)" + partial + "y = ";
                        std::cout << dv2dy << "\n";
                    }

                    if(dim > 2) dvwdz = upwind_difference(axis, 2, i, j, k);
                    if constexpr (DEBUG){
                        std::cout << partial + "(vw)" + partial + "z = ";
                        std::cout << dvwdz << "\n";
                    }

                    adv = -duvdx - dv2dy - dvwdz;
                    break;
                case 2:
                    duwdx = upwind_difference(axis, 0, i, j, k); 
                    dvwdy = upwind_difference(axis, 1, i, j, k);
                    dw2dz = upwind_difference(axis, axis, i,j,k);

                    adv = -duwdx - dvwdy - dw2dz;
                    break;
                default:
                    throw std::invalid_argument("Invalid dimension.");
            }

            return adv;
        }

        /**
         * Old diffusion operator. Save until new one is verified.
        */
        /*T diffusion(int axis, int i, int j=0, int k=0){

            constexpr bool DEBUG = false;
            const std::string delsquared = "\u2207\u00b2";
            if constexpr (DEBUG){
                

                std::cout << "*** Computing " << delsquared << " for ";
                if(axis == 0) std::cout << "U ";
                if(axis == 1) std::cout << "V ";
                if(axis == 2) std::cout << "W ";

                std::cout << "(" << i << ", " << j << ", " << k << ") ***\n";
            }


            auto dim = _problem.geometry().dimension();
            auto dh = _problem.geometry().cell_sizes();

            auto& UVW = _solution.velocity_component(axis);

            T ans = 0.0;            
            T num = 0.0;

            for(int d=0; d<dim; d++){
                switch(d){
                    case 0:
                        // d^2x
                        num = (UVW(i+1,j,k) + 2*UVW(i,j,k) + UVW(i-1,j,k));

                        if constexpr (DEBUG){
                            std::cout << delsquared << "x = " << num << "\n";
                        }
                        break;

                    case 1:
                        // d^2y
                        num = (UVW(i,j+1,k) + 2*UVW(i,j,k) + UVW(i,j-1,k));

                        if constexpr (DEBUG){
                            std::cout << delsquared << "y = " << num << "\n";
                        }
                
                        break;

                    case 2:
                        // d^2z
                        num = (UVW(i,j,k+1) + 2*UVW(i,j,k) + UVW(i,j,k-1));
                        if constexpr (DEBUG){
                            std::cout  << delsquared << "z = " << num << "\n";
                        }
                        break;
                    default:
                        break;
                }

                if constexpr (DEBUG){
                    std::cout << "h = " << dh[d] << " h\u00b2 = " << dh[d]*dh[d] << "\n";
                }
                // divide by dx^2, dy^2 or dz^2, as appropriate.
                ans += num/(dh[d]*dh[d]);
            }

            if constexpr (DEBUG){
                std::cout << delsquared << "x ";
                if (dim > 1) std::cout << "+ " << delsquared << "y ";
                if (dim > 2) std::cout << "+ " << delsquared << "z ";

                std::cout << "= " << ans << "\n\n";
            }

            return ans;
       }*/

        /**
         * Apply the diffusion (Laplacian) operator to the solution velocity
         * field.
        */
        inline vector3d_t& diffusion(int i, int j, int k){

            // We need the velocity field and cell sizes to compute
            // the Laplacian.
            const auto& UVW = _solution.velocity;

            const auto& h = _problem.geometry().cell_sizes();
            const auto& dim = _problem.geometry().dimension();

            // To simplify the looping, we make index arrays to hold the
            // indices we need for central differencing.
            std::array<int, 3> base_idx = {i,j,k};
            auto fwd_idx = base_idx;
            auto bwd_idx = base_idx;

            // For each velocity component, accumulate the Laplacian terms
            vector3d_t delsquared(0.0, 0.0, 0.0);
            
            for(int d=0; d<dim; d++){

                // Set the forward and backward indices needed for the
                // central differencing
                fwd_idx[d] += 1;
                bwd_idx[d] -= 1;

                // Compute d^2/dh^2 for the component, where h = x,y,or z
                delsquared[d] = (UVW(fwd_idx)[d] + \
                                 2*UVW(base_idx)[d] + \
                                 UVW(bwd_idx)[d]) / (h[d]*h[d]);

                // Reset the indices for the next iteration
                fwd_idx = base_idx;
                bwd_idx = base_idx;
            }
            return delsquared;
       }

        /**
         * 
        */
        void compute_momenta(){

            auto dim = _problem.geometry().dimension();

            // Start by copying the velocity components into the FGH arrays
            for(int d=0; d<dim; d++){

                auto& UVW = _solution.velocity_component(d);
                auto& FGH = _solution.tentative_momentum(d);

                for(int i=0; i<UVW.data().size(); i++){
                    FGH.data()[i] = UVW.data()[i];
                }
            }

            // Compute the starting momenta in the internal region of the problem
            auto Re = _problem.flow_parameters().Re();
            auto inv_Re = 1/Re;
            auto dt = _problem.timestepper().dt();

            auto imax = _problem.geometry().ncells()[0]-1;
            auto jmax = (dim == 2) ? _problem.geometry().ncells()[1]-1 : 1;
            auto kmax = (dim == 3) ? _problem.geometry().ncells()[2]-1 : 1;

            for(int d=0; d<dim; d++){
                auto& UVW = _solution.velocity_component(d);
                auto& FGH = _solution.tentative_momentum(d);
                auto& body_forces = _problem.flow_parameters().body_forces();

                int k = (dim > 2) ? 1 : 0;
                do{
                    int j = (dim > 1) ? 1 : 0;
                    do{
                        for(int i=1; i<imax; i++){

                            auto diff = diffusion(d,i,j,k);
                            auto advect = advection(d,i,j,k);
                            
                            FGH(i,j,k) += dt*((inv_Re)*diff + advect + body_forces[d]);
                        }
                        j++;
                    }while(j<jmax);
                    k++;
                }while(k<kmax);
            }
        }

        /**
         * Compute the intermediate velocity terms needed for subsequent
         * pressure iterations. This is done by taking the current velocity
         * and applying the advection-diffusion operator.
        */
        void compute_intermediate_velocity(){

            const auto& UVW = _solution.velocity;
            auto& FGH = _solution.intermediate_velocity;

            // Precompute and cache constants needed in later steps 
            const auto& Re = _problem.flow_parameters().Re();
            const auto& inv_Re = 1/Re;
            const auto& dt = _problem.timestepper().dt();
            const auto& bforce = _problem.flow_parameters().body_forces();

            // Precompute the max indices for each dimension of the problem
            const auto& imax = _problem.geometry().ncells()[0]-1;
            const auto& jmax = _problem.geometry().ncells()[1]-1;
            const auto& kmax = _problem.geometry().ncells()[2]-1;

            // Apply the advection-diffusion operator to each non-boundary
            // element in the problem.

            // Each layer in the block
            for(int k=1; k<kmax; k++){
                // Each row in a layer
                for(int j=1; j<jmax; j++){
                    // Each column in a row
                    for(int i=1; i<imax; i++){

                        auto diff = diffusion(i,j,k);
                        auto advect = advection(i,j,k);
                        auto op = dt*(inv_Re*(diff) - advect + bforce);

                        FGH(i,j,k) = UVW(i,j,k) + op;
                    }
                }
            }

            // For the boundary cells, copy in the velocities of the adjacent
            // interior cells. 
            


        }

        /**
         * Pre-compute the right-hand side of the pressure equation
        */
        void compute_rhs(){

            auto dt = _problem.timestepper().dt();
            auto dh = _problem.geometry().cell_sizes();
            auto dim = _problem.geometry().dimension();

            auto& RHS = _solution.rhs;
            RHS.zeros();

            auto imax = _problem.geometry().ncells()[0]-1;
            auto jmax = (dim == 2) ? _problem.geometry().ncells()[1]-1 : 1;
            auto kmax = (dim == 3) ? _problem.geometry().ncells()[2]-1 : 1;

            for(int d=0; d<dim; d++){

                auto& FGH = _solution.tentative_momentum(d);
                auto& UVW = _solution.velocity_component(d);

                int k = (dim > 2) ? 1 : 0;
                do{
                    int j = (dim > 1) ? 1 : 0;
                    do{
                        for(int i=1; i<imax; i++){

                            auto rhs_idx = std::vector<size_t>({(size_t)i,(size_t)j,(size_t)k});
                            rhs_idx[d] -= 1;

                            RHS(i,j,k) += (FGH(i,j,k) - FGH(rhs_idx));
                            RHS(i,j,k) /= dh[d];
                            RHS(i,j,k) /= dt;
                        }
                        j++;
                    }while(j<jmax);
                    k++;
                }while(k<kmax);
            }
        }

        /**
         * Perform a single SOR iteration on the pressure field.
        */
        void sor_iteration(){

            _stats.reset_norm();

            auto dt = _problem.timestepper().dt();
            auto dh = _problem.geometry().cell_sizes();
            auto dim = _problem.geometry().dimension();
            auto omega = _settings.upwind_factor;

            auto& RHS = _solution.rhs;
            auto& P = _solution.pressure;

            auto imax = _problem.geometry().ncells()[0]-1;
            auto jmax = (dim == 2) ? _problem.geometry().ncells()[1]-1 : 1;
            auto kmax = (dim == 3) ? _problem.geometry().ncells()[2]-1 : 1;

            auto psum = std::inner_product(P.data().begin(),
                                           P.data().end(),
                                           P.data().begin(),
                                           0.0);

            // Copy the adjacent pressures into the boundary cells
            int k = (dim > 2) ? 1 : 0;
            do{
                int j = (dim > 1) ? 1 : 0;
                do{
                    for(int i=0; i<=imax; i++){
                        P(0,j,k) = P(1,j,k);
                        P(imax,j,k) = P(imax-1,j,k);

                        P(i,0,k) = P(i,1,k);
                        P(i,jmax,k) = P(i,jmax-1,k);

                        P(i,j,0) = P(i,j,1);
                        P(i,j,kmax) = P(i,j,kmax-1);
                    }
                    j++;
                }while(j<=jmax);
                k++;
            }while (k<=kmax);

            // Pre-compute beta = omega/(2/dx^2 + 2/dy^2 +...)
            std::vector<T> dh2;
            std::transform(dh.begin(),
                           dh.end(),
                           std::back_inserter(dh2),
                            [](T h){return h*h;});
            
            T D1 = 0.0;
            for(const auto& val: dh2){
                D1 += 2.0/val;
            }
            T beta = omega/D1;

            // Solve for the pressures in the internal cells
            k = (dim > 2) ? 1 : 0;
            do{
                int j = (dim > 1) ? 1 : 0;
                do{
                    for(int i=0; i<imax; i++){
                        
                        auto base_idx = std::vector<size_t>({(size_t)i,
                                                             (size_t)j,
                                                             (size_t)k});
                        // Apply the overrelaxation term
                        P(i,j,k) *= (1-omega);
                        
                        // Compute pressure spatial differences and residual 
                        T dpdh = 0.0;
                        T res = 0.0;

                        for (int d=0; d<dim; d++){

                            auto fwd_idx = base_idx;
                            auto bwd_idx = base_idx;
                            fwd_idx[d] += 1;
                            bwd_idx[d] -= 1;

                            dpdh += (P(fwd_idx) + P(bwd_idx))/dh2[d];
                            res += (P(fwd_idx) - P(bwd_idx))/dh2[d];

                        }

                        P(i,j,k) += beta*(dpdh-RHS(i,j,k));

                        // Update residuals
                        res -= RHS(i,j,k);
                        res = fabs(res);

                        auto& res_sq = _stats.residual_sumsquares;
                        res_sq += res*res;

                        auto& res_max = _stats.residual_max;
                        res_max = res > res_max ? res : res_max;
                    }
                    j++;
                }while(j<jmax);
                k++;
            }while (k<kmax);


            compute_norms();

            _stats.l2_norm /= psum;
            _stats.linf_norm /= psum;
        }

        /**
         * Compute the L2 and L-inf norms.
        */
        void compute_norms(){

            auto nelem = _problem.geometry().total_cells();

            const auto& res_sq = _stats.residual_sumsquares;
            const auto& res_max = _stats.residual_max;

            auto& l2 = _stats.l2_norm;
            l2 = sqrt(res_sq/nelem);

            auto& l_inf = _stats.linf_norm;
            l_inf = _stats.residual_max;
        }

        /**
         * Solve for the pressure field using the successive over-relaxation
         * (SOR) method. Iterate until the L2 norm reaches a tolerance
         * threshold, or until the solver reaches a maximum number of
         * iterations.
        */
        void solve_pressure(){

            _stats.reset();

            auto& iter = _stats.itercount;
            const auto& max_iters = _settings.max_iters;

            do{
                std::cout << "iter = " << iter << " of " << max_iters << " ";
                sor_iteration();

                std::cout << "L_2 = " << _stats.l2_norm << "\tL_inf = " << _stats.linf_norm << "\n";
                iter++;
                
            }while ((iter < max_iters) && 
                    (_stats.l2_norm > _settings.tolerance));

            //std::cout << "*** P (iter = " << iter << ")***\n";
            //_solution.pressure.print_field2d(0);

        }

        /**
         * Update the velocity field based on the computed pressure field.
        */
        void update_velocity(){
            
            auto dt = _problem.timestepper().dt();
            auto dh = _problem.geometry().cell_sizes();
            auto dim = _problem.geometry().dimension();

            auto& P = _solution.pressure;

            auto imax = _problem.geometry().ncells()[0]-1;
            auto jmax = _problem.geometry().ncells()[1]-1;
            auto kmax = _problem.geometry().ncells()[2]-1;

            for(int d=0; d<dim; d++){

                if (d == 0) imax = imax - 1;
                else imax = _problem.geometry().ncells()[0]-1;
                if (d == 1) jmax = jmax - 1;
                else jmax = (dim == 2) ? _problem.geometry().ncells()[1]-1 : 1;
                if (d == 2) kmax = kmax - 1;
                else kmax = (dim == 3) ? _problem.geometry().ncells()[2]-1 : 1;

                auto& UVW = _solution.velocity_component(d);
                auto& FGH = _solution.tentative_momentum(d);

                int k = (dim > 2) ? 1 : 0;
                do{
                    int j = (dim > 1) ? 1 : 0;
                    do{
                        for(int i=1; i<=imax; i++){

                            auto base_idx = std::vector<size_t>({(size_t)i,
                                                                (size_t)j,
                                                                (size_t)k});
                            auto fwd_idx = base_idx;
                            fwd_idx[d] += 1;

                            auto dP = P(fwd_idx)-P(base_idx);

                            UVW(i,j,k) = FGH(i,j,k) - (dt/dh[d])*(dP);
                        }
                        j++;
                    }while(j<=jmax);
                    k++;
                }while (k<=kmax);
            }
        }

        /**
         * Perform a single timestep. This consists of applying boundary
         * conditions and motions, computing the intermediate velocity field
         * based on an advection-diffusion operator, converging the pressure
         * field using the SOR method, then finally updating the velocity
         * field. 
        */
        void step(){

            // Apply boundary conditions
            _problem.boundaries().apply_conditions(_solution);
            _problem.boundaries().apply_motions(_solution);

            /*std::cout << "*** U ***\n";
            _solution.U.print_field2d(0);
            std::cout << "*** V ***\n";
            _solution.V.print_field2d(0);*/


            // Compute the intermediate velocities
            compute_intermediate_velocity();

           /* std::cout << "*** F ***\n";
            _solution.F.print_field2d(0);
            std::cout << "*** G ***\n";
            _solution.G.print_field2d(0);*/

            // Construct the right-hand side of the pressure equation
            compute_rhs();

            /*std::cout << "*** RHS ***\n";
            _solution.rhs.print_field2d(0);*/

            // Solve for the pressure
            solve_pressure();

            // Update the velocity components
            update_velocity();

            std::cout << "*** U ***\n";
            _solution.U.print_field2d(0);
            std::cout << "*** V ***\n";
            _solution.V.print_field2d(0);
        }

        /**
         * Solve the time-dependent Navier-Stokes equations. We use an explicit
         * Euler method for time stepping, and adaptively advance the time step
         * based on sability criteria. At each step, we converge the pressure
         * field using the successive over-relaxation (SOR) method then project
         * this solution to determine the velocity field.
        */
        void solve(){

            auto tmax = _problem.timestepper().max_time();

            int step_count = 0;
            do{
                step_count++;

                std::cout << "Timestep " << step_count << " (t = "
                          << _problem.timestepper().current_time()
                          << ")\n";

                // Solve the current timestep
                step();
                //if (step_count == 1)
                //            throw std::runtime_error("Stopping iteration.");

                // Adaptively advance the timestepper
                auto dim = _problem.geometry().dimension();
                auto cell_sizes = _problem.geometry().cell_sizes();
                auto Re = _problem.flow_parameters().Re();
                auto max_vels = _solution.max_velocity_components(dim);

                std::cout << "advancing timestep\n";
                _problem.timestepper().advance(cell_sizes, Re, max_vels);

                // Save Timestep Solution
                //if(step_count % write_every == 0){
                //  _solution.to_vtk(fname);
                //}

            } while (_problem.timestepper().current_time() <= tmax);

            std::cout << "Solution complete! (t = " 
              << _problem.timestepper().current_time() << ")\n";
        }
};

#endif