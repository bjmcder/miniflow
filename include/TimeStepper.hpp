#ifndef __TIMESTEPPER_HPP
#define __TIMESTEPPER_HPP

#include <limits>

/**
 * Timestepper class. Controls the tracking of time dependence in the problem,
 * and adaptively advances to the next timestep based on the the
 * Courant-Friedrichs-Levy criterion. 
*/
template<typename T>
class TimeStepper{

    private:
        T _current_time;
        T _dt;
        T _dt0;
        T _max_time;
        T _step_factor;

    public:

        /**
         * Default constructor (does nothing).
        */
       TimeStepper(){}

        /**
         * Preferred constructor from user specified time limit, adaptive
         * safety factor, and initial step size.
        */
        TimeStepper(const T& t_max, const T& tau, const T& dt_0 = 0){

            _max_time = t_max;
            _step_factor = tau;
            _dt = dt_0;
            _dt0 = dt_0;
            _current_time = 0;
        }

        /**
         * Return the current step size.
        */
        const T& dt(){return _dt;}

        /**
         * Return the current time.
        */
        const T& current_time(){return _current_time;}

        /**
         * Return the time limit of the simulation.
        */
        const T& max_time(){return _max_time;}
 
        /**
         * Advance the timestep adaptively by choosing the minimum among
         * several stability criteria.
        */
        void advance(const std::array<T, 3>& cell_sizes,
                     const T& reynolds,
                     const std::array<T, 3>& max_vels){

            auto tau = _step_factor;
            auto h = cell_sizes;
            auto Re = reynolds;

            // First, we compute the characteristic diffusion time, which
            // places a lower bound on the timestep due to diffusion across
            // mesh cells. This term is given by:
            //
            // t_d = (Re/2)*1/(1/dx^2 + 1/dy^2 + 1/dz^2)
            //

            // Precompute the Reynolds/2 term.
            auto Re_2 = 0.5*Re;

            // Precompute the squares of the cell sizes in each dimension.
            std::array<T,3> h2;
            for(int i=0; i<h.size(); i++){
                h2[i] = h[1]*h[i];
            }

            // Precompute the sum of the reciprocals of the squares 
            auto recip_sum = 0.0;
            for (auto& sq: h2){
                recip_sum += (1/sq);
            }

            // Take the reciprocal of the sum of the reciprocals
            auto A2 = 1/recip_sum;

            // Now put it all together to compute the characteristic
            // diffusion time.
            auto diff_time = Re_2*A2;

            // Now, we add the criteria that the timestep can be no greater
            // than the shortest time it takes a hypothetical particle to advect
            // across a mesh cell in any of the three Cartesian directions.
            // This is the Courant-Friedrichs-Levy criterion.

            // Collect the ratio of the mesh size to the absolute value of 
            // the maximum velocity component in each Cartesian direction:
            //
            // t_x = dx/|u_max|, t_y = dy/|v_max|, t_z=dz/|w_max|
            //
            // Note: the absolute value of the components is assumed to be 
            // pre-computed by the solver.
            std::vector<T> comparands;
            auto max_components = max_vels;

            for (int i=0; i<max_components.size(); i++){

                // If we're going to wind up dividing by zero, store the maximum
                // finite value of that numeric type. 
                if(max_components[i] == 0.0){
                    comparands.push_back(std::numeric_limits<T>::max());
                }
                else comparands.push_back(cell_sizes[i]/max_components[i]);
            }
            
            // Add the diffusion time we computed in the first step to our
            // collection of comparison criteria.
            comparands.push_back(diff_time);

            // Finally, compute the new timestep size by choosing the minimum
            // of all the criteria we computed, multiplied by a safety factor:
            //
            // dt = tau*min(t_x, t_y, t_z, t_d)
            //
            auto dt_new = tau * (*std::min_element(comparands.begin(),
                                                   comparands.end()));

            // If the next timestep would get set to zero or a negative value, 
            // then don't update it. Keep the user-defined one instead.
            if(dt_new <= 0.0) dt_new = _dt0;

            // Store the updated timestep size and advance the current time.
            _dt = dt_new;
            _current_time += _dt; 
        }
};

#endif