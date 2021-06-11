#ifndef SOLUTION_HPP
#define SOLUTION_HPP

#include <array>
#include <cmath>
#include <vector>

#include "Types.hpp"
#include "Problem.hpp"

// Forward declarations
template<typename T>
class Problem;

/******************************************************************************
 * Solution class. This is a container for all the data generated during the
 * solution of a Navier-Stokes problem.
******************************************************************************/
template<typename T>
struct Solution{

    std::array<int,3> shape;

    // Primary quantities of interest
    scalar_field_t pressure;
    vector_field_t velocity;
    vector_field_t vorticity;

    //scalar_field_t temperature

    // Intermediate solution fields
    vector_field_t intermediate_velocity;
    scalar_field_t rhs;

    /**************************************************************************
     * Default constructor (does nothing).
    **************************************************************************/
    Solution(){}

    /**
     * 
    */
    Solution(Problem<T>& problem){

        this->shape = problem.geometry().ncells();

        std::vector<size_t> arraysizes;
        arraysizes.resize(this->shape.size());

        // Need to convert int to size_t
        std::transform(this->shape.begin(), 
                       this->shape.end(),
                       arraysizes.begin(),
                        [](int x) { return (size_t)x;});

        set_array_sizes(arraysizes);

        initialize_solution(problem);
    }

    /**************************************************************************
     * Set the size of the solution fields.
    **************************************************************************/
    void set_array_sizes(const std::vector<size_t>& shape){

        velocity = vector_field_t(shape);
        vorticity = vector_field_t(shape);
        intermediate_velocity = vector_field_t(shape);

        pressure = scalar_field_t(shape);
        rhs = scalar_field_t(shape);
    }

    /**************************************************************************
     * Initialize all solution fields to their starting values.
    **************************************************************************/
    void initialize_solution(Problem<T>& problem){
        
        pressure.zeros();
        velocity.fill(problem.flow_parameters().initial_velocities());
        vorticity.zeros();

        intermediate_velocity.zeros();
        rhs.zeros();
    }

    /**************************************************************************
     * Get a field of single vector direction components from a vector field.
    **************************************************************************/
    scalar_field_t component_to_field(vector_field_t& field, int dim){
        auto result = scalar_field_t(field.shape());

        // Since the fields have the same shapes and strides, we can just
        // copy the components into the data array directly without nested
        // looping.
        for(int i=0; i<field.size(); i++){
            result.data()[i] = field.data()[i][dim];
        }
        return result;
    }

    /**************************************************************************
     * Get the field of velocity components in the x-direction (U).
    **************************************************************************/
    scalar_field_t U(){return component_to_field(velocity, 0);}

    /**************************************************************************
     * Get the field of velocity components in the y-direction (V).
    **************************************************************************/
    scalar_field_t V(){return component_to_field(velocity, 1);}

    /**************************************************************************
     * Get the field of velocity components in the z-direction (W).
    **************************************************************************/
    scalar_field_t W(){return component_to_field(velocity, 2);}

    /**************************************************************************
     * Get the field of intermediate velocity components in the 
     * x-direction (F).
    **************************************************************************/
    scalar_field_t F(){return component_to_field(intermediate_velocity, 0);}

    /**************************************************************************
     * Get the field of intermediate velocity components in the 
     * y-direction (G).
    **************************************************************************/
    scalar_field_t G(){return component_to_field(intermediate_velocity, 1);}

    /**************************************************************************
     * Get the field of intermediate velocity components in the 
     * y-direction (H).
    **************************************************************************/
    scalar_field_t H(){return component_to_field(intermediate_velocity, 2);}

    /**************************************************************************
     * Get the field of velocity components in the specified direction.
    **************************************************************************/
    scalar_field_t velocity_component(int dim){
        switch(dim){
            case 0:
                return U();
            case 1:
                return V();
            case 2:
                return W();
            default:
                throw std::invalid_argument("Invalid dimension.");
        }
    }

    /**************************************************************************
     * Get the field of intermediate velocity components in the specified
     *  direction.
    **************************************************************************/
    scalar_field_t intermediate_velocity_component(int dim){
        switch(dim){
            case 0:
                return F();
            case 1:
                return G();
            case 2:
                return H();
            default:
                throw std::invalid_argument("Invalid dimension.");
        }
    }

    /**************************************************************************
     * Get the absolute value of the maximum velocity component in each
     * direction.
    **************************************************************************/
    std::array<T,3> max_velocity_components(){

        std::array<T,3> max_components = {0,0,0};

        // Get each velocity vector
        for(const auto& vec: velocity.data()){

            // Check each component's absolute value. Store it if it's larger
            // than what we've accumulated so far.
            for(int j=0; j<3; j++){
                T val = fabs(vec[j]);
                if(val > max_components[j]) max_components[j] = val;
            }
        }
        return max_components;
    }
};

#endif