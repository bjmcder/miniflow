#ifndef FLOW_PARAMETERS_HPP
#define FLOW_PARAMETERS_HPP

#include "Types.hpp"

template<typename T>
class FlowParameters{
    
    private:
        T _reynolds;
        vector3d_t _initial_velocities;
        vector3d_t _body_forces;

    public:

    // Default constructor (does nothing)
    FlowParameters(){}

    /**
     * Preferred constructor. Sets the Reynolds number, initial bulk
     * velocities and bulk body forces.
    */
    FlowParameters(T& reynolds,
                   std::vector<T>& initial_velocities,
                   std::vector<T>& body_forces){
                       
        _reynolds = reynolds;
        _initial_velocities = initial_velocities;
        _body_forces = body_forces;
    }

    /**
     * Return the Reynolds number defined for this problem.
    */
    const T& Re(){return _reynolds;}

    /**
     * Return the initial bulk velocities defined for this problem.
    */
    const vector3d_t& initial_velocities(){return _initial_velocities;}

    /**
     * Return the bulk body force vector for this problem.
    */
    const vector3d_t& body_forces(){return _body_forces;}
};

#endif