#ifndef PROBLEM_HPP
#define PROBLEM_HPP

#include "BoundaryConditions.hpp"
#include "Geometry.hpp"
#include "TimeStepper.hpp"
#include "FlowParameters.hpp"

// Forward declarations
template<typename T>
class BoundaryConditions;

template<typename T>
class Geometry;

template<typename T>
class TimeStepper;

template<typename T>
struct FlowParameters;

/**
 * Problem class. This is a container for all the data and parameters needed
 * to define a Navier-Stokes flow problem.
*/
template<typename T>
class Problem{

    private:
        Geometry<T> _geom;
        TimeStepper<T> _tstepper;
        BoundaryConditions<T> _bcs;
        FlowParameters<T> _flowparams;

    public:
        /**
         * Default constructor (does nothing).
        */
        Problem(){}

        /**
         * Preferred constructor. Defines a 3D Navier-Stokes problem from
         * the geometry, timestepper, boundary conditions and flow parameters.
        */
        Problem(const Geometry<T>& geom,
                const TimeStepper<T>& tstepper,
                const BoundaryConditions<T>& bcs,
                const FlowParameters<T>& flowparams){

            _geom = geom;
            _tstepper = tstepper;
            _bcs = bcs;
            _flowparams = flowparams;
        }

        /**
         * Return the geometry object.
        */
        Geometry<T>& geometry(){return _geom;}

        /**
         * Return the timestepper object.
        */
        TimeStepper<T>& timestepper(){return _tstepper;}

        /**
         * Return the boundary conditions object.
        */
        BoundaryConditions<T>& boundaries(){return _bcs;}

        /**
         * Return the flow parameters object.
        */
        FlowParameters<T>& flow_parameters(){return _flowparams;}
};

#endif