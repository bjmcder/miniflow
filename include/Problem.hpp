#ifndef __PROBLEM_HPP
#define __PROBLEM_HPP

#include "BoundaryConditions.hpp"
#include "Geometry.hpp"
#include "TimeStepper.hpp"
#include "FlowParameters.hpp"

template<typename T>
class Problem{

    private:
    
        Geometry<T> _geom;
        TimeStepper<T> _tstepper;
        BoundaryConditions<T> _bcs;
        FlowParameters<T> _flowparams;

    public:

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
     * 
    */
    Geometry<T>& geometry() const{
        return _geom;
    }

    /**
     * 
    */
    TimeStepper<T>& timestepper() const{
        return _tstepper;
    }

    /**
     * 
    */
    BoundaryConditions<T>& boundaries() const{
        return _bcs;
    }

    /**
     * 
    */
    FlowParameters<T>& flow_parameters() const{
        return _flowparams;
    }

};

#endif