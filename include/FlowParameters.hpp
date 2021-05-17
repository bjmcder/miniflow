#ifndef __FLOW_PARAMETERS_HPP
#define __FLOW_PARAMETERS_HPP

template<typename T>
struct FlowParameters{
    
    T _reynolds;
    std::vector<T> _initial_velocities;
    std::vector<T> _body_forces;

    FlowParameters(){}

    /**
     * 
    */
    FlowParameters(int dim){
        _reynolds = 0;
        _initial_velocities.resize(dim);
        _body_forces.resize(dim);
    }

    /**
     * 
    */
    T& Re() const{
        return _reynolds;
    }

};

#endif