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
    FlowParameters(int dim,
                   T& reynolds,
                   std::vector<T>& initial_velocities,
                   std::vector<T>& body_forces){
                       
        _reynolds = reynolds;
        _initial_velocities= initial_velocities;
        _body_forces.resize(dim);
    }

    /**
     * 
    */
    const T& Re(){
        return _reynolds;
    }

    /**
     * 
    */
    const std::vector<T>& initial_velocities(){
        return _initial_velocities;
    }

    const std::vector<T>& body_forces(){
        return _body_forces;
    }

};

#endif