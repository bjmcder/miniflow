#ifndef __BOUNDARY_CONDITIONS_HPP
#define __BOUNDARY_CONDITIONS_HPP

#include <vector>

template<typename T>
class BoundaryConditions{

    private:
        std::vector<int> _conditions;
        std::vector<T> _motion;

    public:

        BoundaryConditions(){}

        BoundaryConditions(std::vector<int>& conditions, 
                           std::vector<T>& motion){

            _conditions = conditions;
            _motion = motion;
        }

        void apply_noslip(){}
        void apply_freeslip(){}
        void apply_periodic(){}
        void apply_outflow(){}
        void apply_inflow(){}

        void apply_conditions(){}

        void apply_motions(){}
};

#endif