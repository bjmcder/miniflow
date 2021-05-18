#ifndef __BOUNDARY_CONDITIONS_HPP
#define __BOUNDARY_CONDITIONS_HPP

#include <vector>

#include "NDArray.hpp"
#include "Solution.hpp"

template<typename T>
struct Solution;

template<typename T>
class BoundaryConditions{

    private:
        std::vector<int> _conditions;
        std::vector<T> _motion;

    public:

        enum BOUNDARY{WEST=0, EAST, SOUTH, NORTH, DOWN, UP};
        enum BOUNDARY_TYPES{FREE_SLIP=0, NO_SLIP, FLOW, PERIODIC};

        /**
         * 
        */
        BoundaryConditions(){}

        /**
         * 
        */
        BoundaryConditions(std::vector<int>& conditions, 
                           std::vector<T>& motion){

            _conditions = conditions;
            _motion = motion;
        }


        /**
         * 
        */
        void apply_noslip(NDArray<T>& field, int boundary){

            switch(boundary){

                case WEST:
                    for(int i=0; i<field.shape[0]; i++){
                        
                    }
                case EAST:
                case SOUTH:
                case NORTH:
                case DOWN:
                case UP:
                default:
                    break;
            }

        }

        /**
         * 
        */
        void apply_freeslip(NDArray<T>& field){}

        /**
         * 
        */
        void apply_periodic(NDArray<T>& field){}

        /**
         * 
        */
        void apply_outflow(NDArray<T>& field){}

        /**
         * 
        */
        void apply_inflow(NDArray<T>& field){}

        /**
         * 
        */
        void apply_conditions(Solution<T>& sol){
            
            auto dim = _conditions.size()/2; 

            // Check each boundary
            for (int i=0; i<_conditions.size(); i++){
                switch(_conditions[i]){
                    
                    case FREE_SLIP:
                        break;

                    case NO_SLIP:
                        apply_noslip(sol.U, i);
                        if(dim > 1) apply_noslip(sol.V, i);
                        if(dim > 2) apply_noslip(sol.W, i);
                        break;

                    case FLOW:
                        break;

                    case PERIODIC:
                        break;
                      

                    default:
                        throw std::invalid_argument("apply_conditions(): " 
                                                    "Invalid boundary type.\n");
                }
            }
        }

        /**
         * 
        */
        void apply_motions(Solution<T>& sol){

        }
};

#endif