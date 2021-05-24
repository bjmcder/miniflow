#ifndef __BOUNDARY_CONDITIONS_HPP
#define __BOUNDARY_CONDITIONS_HPP

#include <vector>

#include "NDArray.hpp"
#include "Solution.hpp"

template<typename T>
struct Solution;

template<typename T>
class NDArray;

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
        void apply_noslip(Solution<T>& field, int boundary){

            int dim = field.shape.size();
            int kmax = 0;

            switch(boundary){

                case WEST:

                    switch(dim){
                        case 1:
                            field.U(0) = -field.U(1);
                            break;
                        case 2:
                            for(int j=0; j<field.shape[1]; j++){
                                field.U(0,j) = -field.U(1,j);
                                field.V(0,j) = 0.0;
                            }
                            break;
                        case 3:
                            for(int i=0; i<field.shape[2]; i++){
                                for(int j=0; j<field.shape[1]; j++){
                                    field.U(0,j,i) = -field.U(1,j,i);
                                    field.V(0,j,i) = 0.0;
                                    field.W(0,j,i) = 0.0;
                                }
                            }
                            break;
                        default:
                            throw std::invalid_argument("Invalid dimension");
                    }
                    break;
                    
                case EAST:

                    kmax = field.shape[0]-1;

                    switch(dim){

                        case 1:
                            field.U(kmax) = -field.U(kmax-1);
                            break;

                        case 2:
                            for(int j=0; j<field.shape[1]; j++){
                                field.U(kmax, j) = -field.U(kmax-1,j);
                                field.V(kmax, j) = 0.0; 
                            }
                            break;

                        case 3:
                            for(int i=0; i<field.shape[2]; i++){
                                for(int j=0; j<field.shape[1]; j++){
                                    field.U(kmax,j,i) = -field.U(kmax-1,j,i);
                                    field.V(kmax,j,i) = 0.0;
                                    field.W(kmax,j,i) = 0.0;
                                }
                            }
                            break;

                        default:
                            throw std::invalid_argument("Invalid dimension");
                    }
                    break;

                case SOUTH:

                    switch(dim){
                        case 2:
                            for(int j=0; j<field.shape[0]; j++){
                                field.U(j,0) = 0.0;
                                field.V(j,0) = -field.V(j,1);
                            }
                            break;

                        case 3:
                            for(int i=0; i<field.shape[2]; i++){
                                for(int j=0; j<field.shape[0]; j++){
                                    field.U(j,0,i) = 0.0;
                                    field.V(j,0,i) = -field.V(j,1,i);
                                    field.W(j,0,i) = 0.0;
                                }
                            }
                            break;

                        default:
                            throw std::invalid_argument("Invalid dimension");
                    }
                    break;

                case NORTH:
                    kmax = field.shape[1]-1;

                    switch(dim){
                        case 2:
                            for(int j=0; j<field.shape[0]; j++){
                                field.U(j,kmax) = 0.0;
                                field.V(j,kmax) = -field.V(j,kmax-1);
                            }
                            break;

                        case 3:
                            for(int i=0; i<field.shape[1]; i++){
                                for(int j=0; j<field.shape[0]; j++){
    
                                    field.U(j,kmax,i) = 0.0;
                                    field.V(j,kmax,i) = -field.V(j,kmax-1,i);
                                    field.W(j,kmax,i) = 0.0;
                                }
                            }
                            break;

                        default:
                            throw std::invalid_argument("Invalid dimension");
                    }
                    break;
                case DOWN:

                    switch (dim){
                        case 3:
                            for(int i=0; i<field.shape[2]; i++){
                                for(int j=0; j<field.shape[0]; j++){

                                    field.U(j,i,0) = 0.0;
                                    field.V(j,i,0) = 0.0;
                                    field.W(j,i,0) = -field.W(j,i,0);

                                }   
                            }
                        break;
                    
                    default:
                        throw std::invalid_argument("Invalid dimension");
                    }
                    break;

                case UP:
                    kmax = field.shape[2]-1;
                    switch (dim){
                        case 3:
                            for(int i=0; i<field.shape[2]; i++){
                                for(int j=0; j<field.shape[0]; j++){

                                    field.U(j,i,kmax) = 0.0;
                                    field.V(j,i,kmax) = 0.0;
                                    field.W(j,i,kmax) = -field.W(j,i,kmax-1);

                                }   
                            }
                        break;
                    
                    default:
                        throw std::invalid_argument("Invalid dimension");
                    }
                    break;
                default:
                    throw std::invalid_argument("Invalid boundary direction.");
            }

        }

        /**
         * 
        */
        void apply_freeslip(Solution<T>& field, int boundary){
            throw std::invalid_argument("Free slip boundary not implemented.");
        }

        /**
         * 
        */
        void apply_periodic(Solution<T>& field, int boundary){
            throw std::invalid_argument("Periodic boundary not implemented.");
        }

        /**
         * 
        */
        void apply_flow(Solution<T>& field, int boundary){
            throw std::invalid_argument("Flow boundary not implemented.");
        }

        /**
         * 
        */
        void apply_conditions(Solution<T>& sol){
            
            auto dim = _conditions.size()/2; 

            // Check each boundary
            for (int i=0; i<_conditions.size(); i++){

                switch(_conditions[i]){
                    
                    case FREE_SLIP:
                        apply_freeslip(sol, i);
                        break;

                    case NO_SLIP:
                        apply_noslip(sol, i);
                        break;

                    case FLOW:
                        apply_flow(sol,i);
                        break;

                    case PERIODIC:
                        apply_periodic(sol, i);
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
        void apply_motions(Solution<T>& field){
            
            int dim = _motion.size()/2;

            for (int i=0; i<_motion.size(); i++){

                int kmax = 0;
                
                if (_motion[i] == 0.0) continue;

                T vel = 2.0*_motion[i];

                switch(i){

                    case WEST:

                        switch(dim){
                            case 1:
                                break;
                            case 2:
                                for(int j=0; j<field.shape[1]; j++){
                                    field.V(0,j) += vel;
                                }
                                break;
                            case 3:
                                for(int i=0; i<field.shape[2]; i++){
                                    for(int j=0; j<field.shape[1]; j++){
                                        field.W(0,j,i) += vel;
                                    }
                                }
                                break;
                            default:
                                throw std::invalid_argument("Invalid dimension");
                        }
                        break;
                        
                    case EAST:

                        kmax = field.shape[0]-1;

                        switch(dim){

                            case 1:
                                break;

                            case 2:
                                for(int j=0; j<field.shape[1]; j++){
                                    field.V(kmax, j) += vel; 
                                }
                                break;

                            case 3:
                                for(int i=0; i<field.shape[2]; i++){
                                    for(int j=0; j<field.shape[1]; j++){
                                        field.V(kmax,j,i) += vel;
                                    }
                                }
                                break;

                            default:
                                throw std::invalid_argument("Invalid dimension");
                        }
                        break;

                    case SOUTH:

                        switch(dim){
                            case 2:
                                for(int j=0; j<field.shape[0]; j++){
                                    field.U(j,0) += vel;
                                }
                                break;

                            case 3:
                                for(int i=0; i<field.shape[2]; i++){
                                    for(int j=0; j<field.shape[0]; j++){
                                        field.U(j,0,i) += vel;
                                    }
                                }
                                break;

                            default:
                                throw std::invalid_argument("Invalid dimension");
                        }
                        break;

                    case NORTH:
                        kmax = field.shape[1]-1;

                        switch(dim){
                            case 2:
                                for(int j=0; j<field.shape[0]; j++){
                                    field.U(j,kmax) += vel;
                                }
                                break;

                            case 3:
                                for(int i=0; i<field.shape[1]; i++){
                                    for(int j=0; j<field.shape[0]; j++){
        
                                        field.U(j,kmax,i) += vel;
                                    }
                                }
                                break;

                            default:
                                throw std::invalid_argument("Invalid dimension");
                        }
                        break;
                    case DOWN:

                        switch (dim){
                            case 3:
                                for(int i=0; i<field.shape[2]; i++){
                                    for(int j=0; j<field.shape[0]; j++){

                                        field.U(j,i,0) += vel;
                                    }   
                                }
                            break;
                        
                        default:
                            throw std::invalid_argument("Invalid dimension");
                        }
                        break;

                    case UP:
                        kmax = field.shape[2]-1;
                        switch (dim){
                            case 3:
                                for(int i=0; i<field.shape[2]; i++){
                                    for(int j=0; j<field.shape[0]; j++){

                                        field.U(j,i,kmax) += vel;
                                    }   
                                }
                            break;
                        
                        default:
                            throw std::invalid_argument("Invalid dimension");
                        }
                        break;
                    default:
                        throw std::invalid_argument("Invalid boundary direction.");
                }
            } 

        }
};

#endif