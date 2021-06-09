#ifndef BOUNDARY_CONDITIONS_HPP
#define BOUNDARY_CONDITIONS_HPP

#include <vector>

#include "Types.hpp"
#include "Solution.hpp"

// Forward Declarations
template<typename T>
struct Solution;

/******************************************************************************
 * BoundaryConditions - A class for defining and applying boundary conditions
 * to a flow solution field.
 * 
 * Template Parameters
 * -------------------
 * T : Numeric type (float or double recommended)
******************************************************************************/
template<typename T>
class BoundaryConditions{

    private:
        std::array<int, 6> _conditions;
        std::array<vector3d_t, 6> _motion;

    public:

        enum BOUNDARY{WEST=0, EAST, SOUTH, NORTH, DOWN, UP};
        enum BOUNDARY_TYPES{FREE_SLIP=0, NO_SLIP, INFLOW, OUTFLOW, PERIODIC};

        /**********************************************************************
         * Default constructor (Does nothing).
        **********************************************************************/
        BoundaryConditions(){}

        /**********************************************************************
         * Preferred constructor from user-specified boundary conditions
         * and boundary motion vectors.
         * 
         * Parameters
         * ----------
         * condtions : std::vector<int>&
         *  A vector of integer-specified boundary conditions in canonical
         *  order (W,E,S,N,D,U).
         * 
         * motion : std::vector<std::array<T,3>>&
         *  A vector of 3-element arrays specifying velocities of boundaries
         * (e.g. for lid-driven cavity problems).
        **********************************************************************/
        BoundaryConditions(std::vector<int>& conditions, 
                           std::vector<std::array<T,3>>& motion){

            std::copy(conditions.begin(),
                      conditions.end(),
                     _conditions.begin());
            
            std::copy(motion.begin(), motion.end(), _motion.begin());
        }


        /**********************************************************************
         * No-slip boundary condition. Enforces the condition that the
         * velocity component normal to the boundary equals zero. This is done
         * by setting the velocity component in the boundary cell equal to the
         * negative of the component in the immediately adjacent cell.
         * 
         * Example - West boundary:
         * U(0,j,k) = -U(1,j,k)
         * 
         * Parameters
         * ----------
         * field : Solution<T>&
         *  The solution field that the boundary conditions are being applied
         *  to.
         * 
         * boundary : int
         *  The boundary along which we are applying the condition.
         *  
        **********************************************************************/
        void apply_noslip(Solution<T>& field, int boundary){

            auto& UVW = field.velocity;

            const auto& imax = field.shape[0]-1;
            const auto& jmax = field.shape[1]-1;
            const auto& kmax = field.shape[2]-1;

            switch(boundary){

                case WEST:
                    for(int k=0; k<=kmax; k++){
                        for(int j=0; j<=jmax; j++){

                            UVW(0,j,k).u() = -UVW(1,j,k).u();
                            UVW(0,j,k).v() = 0.0;
                            UVW(0,j,k).w() = 0.0; 
                        }
                    }  
                    break;
                    
                case EAST:
                    for(int k=0; k<=kmax; k++){
                        for(int j=0; j<=jmax; j++){

                            UVW(imax,j,k).u() = -UVW(imax-1,j,k).u();
                            UVW(imax,j,k).v() = 0.0;
                            UVW(imax,j,k).w() = 0.0;
                        }
                    }
                    break;

                case SOUTH:
                    for(int k=0; k<=kmax; k++){
                        for(int i=0; i<=imax; i++){

                            UVW(i,0,k).u() = 0.0;
                            UVW(i,0,k).v() = -UVW(i,1,k).v();
                            UVW(i,0,k).w() = 0.0;
                        }
                    }
                    break;

                case NORTH:
                    for(int k=0; k<=kmax; k++){
                        for(int i=0; i<=imax; i++){

                            UVW(i,jmax,k).u() = 0.0;
                            UVW(i,jmax,k).v() = -UVW(i,jmax-1,k).v();
                            UVW(i,jmax,k).w() = 0.0;
                        }
                    }
                    break;

                case DOWN:
                    for(int j=0; j<=jmax; j++){
                        for(int i=0; i<=imax; i++){

                            UVW(i,j,0).u() = 0.0;
                            UVW(i,j,0).v() = 0.0;
                            UVW(i,j,0).w() = -UVW(i,j,0).w();
                        }   
                    }
                    break;

                case UP:
                    for(int j=0; j<=jmax; j++){
                        for(int i=0; i<=imax; i++){

                            UVW(i,j,kmax).u() = 0.0;
                            UVW(i,j,kmax).v() = 0.0;
                            UVW(i,j,kmax).w() = -UVW(i,j,kmax-1).w();
                        }   
                    }
                    break;

                default:
                    throw std::invalid_argument("Invalid boundary direction.");
            }
        }

        /**********************************************************************
         * Free-slip boundary condition.
        **********************************************************************/
        void apply_freeslip(Solution<T>& field, int boundary){
            throw std::invalid_argument("Free slip boundary not implemented.");
        }

        /**********************************************************************
         * Periodic boundary condition.
        **********************************************************************/
        void apply_periodic(Solution<T>& field, int boundary){
            throw std::invalid_argument("Periodic boundary not implemented.");
        }

        /**********************************************************************
         * Inflow boundary condition.
        **********************************************************************/
        void apply_inflow(Solution<T>& field, int boundary){
            throw std::invalid_argument("Inflow boundary not implemented.");
        }

        /**********************************************************************
         * Outflow boundary condition.
        **********************************************************************/
        void apply_outflow(Solution<T>& field, int boundary){
            throw std::invalid_argument("Outflow boundary not implemented.");
        }

        /**********************************************************************
         * Apply all boundary conditions to the solution field.
         * 
         * Parameters
         * ----------
         * field : Solution<T>&
         *  The object containg the solution field.
        **********************************************************************/
        void apply_conditions(Solution<T>& field){

            // Check each boundary
            for (int i=0; i<_conditions.size(); i++){

                switch(_conditions[i]){
                    
                    case FREE_SLIP:
                        apply_freeslip(field, i);
                        break;

                    case NO_SLIP:
                        apply_noslip(field, i);
                        break;

                    case INFLOW:
                        apply_inflow(field,i);
                        break;

                    case OUTFLOW:
                        apply_outflow(field,i);
                        break;

                    case PERIODIC:
                        apply_periodic(field, i);
                        break;
                      
                    default:
                        throw std::invalid_argument("apply_conditions(): " 
                                                    "Invalid boundary type.\n");
                }
            }
        }

        /**********************************************************************
         * Applies a tangential velocity to each boundary, if specified in the
         * user's input.
         * 
         * Parameters
         * ----------
         * field : Solution<T>&
         *  The object containg the solution field.
        **********************************************************************/
        void apply_motions(Solution<T>& field){

            auto& UVW = field.velocity;

            const auto& imax = field.shape[0]-1;
            const auto& jmax = field.shape[1]-1;
            const auto& kmax = field.shape[2]-1;

            for(int b=0; b<_motion.size(); b++){

                // If the boundary motion is zero in all components, don't 
                // bother going any further.
                if (_motion[b] == 0.0) continue;

                // Since the input defines the boundary velocity as a
                // free-stream average, we invert the average here to get
                // the boundary region velocity.
                auto vel = 2.0*_motion[b];

                switch(b){

                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // West & East boundaries: Only the YZ components of the 
                    // boundary motion are considered.
                    case WEST:
                        for(int k=0; k<=kmax; k++){
                            for(int j=0; j<=jmax; j++){

                                UVW(0,j,k).v() = vel[1] - UVW(1,j,k).v();
                                UVW(0,j,k).w() = vel[2] - UVW(1,j,k).w(); 
                            }
                        }  
                        break;

                    case EAST:
                        for(int k=0; k<=kmax; k++){
                            for(int j=0; j<=jmax; j++){

                                UVW(imax,j,k).v() = vel[1] - UVW(imax-1,j,k).v();
                                UVW(imax,j,k).w() = vel[2] - UVW(imax-1,j,k).w();
                            }
                        }
                        break;

                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // South & North boundaries: only the XZ components of the
                    // boundary motion are considered.
                    case SOUTH:
                        for(int k=0; k<=kmax; k++){
                            for(int i=0; i<=imax; i++){

                                UVW(i,0,k).u() = vel[0] - UVW(i,1,k).u();
                                UVW(i,0,k).w() = vel[2] - UVW(i,1,k).w();
                            }
                        }
                        break;

                    case NORTH:
                        for(int k=0; k<=kmax; k++){
                            for(int i=0; i<=imax; i++){
                                
                                UVW(i,jmax,k).u() = vel[0] - UVW(i,jmax-1,k).u();
                                UVW(i,jmax,k).w() = vel[2] - UVW(i,jmax-1,k).w();
                            }
                        }
                        break;

                    // - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    // Down & Up boundaries: only the XY components of the
                    // boundary motion are considered.
                    case DOWN:
                        for(int j=0; j<=jmax; j++){
                            for(int i=0; i<=imax; i++){

                                UVW(i,j,0).u() = vel[0] - UVW(i,j,1).u();
                                UVW(i,j,0).v() = vel[1] - UVW(i,j,1).v();
                            }   
                        }
                        break;

                    case UP:
                        for(int j=0; j<=jmax; j++){
                            for(int i=0; i<=imax; i++){

                                UVW(i,j,kmax).u() = vel[0] - UVW(i,j,kmax-1).u();
                                UVW(i,j,kmax).v() = vel[1] - UVW(i,j,kmax-1).v();
                            }   
                        }
                        break;

                    default:
                        throw std::invalid_argument("Invalid boundary direction.");
                }
            }
        }
};

#endif