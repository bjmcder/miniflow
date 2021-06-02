#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <array>
#include <numeric>
#include <vector>

#include "Types.hpp"

template<typename T>
class Geometry{

    private:

        int _dimension = 3;

        std::array<T,3> _length;
        std::array<int,3> _ncells;
        std::array<T,3> _cellsize;

        /**
         * Compute the mesh cell size. Used in the class constructor.
        */
        inline void _compute_cellsize(){
            for (int i=0; i<_dimension; i++){
                _cellsize[i] = _length[i]/((T)_ncells[i]);
            }
        }

    public:

        /**
         * Default constructor (does nothing)
        */
        Geometry(){}

        /**
         * Preferred Constructor (3D).
        */
        Geometry(const T& Lx, 
                 const T& Ly,
                 const T& Lz,
                 const int& nx,
                 const int& ny,
                 const int& nz){

            _length = {Lx, Ly, Lz};
            _ncells = {nx, ny, nz};

            _compute_cellsize(); 
        }

        /**
         * Return the cell size in the x-dimension.
        */
        inline T dx(){
            return _cellsize[0];
        }

        /**
         * Return the cell size in the y-dimension.
        */
        inline T dy(){
            return _cellsize[1];
        }

        /**
         * Return the cell size in the z-dimension.
        */
        inline T dz(){
            return _cellsize[2];
        }

        /**
         * Return the cell count in the x-dimension.
        */
        inline size_t nx(){
            return _ncells[0];
        }

        /**
         * Return the cell count in the y-dimension.
        */
        inline size_t ny(){
            return _ncells[1];
        }

        /**
         * Return the cell count in the x-dimension.
        */
        inline size_t nz(){
            return _ncells[2];
        }

        /**
         * 
        */
        inline std::array<int, 3> ncells(){
            return _ncells;
        }

        /**
         * Return the total number of mesh cells in the problem.
        */
        inline size_t total_cells() {
            auto val = std::accumulate(_ncells.begin(),
                                       _ncells.end(),
                                       0,
                                       std::plus<>());
            return (size_t)val;
        }

        inline std::array<T, 3>& cell_sizes(){
            return _cellsize;
        }

        /**
         * Return the dimension of the geometry.
        */
        inline size_t dimension(){return _dimension;}
};


#endif