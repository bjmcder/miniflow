#ifndef __GEOMETRY_HPP
#define __GEOMETRY_HPP

#include <vector>

#include "NDArray.hpp"

template<typename T>
class Geometry{
    protected:

        int _dimension;
        std::vector<T> _length;
        std::vector<int> _ncells;
        std::vector<T> _cellsize;

        /**
         * Compute the mesh cell size. Used in the class constructor.
        */
        inline void _compute_cellsize(){
            for (int i=0; i<_dimension; i++){
                _cellsize[i] = _length[i]/((T)_ncells[i]);
            }
        }

        /**
         * Setup the array sizes according to the problem dimension.
        */
        inline void _size_arrays(){
            _length.resize(_dimension);
            _ncells.resize(_dimension);
            _cellsize.resize(_dimension);
        }

    public:

        /**
         * Default constructor (does nothing)
        */
        Geometry(){return;}

        /**
         * Return the cell size in the x-dimension.
        */
        inline T dx(){
            return _cellsize[0];
        }

        /**
         * Return the cell count in the x-dimension.
        */
        inline size_t nx(){
            return _ncells[0];
        }

        inline std::vector<int> ncells(){
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

        inline std::vector<T>& cell_sizes(){
            return _cellsize;
        }

        /**
         * Return the dimension of the geometry.
        */
        inline size_t dimension(){return _dimension;}
    
};

/**
 * A class for holding the problem geometry definition.
 * 
 * Template Parameters
 * -------------------
 * T : Type
 *  A numeric type (e.g. double) for defining lengths and cell sizes. 
 * 
*/
template <typename T>
class Geometry1D : public Geometry<T>{
    
    public:

        /**
         * Preferred Constructor (1D).
        */
        Geometry1D(const T& Lx, const int& nx){
            
            Geometry<T>::_dimension = 1;
            Geometry<T>::_size_arrays();

            Geometry<T>::_length = {Lx};
            Geometry<T>::_ncells = {nx};

            Geometry<T>::_compute_cellsize();
        }
};

template<typename T>
class Geometry2D : public Geometry<T>{

    public:

    /**
     * Preferred Constructor (2D).
    */
    Geometry2D(const T& Lx, 
               const T& Ly,
               const int& nx,
               const int& ny){

        Geometry<T>::_dimension = 2;
        Geometry<T>::_size_arrays();

        Geometry<T>::_length = {Lx, Ly}; 
        Geometry<T>::_ncells = {nx, ny};

        Geometry<T>::_compute_cellsize();
    }

    /**
     * Return the cell size in the y-dimension. (Dim > 1).
    */
    inline T dy(){
        return Geometry<T>::_cellsize[1];
    }

    /**
     * Return the cell count in the y-dimension. (Dim > 1).
    */
    inline size_t ny(){
        return Geometry<T>::_ncells[1];
    }
};

template<typename T>
class Geometry3D : public Geometry<T>{

    public:

        /**
         * Preferred Constructor (3D).
        */
        Geometry3D(const T& Lx, 
                   const T& Ly,
                   const T& Lz,
                   const int& nx,
                   const int& ny,
                   const int& nz){
                           
            Geometry<T>::_dimension = 3;
            Geometry<T>::_size_arrays();

            Geometry<T>::_length = {Lx, Ly, Lz};
            Geometry<T>::_ncells = {nx, ny, nz};

            Geometry<T>::_compute_cellsize(); 
        }

        /**
         * Return the cell size in the y-dimension. (Dim > 1).
        */
        inline T dy(){
            return Geometry<T>::_cellsize[1];
        }

        /**
         * Return the cell size in the z-dimension (Dim > 2).
        */
        inline T dz(){
            return Geometry<T>::_cellsize[2];
        }

        /**
         * Return the cell count in the y-dimension. (Dim > 1).
        */
        inline size_t ny(){
            return Geometry<T>::_ncells[1];
        }

        /**
         * Return the cell count in the z-dimension (Dim > 2).
        */
        inline size_t nz(){
            return Geometry<T>::_ncells[2];
        }
};

#endif