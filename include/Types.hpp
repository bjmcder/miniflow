#ifndef TYPES_HPP
#define TYPES_HPP

#include "Vec3D.hpp"
#include "Field3D.hpp"

// Defines the real number type used in the application.
// Usually float or double.
using scalar_t = double;

// Multi-dimensional numerical types are derived from the base scalar type.
using vector3d_t = Vec3<scalar_t>;
using vector_field_t = Field3<vector3d_t>;
using scalar_field_t = Field3<scalar_t>; 

#endif