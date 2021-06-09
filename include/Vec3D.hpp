#ifndef VEC3D_HPP
#define VEC3D_HPP

#include <array>
#include <iostream>
#include <vector>

template<typename T>
class Vec3{

    private:
        std::array<T, 3> _data;
    public:

        /**
         * Default constructor: Initialize all components to zero.
        */
        Vec3(): _data({0,0,0}){}

        /**
         * Constructor from 3 individual components
        */
        Vec3(const T u, const T v=0, const T w=0): _data({u,v,w}){}

        /**
         * Returns the x-component.
        */
        T& x(){return _data[0];}

        /**
         * Returns the y-component.
        */
        T& y(){return _data[1];}

        /**
         * Returns the x-component.
        */
        T& z(){return _data[2];}

        /**
         * Returns the x-component (Alternate naming convention).
        */
        T& u(){return _data[0];}

        /**
         * Returns the y-component (Alternate naming convention).
        */
        T& v(){return _data[1];}

        /**
         * Returns the x-component (Alternate naming convention).
        */
        T& w(){return _data[2];}

        /**
         * Index into the vector directly. 
        */
        T& operator()(int i){return _data[i];}

        /**
         * Index into the vector directly (read only).
        */
        T operator[](int i) const {return _data[i];};

        /**
         * Index into the vector directly (with bounds checking).
        */
        T& at(int i){return _data.at(i);}

        /**
         * Return a reference to the raw data container.
        */
        std::array<T, 3>& data(){return _data;}

        /**
         * Assignment operator for scalars. Sets all components to equal
         * a single scalar value.
        */
        Vec3<T>& operator=(T scalar){
            _data[0] = scalar;
            _data[1] = scalar;
            _data[2] = scalar;

            return (*this);
        }

        /**
         * Assignment operator for STL vectors. Sets components equal to
         * each element of the input vector.
        */
        Vec3<T>& operator=(std::vector<T>& vec){

            if(vec.size() != 3){
                throw std::length_error("Input vector size "
                                        "must be equal to 3\n");
            }

            _data[0] = vec[0];
            _data[1] = vec[1];
            _data[2] = vec[2];

            return (*this);
        }

        /**
         * Assignment operator for STL arrays. Sets components equal to
         * each element of the input array.
        */
        Vec3<T>& operator=(std::array<T,3>& vec){

            _data[0] = vec[0];
            _data[1] = vec[1];
            _data[2] = vec[2];

            return (*this);
        }

        /**
         * Assignment operator for C-style arrays. Sets components equal to
         * each element of the input array.
        */
        Vec3<T>& operator=(T* vec){

            _data[0] = vec[0];
            _data[1] = vec[1];
            _data[2] = vec[2];

            return (*this);
        }

        /**
         * Increment operator. Increment this vector by the element-wise
         * values of another vector.
        */
        Vec3<T>& operator+=(Vec3<T> val){
            _data[0] += val[0];
            _data[1] += val[1];
            _data[2] += val[2];
            return (*this);
        }

        /**
         * Decrement operator. Decrement this vector by the element-wise
         * values of another vector.
        */
        Vec3<T>& operator-=(Vec3<T> val){
            _data[0] -= val[0];
            _data[1] -= val[1];
            _data[2] -= val[2];
            return (*this);
        }

        /**
         * Multiplicative increment operator. Increment this vector by
         * multiplying th element-wise values of another vector.
        */
        Vec3<T>& operator*=(Vec3<T> val){
            _data[0] *= val[0];
            _data[1] *= val[1];
            _data[2] *= val[2];
            return (*this);
        }

        /**
         * Divide decrement operator. Decrement this vector by dividing
         * the element-wise values of another vector.
        */
        Vec3<T>& operator/=(Vec3<T> val){
            _data[0] /= val[0];
            _data[1] /= val[1];
            _data[2] /= val[2];
            return (*this);
        }

        /**
         * Unary negation operator. Negate each component of the vector.
        */
        Vec3<T>& operator-(){
            _data[0] = -_data[0];
            _data[1] = -_data[1];
            _data[2] = -_data[2];
            return (*this);
        }

        /**
         * Return the sum of the vector's components.
        */
        T sum(){ return _data[0] + _data[1] + _data[2];}

        /**
         * Return the sum of the squares of the vector components.
        */
        T sumsquares(){
            return _data[0]*_data[0] + _data[1]*_data[1] + _data[2]*_data[2];
        }

        /**
         * Return the L2 norm (Euclidean length) of this vector.
        */
        T norm(){return sqrt(sumsquares());}
};

/**
 * Add two vectors together element-wise.
*/
template <typename T>
Vec3<T> operator+(const Vec3<T>& a, const Vec3<T>& b){
    return Vec3<T>(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

/**
 * Subtract two vectors element-wise.
*/
template <typename T>
Vec3<T> operator-(const Vec3<T>& a, const Vec3<T>& b){
    return Vec3<T>(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}

/**
 * Multiply two vectors together element-wise.
*/
template <typename T>
Vec3<T> operator*(const Vec3<T>& a, const Vec3<T>& b){
    return Vec3<T>(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}

/**
 * Divide two vectors element-wise.
*/
template <typename T>
Vec3<T> operator/(const Vec3<T>& a, const Vec3<T>& b){
    return Vec3<T>(a[0]/b[0], a[1]/b[1], a[2]/b[2]);
}

/**
 * Compute the dot-product of two vectors.
*/
template <typename T>
Vec3<T> dot(const Vec3<T>& a, const Vec3<T>& b){
    return sum(a*b);
}

/**
 * Compute the cross-product of two vectors.
*/
template <typename T>
Vec3<T> cross(const Vec3<T>& a, const Vec3<T>& b){
    auto cx = a.y()*b.z() - a.z()*b.y();
    auto cy = a.z()*b.x() - a.x()*b.z();
    auto cz = a.x()*b.y() - a.y()*b.x();

    return Vec3<T>(cx, cy, cz);
}

/**
 * Scalar multiplication.
*/
template <typename T>
Vec3<T> operator*(const Vec3<T>& a, T val){
    return Vec3<T>(val*a[0], val*a[1], val*a[2]);
}

/**
 * Scalar multiplication (reversed argument order).
*/
template <typename T>
Vec3<T> operator*(T val, const Vec3<T>& a){
    return Vec3<T>(val*a[0], val*a[1], val*a[2]);
}

/**
 * Scalar division. Onlt this argument order makes mathematical sense.
*/
template <typename T>
Vec3<T> operator/(const Vec3<T>& a, T val){
    return Vec3<T>(a[0]/val, a[1]/val, a[2]/val);
}

/**
 * Outstream operator.
*/
template <typename T>
std::ostream& operator<<(std::ostream& os, Vec3<T>& vec){
    os << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
    return os;
}

/**
 * Equality operator on scalars. Checks that all elements are equal to a 
 * single value.
*/
template <typename T>
bool operator==(const Vec3<T>& a, T val){

    for(int i=0; i<3; i++){
        if (a[i] == val) continue;
        else return false;
    }
    return true;
}

/**
 * Equality operator on vector. Checks that all elements are equal to a 
 * single value.
*/
template <typename T>
bool operator==(const Vec3<T>& a, const Vec3<T>& b){

    for(int i=0; i<3; i++){
        if (a[i] == b[i]) continue;
        else return false;
    }
    return true;
}

#endif