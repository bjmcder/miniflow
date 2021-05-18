#ifndef __SOLUTION_HPP
#define __SOLUTION_HPP

#include <cmath>
#include <vector>

#include "NDArray.hpp"
#include "Problem.hpp"

template<typename T>
class NDArray; 

template<typename T>
class Problem;

template<typename T>
struct Solution{

    NDArray<T> U;
    NDArray<T> V;
    NDArray<T> W;

    NDArray<T> pressure;

    NDArray<T> F;
    NDArray<T> G;
    NDArray<T> H;
    NDArray<T> rhs;

    //NDArray<T> vorticity;
    //NDArray<T> temperature;

    /**
     * 
    */
    Solution(){}

    /**
     * 
    */
    Solution(Problem<T>& problem){

        auto shape = problem.geometry().ncells();

        std::vector<size_t> arraysizes;
        arraysizes.resize(shape.size());

        // Need to convert int to size_t
        std::transform(shape.begin(), 
                       shape.end(),
                       arraysizes.begin(),
                        [](int x) { return (size_t)x;});

        set_array_sizes(arraysizes);
    }

    /**
     * 
    */
    void set_array_sizes(const std::vector<size_t>& shape){

        U = NDArray<T>(shape);
        F = NDArray<T>(shape);
        pressure = NDArray<T>(shape);
        rhs = NDArray<T>(shape);

        if (shape.size()>1){
            V = NDArray<T>(shape);
            G = NDArray<T>(shape);
        }
        if (shape.size() > 2){
            W = NDArray<T>(shape);
            H = NDArray<T>(shape);
        }
    }

    /**
     * 
    */
    std::vector<T> max_velocity_components(int dim){
        std::vector<T> max_components;

        max_components.resize(dim);

        max_components[0] = absmax(U.data());

        if(dim>1){
            max_components[1] = absmax(V.data());
        }

        if(dim>2){
            max_components[2] = absmax(W.data());
        }

        return max_components;
    }

    /**
     * 
    */
    T& absmax(std::vector<T> data) const{
        auto max = 0.0;

        for (const auto& elem: data){
            if (fabs(elem) > max){
                max = std::fabs(elem);
            }
        }
    }
};

#endif