#ifndef __SOLUTION_HPP
#define __SOLUTION_HPP

#include <vector>

#include "NDArray.hpp"
#include "BoundaryConditions.hpp"
#include "Problem.hpp"

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

    Solution(Problem<T>& problem){

        auto shape = problem.geometry().ncells();
        set_array_sizes(shape);
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