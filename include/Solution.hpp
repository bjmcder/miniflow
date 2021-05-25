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

    std::vector<int> shape;

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
     * Default constructor (does nothing).
    */
    Solution(){}

    /**
     * 
    */
    Solution(Problem<T>& problem){

        this->shape = problem.geometry().ncells();

        std::vector<size_t> arraysizes;
        arraysizes.resize(this->shape.size());

        // Need to convert int to size_t
        std::transform(this->shape.begin(), 
                       this->shape.end(),
                       arraysizes.begin(),
                        [](int x) { return (size_t)x;});

        set_array_sizes(arraysizes);

        initialize_solution(problem);
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
    void initialize_solution(Problem<T>& problem){

        std::fill(U.data().begin(),
                  U.data().end(), 
                  problem.flow_parameters().initial_velocities()[0]);

        F.zeros();
        rhs.zeros();
        pressure.zeros();

        if (shape.size() > 1){
            std::fill(V.data().begin(),
                  V.data().end(), 
                  problem.flow_parameters().initial_velocities()[1]);

            G.zeros();
        }

        if (shape.size() > 2){
            std::fill(W.data().begin(),
                  W.data().end(), 
                  problem.flow_parameters().initial_velocities()[2]);

            H.zeros();
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
    T absmax(std::vector<T> data) const{
        auto max = 0.0;
        for (const auto& elem: data){
            if (fabs(elem) > max){
                max = std::fabs(elem);
            }
        }
        return max;
    }

    /**
     * 
    */
    NDArray<T>& velocity_component(int dim){
        switch(dim){
            case 0:
                return U;
            case 1:
                return V;
            case 2:
                return W;
            default:
                throw std::invalid_argument("Invalid dimension.");
        }
    }

    /**
     * 
    */
    NDArray<T>& tentative_momentum(int dim){
        switch(dim){
            case 0:
                return F;
            case 1:
                return G;
            case 2:
                return H;
            default:
                throw std::invalid_argument("Invalid dimension.");
        }
    }
};

#endif