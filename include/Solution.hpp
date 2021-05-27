#ifndef __SOLUTION_HPP
#define __SOLUTION_HPP

#include <array>
#include <cmath>
#include <vector>

#include "Field3D.hpp"
#include "Problem.hpp"

template<typename T>
class Field3; 

template<typename T>
class Problem;

template<typename T>
struct Solution{

    std::vector<int> shape;

    Field3<T> U;
    Field3<T> V;
    Field3<T> W;

    Field3<T> pressure;

    Field3<T> F;
    Field3<T> G;
    Field3<T> H;
    Field3<T> rhs;

    //Field3<T> vorticity;
    //Field3<T> temperature;

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

        U = Field3<T>(shape);
        V = Field3<T>(shape);
        W = Field3<T>(shape);

        F = Field3<T>(shape);
        G = Field3<T>(shape);
        H = Field3<T>(shape);

        pressure = Field3<T>(shape);
        rhs = Field3<T>(shape);
    }

    /**
     * 
    */
    void initialize_solution(Problem<T>& problem){
        
        U.fill(problem.flow_parameters().initial_velocities()[0]);
        V.fill(problem.flow_parameters().initial_velocities()[1]);
        W.fill(problem.flow_parameters().initial_velocities()[2]);

        F.zeros();
        G.zeros();
        H.zeros();

        rhs.zeros();
        pressure.zeros();
    }

    /**
     * 
    */
    std::array<T,3> max_velocity_components(){
        std::array<T,3> max_components;

        max_components[0] = absmax(U.data());
        max_components[1] = absmax(V.data());
        max_components[2] = absmax(W.data());

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
    Field3<T>& velocity_component(int dim){
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
    Field3<T>& tentative_momentum(int dim){
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