#include <iostream>

#include "include/Problem.hpp"
#include "include/ext/toml.h"
#include "include/BoundaryConditions.hpp"
#include "include/Input.hpp"
#include "include/OutputSettings.hpp"
#include "include/Geometry.hpp"
#include "include/Types.hpp"
#include "include/Solver.hpp"


int main(int argc, char** argv){

    if (argc < 2){
        std::cout << "Usage: miniflow <input-file>.toml" << std::endl;
        return 1; 
    }

    std::cout << std::setprecision(12);

    // Load the input TOML file
    auto indat = Input<scalar_t>(argv[1]);

    // Build the problem
    auto geom = indat.build_geometry();
    auto tstepper = indat.build_timestepper();
    auto bcs = indat.build_boundary_conditions();
    auto flow_params = indat.build_flow_params();
    auto solver_settings = indat.build_solver_settings();
    auto output_settings = indat.build_output_settings();


    auto problem = Problem<scalar_t>(geom, tstepper, bcs, flow_params);

    auto solver = Solver<scalar_t>(problem, solver_settings, output_settings);

    solver.solve();
    
    return 0;
}