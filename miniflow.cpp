#include <iostream>

#include "include/ext/toml.h"
#include "include/BoundaryConditions.hpp"
#include "include/Input.hpp"
#include "include/Geometry.hpp"
#include "include/Types.hpp"
#include "include/Solver.hpp"
#include "include/Problem.hpp"

int main(int argc, char** argv){

    if (argc < 2){
        std::cout << "Usage: miniflow <input-file>.toml" << std::endl;
        return 1; 
    }

    // Load the input TOML file
    auto indat = Input<real_t>(argv[1]);

    // Build the problem
    auto geom = indat.build_geometry();
    auto tstepper = indat.build_timestepper();
    auto bcs = indat.build_boundary_conditions();
    auto flow_params = indat.build_flow_params();

    auto problem = Problem<real_t>(geom, tstepper, bcs, flow_params);

    //auto solver = Solver<real_t>(problem);

    // Solve the problem
    //solver.solve();

    
    return 0;
}