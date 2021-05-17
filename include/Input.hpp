#ifndef __INPUT_HPP
#define __INPUT_HPP

#include <string>
#include <memory>

#include "ext/toml.h"
#include "BoundaryConditions.hpp"
#include "FlowParameters.hpp"
#include "Geometry.hpp"
#include "Problem.hpp"
#include "TimeStepper.hpp"

template<typename T>
class Input{

    private:

        toml::Value _toml_dat;
        int _dim;

    public:

    /**
     * Default constructor (does nothing).
    */
    Input(){}

    /**
     * Preferred constructor, load a TOML file and process the contents in
     * the constructed object.
     * 
     * Parameters
     * ----------
     * path : std::string
     *  Path to a TOML file. 
    */
    Input(const std::string& path){
        load_toml(path);
        if (this->valid()){
            this->set_dimension();
        }
        else{
            throw std::invalid_argument("Input validation failed. Check input.\n");
        }
    }

    /**
     * Load a TOML file from a specified path.
    */
    void load_toml(const std::string& path){
        auto raw_parsed = toml::parseFile(path);

        if(raw_parsed.valid()){
            _toml_dat = raw_parsed.value;
        }
    }

    /**
     * Validate the file by checking the minimum runnable parameters are
     * included.
    */
    bool valid(){

       // Geometry
        if(!valid_geometry()) return false;

       // Time

       // Solver

       // Problem

       return true;
    }

    /**
     * Validate the geometry section in the input file.
    */
    bool valid_geometry(){

        std::string section = "geometry";

        // Verify the main section exists
        if (!_toml_dat.has(section)){
            std::cerr << "Input file missing [geometry] section.\n";
            return false;
        }
        
        // Verify the required subsections exist
        int subsections_ok = 0;

        subsections_ok += _toml_dat.has(section + ".dimension");
        subsections_ok += _toml_dat.has(section + ".domain_size");
        subsections_ok += _toml_dat.has(section + ".num_cells");

        if(subsections_ok < 3) {
            std::cerr << "Input file missing one or more required" \
            " subsections in [geometry]: (dimension, domain_size, num_cells)."
            << std::endl;
            return false;
        }

        // Verify the subsection dimensions match
        auto dim = _toml_dat.get<int>("geometry.dimension");
        bool ok;

        ok = _toml_dat.get<toml::Array>("geometry.domain_size").size() == dim;
        if(!ok) {
            std::cerr << "[geometry] domain_size dimensions do not match" \
            " (dimension = " << dim << ")\n";
            return false;
        }

        ok = _toml_dat.get<toml::Array>("geometry.num_cells").size() == dim;
        if(!ok) {
            std::cerr << "[geometry] num_cells dimensions do not match" \
            " (dimension = " << dim << ")\n";
            return false;
        }
        return true;
    }

    /**
     * Set the geometric dimenstion based on the parameter defined in
     * the input file.
    */
    void set_dimension(){
        _dim = _toml_dat.get<int>("geometry.dimension");
    }

    /**
     * Return the geometric dimension of the problem.
    */
    int dimension(){
        return _dim;
    }

    /**
     * 
    */
    auto build_geometry1d(){
        auto Lx = \
            _toml_dat.get<toml::Array>("geometry.domain_size")[0].as<T>();
        auto nx = \
            _toml_dat.get<toml::Array>("geometry.num_cells")[0].as<int>();

        return Geometry1D<T>(Lx, nx);
    }

    /**
     * 
    */
    auto build_geometry2d(){
        auto Lx = \
            _toml_dat.get<toml::Array>("geometry.domain_size")[0].as<T>();
        auto Ly = \
            _toml_dat.get<toml::Array>("geometry.domain_size")[1].as<T>();

        auto nx = \
            _toml_dat.get<toml::Array>("geometry.num_cells")[0].as<int>();
        auto ny = \
            _toml_dat.get<toml::Array>("geometry.num_cells")[1].as<int>();

        return Geometry2D<T>(Lx, Ly, nx, ny);
    } 

    /**
     * 
    */
    auto build_geometry3d(){
        auto Lx = \
            _toml_dat.get<toml::Array>("geometry.domain_size")[0].as<T>();
        auto Ly = \
            _toml_dat.get<toml::Array>("geometry.domain_size")[1].as<T>();
        auto Lz = \
            _toml_dat.get<toml::Array>("geometry.domain_size")[2].as<T>();

        auto nx = \
            _toml_dat.get<toml::Array>("geometry.num_cells")[0].as<int>();
        auto ny = \
            _toml_dat.get<toml::Array>("geometry.num_cells")[1].as<int>();
        auto nz = \
            _toml_dat.get<toml::Array>("geometry.num_cells")[2].as<int>();

        return Geometry3D<T>(Lx, Ly, Lz, nx, ny, nz);
    }

    /**
     * 
    */
    Geometry<T> build_geometry(){
        switch (_dim){
        case 1:
            return build_geometry1d();
            break;
        case 2:
            return build_geometry2d();
            break;
        case 3:
            return build_geometry3d();
            break;
        default:
            throw std::invalid_argument("Could not build geometry");
        }
    }

    /**
     * 
    */
    TimeStepper<T> build_timestepper(){

        auto max_time = _toml_dat.get<T>("time.max_time");
        auto tau = _toml_dat.get<T>("time.step_factor");
        auto dt0 = 0;
        if (_toml_dat.has("time.dt_start")){
            dt0 = _toml_dat.get<T>("time.dt_start");
        }

        return TimeStepper<T>(max_time, tau, dt0);
    }


    /**
     * 
    */
    BoundaryConditions<T> build_boundary_conditions(){

        auto conditions = _toml_dat.get<toml::Table>("boundary.conditions");
        auto motion = _toml_dat.get<toml::Table>("boundary.motion");

        std::cout << conditions << std::endl;
        std::cout << motion << std::endl;

        // We don't loop here because the ordering of the boundaries in the
        // input is not guaranteed. We enforce the ordering explicitly
        // by populating the vector in our ordering convention.
        //
        // Enforcing the naming conventions would be straightforward with
        // regular expressions, but TinyTOML doesn't support that (yet).
        std::vector<std::string> string_bcs;

        string_bcs.push_back(conditions["west"].as<std::string>());
        string_bcs.push_back(conditions["east"].as<std::string>());

        if (_dim > 1){
            string_bcs.push_back(conditions["south"].as<std::string>());
            string_bcs.push_back(conditions["north"].as<std::string>());
        }
        if (_dim > 2){ 
            string_bcs.push_back(conditions["down"].as<std::string>());
            string_bcs.push_back(conditions["up"].as<std::string>());
        }
        
        // Convert string bcs to integers
        std::vector<int> bcs;
        for (const auto& b: string_bcs){
            bcs.push_back(0);
        }

        // Similarly, with boundary motion values, we need to enforce the
        // The ordering.
        std::vector<T> motion_vals;

        T mval = 0.0;

        // If we have the value defined in the input, add it to the motion_vals
        // vector, otherwise, add a zero.
        mval = !motion["west"].empty() ? motion["west"].as<T>() : 0.0;
        motion_vals.push_back(mval);

        mval = !motion["east"].empty() ? motion["east"].as<T>() : 0.0;
        motion_vals.push_back(mval);

        if(_dim > 1){
            mval = !motion["south"].empty() ? motion["south"].as<T>() : 0.0;
            motion_vals.push_back(mval);

            mval = !motion["north"].empty() ? motion["north"].as<T>() : 0.0;
            motion_vals.push_back(mval);   
        }

        if(_dim > 2){
            mval = !motion["down"].empty() ? motion["down"].as<T>() : 0.0;
            motion_vals.push_back(mval);

            mval = !motion["up"].empty() ? motion["up"].as<T>() : 0.0;
            motion_vals.push_back(mval);      
        }

        return BoundaryConditions<T>(bcs, motion_vals);
    }

    /**
     * 
    */
    FlowParameters<T> build_flow_params(){
        auto reynolds = _toml_dat.get<T>("flow.reynolds");
        auto initial_velocities = \
            _toml_dat.get<toml::Array>("flow.initial_velocity");
        auto body_forces = \
            _toml_dat.get<toml::Array>("flow.body_force");

        return FlowParameters<T>(_dim);
    }

    /**
     * Print the raw TOML input to stdout.
    */
    void print_raw(){
        std::cout << _toml_dat << std::endl;
    }
    
};

#endif