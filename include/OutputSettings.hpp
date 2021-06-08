#ifndef OUTPUT_SETTINGS_HPP
#define OUTPUT_SETTINGS_HPP

#include <string>

struct OutputSettings{

    int write_every;
    std::string base_name;
    std::string format;

    /**************************************************************************
     * Default constructor. Sets default values for Output behavior.
    **************************************************************************/
    OutputSettings(): write_every(0),
              base_name("tstep"),
              format("binary"){}

    /**************************************************************************
     * Preferred constructor. Construct from user-defined output settings.
     * 
     * Parameters
     * ----------
     * writeevery : int
     *  Interval (in timesteps) in which to write output files.
     * 
     * basename : std::string
     *  Base file name for output files. The timestep number is appended to
     *  the base name.
     * 
     * fformat : std::string
     *  The dataset format (ascii or binary.)
    **************************************************************************/
    OutputSettings(int writeevery, std::string basename, std::string fformat){
        write_every = writeevery;
        base_name = basename;
        format = fformat;
    }
};




#endif