#ifndef OUTPUT_SETTINGS_HPP
#define OUTPUT_SETTINGS_HPP

#include <string>

struct OutputSettings{

    int write_every;
    std::string base_name;
    std::string format;

    /**
     * Default constructor.
    */
    OutputSettings(): write_every(0),
              base_name("tstep"),
              format("binary"){}

    /**
     * Preferred constructor. Construct from user-defined output settings.
    */
    OutputSettings(int writeevery, std::string basename, std::string fformat){
        write_every = writeevery;
        base_name = basename;
        format = fformat;
    }
};




#endif