#ifndef VTK_OUTPUT_HPP
#define VTK_OUTPUT_HPP

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include "ext/pugixml.hpp"
#include "ext/base64.hpp"
#include "Types.hpp"
#include "Problem.hpp"
#include "Solution.hpp"
#include "OutputSettings.hpp"

template<typename T>
class VTKFile{

    private:
        
        pugi::xml_document _vtk_doc;

        std::string _write_path;
        std::string _vtk_type;
        OutputSettings _settings;

    public:

        /**
         * Default constructor. Construct a blank VTI file structure.
        */
        VTKFile(OutputSettings settings,
                std::string fpath="miniflow.vti",
                std::string vtk_type="ImageData"){

            _settings = settings;
            _write_path = fpath;
            _vtk_type = vtk_type;

            create_document();
        }

        /**
         * 
        */
        void create_document(){
            
            // Create the XML declaration
            create_xml_declaration();

            // Create the root VTK section
            create_vtk_section(_vtk_type);

            // Create the dataset section
            create_dataset(_vtk_type);
        }

        /**
         * Set the extent attributes in the VTK file from the problem
         * geometry parameters.
        */
        void set_geometry(Problem<T>& problem){

            auto& geom = problem.geometry();

            // Set the global dataset attributes
            auto vtk_root = _vtk_doc.child("VTKFile");
            auto vti_node = vtk_root.child("ImageData");
            auto whole_extent = vti_node.attribute("WholeExtent");
            auto spacing = vti_node.attribute("Spacing");

            // Construct a string of the min and max extents and set the
            // attribute
            std::stringstream ss;

            ss << "0 " 
               << geom.nx()-1
               << " 0 " 
               << geom.ny()-1 
               << " 0 " 
               << geom.nz()-1;
            std::string extent_str = ss.str();

            whole_extent.set_value(extent_str.c_str());

            // Re-initialize the string and set the spacings in each direction
            ss = std::stringstream();

            ss << geom.dx() << " " << geom.dy() << " " << geom.dz();
            std::string spacing_str = ss.str();

            spacing.set_value(spacing_str.c_str());

            // Set the piece-level extents (same as global for a serial case)
            auto piece_node = vti_node.child("Piece");
            auto piece_extent = piece_node.attribute("Extent");
            piece_extent.set_value(extent_str.c_str());
        }

        /**
         * 
        */
        void store_velocity(Solution<T>& solution){

            // Get the raw solution vector. Paraview can infer the strides
            // from the dataset attributes.
            auto& vel = solution.velocity.data();

            int dim = solution.shape.size();

            // Create the dataset and attributes
            auto vtk_root = _vtk_doc.child("VTKFile");
            auto vti_node = vtk_root.child("ImageData");
            auto piece_node = vti_node.child("Piece");
            auto pdata_node = piece_node.child("PointData");
            auto data_array = pdata_node.append_child("DataArray");

            data_array.append_attribute("type") = "Float64";
            data_array.append_attribute("Name") = "Velocity";
            data_array.append_attribute("NumberOfComponents") = 3;
            
            auto format = _settings.format;
            data_array.append_attribute("format") = format.c_str();

            std::string vtk_buffer("");

            if (format == "binary"){
                // We need to prepend the data sizes in bytes to the main dataset
                size_t data_size = sizeof(vector3d_t)*vel.size();
                size_t total_size = data_size + sizeof(data_size); 

                // Now, create a buffer of bytes to encode...
                auto raw_buffer = \
                    reinterpret_cast<uint8_t*>(malloc(total_size));

                // ... prepend the size header
                auto loc = reinterpret_cast<const uint8_t*>(&data_size);
                std::copy(loc,
                        loc + sizeof(data_size),
                        raw_buffer);

                // ... append the raw data
                memcpy(raw_buffer + sizeof(data_size),
                    (uint8_t*)vel.data(),
                    data_size);

                // Encode and append the solution vector data in base64 to
                // support the XML format
                vtk_buffer = base64_encode(raw_buffer, total_size);

                // Free the byte array so we don't leak memory
                free(raw_buffer);
            }

            
            else if(format == "ascii"){
                std::stringstream ss;
                
                for(int i=0;i<vel.size(); i++){
                    ss << vel[i] << " ";
                }

                vtk_buffer = ss.str();
            }
            else{
                std::cerr << "Warning: Output dataset format not "
                          << "recognized. Data will not be written.\n";
            }

            // Add the encoded dataset to the buffer.
            auto rawdat = data_array.append_child(pugi::node_pcdata);
            rawdat.set_value(vtk_buffer.c_str());
        }

        /**
         * 
        */
        void create_xml_declaration(std::string version="1.0", 
                                    std::string encoding="utf-8"){

            auto declaration_node = \
                _vtk_doc.append_child(pugi::node_declaration);

            declaration_node.append_attribute("version") = version.c_str();
            declaration_node.append_attribute("encoding") = encoding.c_str();
        }

        /**
         * 
        */
        void create_vtk_section(std::string type="ImageData",
                                std::string version="1.0",
                                std::string byte_order="LittleEndian",
                                std::string compressor="",
                                std::string header_type="UInt64"){
            
            auto vtk_root = _vtk_doc.append_child("VTKFile");

            vtk_root.append_attribute("type") = type.c_str();
            vtk_root.append_attribute("version") = version.c_str();
            vtk_root.append_attribute("byte_order") = byte_order.c_str();
            vtk_root.append_attribute("header_type") = header_type.c_str();
            vtk_root.append_attribute("compressor") = compressor.c_str();
            
        }

        /**
         * 
        */
        void create_dataset(std::string vtk_type){

            // ImageData
            if (vtk_type == "ImageData"){
                create_vti_dataset();
            }
            else{
                std::string err_str = "Dataset type ";
                err_str += vtk_type;
                err_str += "is not currently implemented.\n";

                throw std::invalid_argument(err_str);
            }
        }

        /**
         * 
        */
        void create_vti_dataset(){

            // Get the VTK root section
            auto vtk_root = _vtk_doc.child("VTKFile");

            // Append an ImageData section

            auto vti_sec = vtk_root.append_child("ImageData");
            vti_sec.append_attribute("WholeExtent");
            vti_sec.append_attribute("Origin") = "0.0, 0.0, 0.0";
            vti_sec.append_attribute("Spacing");

            // Create the empty VTI structure
            auto vti_piece = vti_sec.append_child("Piece");
            vti_piece.append_attribute("Extent");

            vti_piece.append_child("PointData");
            vti_piece.append_child("CellData");
            
        }

        /**
         * 
        */
        void set_path(std::string path){
            _write_path = path;
        }

        /**
         * Save the data to a VTK file.
        */
        void save_file(){

            _vtk_doc.save_file(_write_path.c_str(), PUGIXML_TEXT("    "));
        }

};

#endif
