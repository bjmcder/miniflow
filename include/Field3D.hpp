#ifndef __FIELD_3D_HPP
#define __FIELD_3D_HPP

#include <array>
#include <vector>

template<typename T>
class Vec3{

    private:
        std::array<T, 3> _data;

    public:

        /**
         * 
        */
        Vec3(){
            _data = {0,0,0};
        }

        /**
         * 
        */
        T& x(){
            return _data[0];
        }

        /**
         * 
        */
        T& y(){
            return _data[1];
        }

        /**
         * 
        */
        T& z(){
            return _data[2];
        }

        /**
        * 
        */
        T& operator(int i){
            return _data[i];
        }

        /**
         * 
        */
        T& at(int i){
            return _data.at(i);
        }
};

template<typename T>
class Field3{

    private:
        std::vector<T>  _data;
        std::array<T,3> _strides;
        size_t          _size;
 
    public:

    /**
     * 
    */
    Field3(){}

    /**
     * 
    */
    Field3(int i=1, int j=1, int k=1){
        _strides = {i,j,k};
        _size = i*j*k;
        _data.resize(_size);
        this->zeros();
    }

    /**
     * Fill an array with zero-valued elements.
    */
    inline void zeros(){
        std::fill(_data.begin(), _data.end(), 0);
    }

    /**
     * Fill an array with one-valued elements.
    */
    inline void ones(){
        std::fill(_data.begin(), _data.end(), 1);
    }

    /**
     * Fill an array with elements of a specified value.
    */
    inline void fill(const T& val){
        std::fill(_data.begin(), _data.end(), val);
    }

    /**
     * 
    */
    inline int get_index(int i, int j, int k){
        return i + j*_strides[1] + k*_strides[2]*_strides[1];
    }

    /**
     * 
    */
    inline int get_index(const std::vector<int>& vals){
        return vals[0] + vals[1]*_strides[1] + vals[2]*_strides[2]*_strides[1];
    }

    /**
     * 
    */
    inline T& at(const std::vector<int>& vals){
        auto idx = get_index(vals);
        return _data.at(idx);
    }

    /**
     * 
    */
    inline T& at(int i, int j, int k){
        auto idx = get_index(i,j,k);
        return _data.at(idx);
    }

    /**
     * 
    */
    inline T& operator()(int i, int j, int k){
        auto idx = get_index(i,j,k);
        return _data[idx];
    }

    /**
     * 
    */
    inline T& operator()(const std::vector<int>& vals){
        auto idx = get_index(vals);
        return _data[idx];
    }

    /**
     * 
    */
    void print_row(int i, int j, int k){

        auto ncols = this->shape()[0];

        auto sep = "    ";
        auto prec = 10;
        auto width = 12;

        std::cout << "[";
        for (int i=0; i<ncols; i++){
            std::cout << 
            std::setw(width) <<
            std::setprecision(prec) <<
            this->at(_indexer.get_index(i, j, k)) <<
            sep;
        }

        std::cout << "]\n";
    }

    /**
     * 
    */
    void print_panel(int i, int j, int k){

        auto nrows = this->shape()[1];

        for (int i=0; i<nrows; i++){
            print_row(i, j, k);
        }
        std::cout << "\n";
    }

    /**
     * 
    */
    void print_field2d(int i, int j, int k){

        auto nrows = this->shape()[1];

        for (int i=nrows-1; i>=0; i--){
            print_row(i, j, k);
        }
        std::cout << "\n";
    }

    /**
     * 
    */
    void print_block(int i, int j, int k){

        auto nlayers = this->shape()[2];

        for (int i=0; i<nlayers; i++){
            print_panel(i, j, k);
        }
        std::cout << "\n";
    }

    /**
     * 
    */
    template<typename ...Args>
    void print_field3d(int i, int j, int k){

        auto nlayers = this->shape()[2];

        for (int i=nlayers-1; i>=0; i--){
            print_field2d(i, j, k);
        }

        std::cout << "\n";
    }
     
};

#endif