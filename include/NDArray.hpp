#ifndef __ARRAY_HPP
#define __ARRAY_HPP

#include <iostream>
#include <iomanip>
#include <numeric>
#include <vector>


/**
 * A helper class for providing N-Dimensional indexing into a 1D array.
*/
class NDIndexer{

private:
    std::vector<size_t> _stride;

public:

    /**
     * Default constructor (does nothing)
    */
    NDIndexer(){return;}

    /**
     * Preferred constructor from vector of dimensions. 
     * 
     * Parameters
     * ----------
     * strides : std::vector<size_t>
     *  A vector of dimensional strides. Length equals the number of
     *  dimensions in the array. 
    */
    NDIndexer(std::vector<size_t> strides){
        _stride = strides;
    }

    /**
     * Compute an index in a 1-D array given an N-D set of indices.
     * 
     * Parameters
     * ----------
     * arg0 : 
    */
    template<typename ...Args>
    inline size_t get_index(size_t const & arg0, Args const & ... args){

        auto ix = std::vector<size_t>({(size_t)arg0, (size_t)args...});

        return get_index(ix);
    }

    /***/
    inline size_t get_index(const std::vector<size_t>& ix){

        size_t index = 0;

        index += ix[0];

        // Here we need to recursively compute the flattened index.
        // l*(sizek*sizej*sizei)+ k*(size_i*size_j) + j*(size_i) + i
        for(size_t i=1; i<ix.size(); i++){

            auto dprod = 1;

            for (size_t d=0; d<i; d++){
                dprod *= _stride[d];
            }

            index += ix[i]*dprod;
        } 
        return index;
    }
};

template<typename T>
class NDArray{

    public:

    std::vector<T> _data;
    std::vector<size_t> _dims;
    size_t _size;
    size_t _ndim;

    NDIndexer _indexer;

    public:

    NDArray(){}

    /**
     * Preferred constructor. Provide a variadic list of arguments specifying
     * the dimensions of the NDArray. There is no default fill value for the
     * array elements.
     * 
     * Parameters
     * ----------
     * 
    */
    template<typename ...Args>
    NDArray(size_t const & arg0, Args const & ... args){

        _dims = {(size_t)arg0, (size_t)args...};

        _size = std::accumulate(_dims.begin(),
                                _dims.end(),
                                1,
                                std::multiplies<double>());

        _ndim = _dims.size();

        _data.resize(_size);

        _indexer = NDIndexer(_dims);
    }

    /**
     * 
    */
    NDArray(const std::vector<size_t>& sizes){
        
        _dims = sizes;

        _size = std::accumulate(_dims.begin(),
                                _dims.end(),
                                1,
                                std::multiplies<double>());

        _ndim = _dims.size();

        _data.resize(_size);

        _indexer = NDIndexer(_dims);
    }

    /**
     * 
    */
    inline size_t size(){
        return _data.size();
    }

    inline std::vector<size_t> shape(){
        return _dims;
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
     * Indexing into the NDArray with bounds checking.
    */
    template<typename ...Args>
    inline T& at(size_t const & arg0, Args const & ... args){
        auto ix =  _indexer.get_index(arg0, args...);
        return _data.at(ix);
    }

    /**
     * Function call operator for convenient indexing. No bounds checking.
    */
    template<typename ...Args>
    inline T& operator()(size_t const & arg0, Args const & ... args){
        auto ix =  _indexer.get_index(arg0, args...);
        return _data[ix];
    }

    /**
     * Function call operator for convenient indexing. No bounds checking.
    */
    template<typename ...Args>
    inline T& operator()(const std::vector<size_t>& vals){
        auto ix =  _indexer.get_index(vals);
        return _data[ix];
    }

    /**
     * 
    */
    inline std::vector<T>& data(){
        return _data;
    }

    /**
     * 
    */
   template<typename ...Args>
    void print_row(size_t const & arg0, Args const & ... args){

        auto ncols = this->shape()[0];

        auto sep = "    ";
        auto prec = 3;
        auto width = 8;

        std::cout << "[";
        for (int i=0; i<ncols; i++){
            std::cout << 
            std::setw(width) <<
            std::setprecision(prec) <<
            this->at(_indexer.get_index(i, arg0, args...)) <<
            sep;
        }

        std::cout << "]\n";

    }

    /**
     * 
    */
   template<typename ...Args>
    void print_panel(size_t const & arg0, Args const & ... args){

        auto nrows = this->shape()[1];

        for (int i=0; i<nrows; i++){
            print_row(i, arg0, args...);
        }

        std::cout << "\n";

    }

    /**
     * 
    */
   template<typename ...Args>
    void print_field2d(size_t const & arg0, Args const & ... args){

        auto nrows = this->shape()[1];

        for (int i=nrows-1; i>=0; i--){
            print_row(i, arg0, args...);
        }

        std::cout << "\n";

    }

    /**
     * 
    */
   template<typename ...Args>
    void print_block(size_t const & arg0, Args const & ... args){

        auto nlayers = this->shape()[2];

        for (int i=0; i<nlayers; i++){
            print_panel(i, arg0, args...);
        }

        std::cout << "\n";
    }

    /**
     * 
    */
   template<typename ...Args>
    void print_field3d(size_t const & arg0, Args const & ... args){

        auto nlayers = this->shape()[2];

        for (int i=nlayers-1; i>=0; i--){
            print_field2d(i, arg0, args...);
        }

        std::cout << "\n";
    }

};

#endif
