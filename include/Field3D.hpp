#ifndef FIELD_3D_HPP
#define FIELD_3D_HPP

#include <array>
#include <iomanip>
#include <string>
#include <vector>

template<typename T>
class Field3{

    private:
        std::vector<T>  _data;
        std::array<size_t,3> _strides;
        size_t          _size;
 
    public:

    /**
     * Default constructor. Create a field with a single element.
    */
    Field3(){
        _strides = {1,1,1};
        _size = 1;
        _data.resize(_size);
        this->zeros();
    }

    /**
     * Preferred constructor. Create a field with dimensions (i,j,k)
     * and initialize all elements to zero.
    */
    Field3(size_t i, size_t j, size_t k){
        _strides = {i,j,k};
        _size = i*j*k;
        _data.resize(_size);
        this->zeros();
    }

    /**
     * Alternate constructor from a std::vector of size parameters.
    */
    Field3(const std::vector<size_t>& sizes){
        _strides = {sizes[0], sizes[1], sizes[2]};
        _size = sizes[0]*sizes[1]*sizes[2];
        _data.resize(_size);
        this->zeros();  
    }

    /**
     * Alternate constructor from a std::array of size parameters.
    */
    Field3(const std::array<size_t,3>& sizes){
        _strides = {sizes[0], sizes[1], sizes[2]};
        _size = sizes[0]*sizes[1]*sizes[2];
        _data.resize(_size);
        this->zeros();  
    }

    /**
     * Reshape a field in 3 dimensions. This operation requires that the 
     * length of the master data array remains unchanged.
    */
    inline void reshape(size_t i, size_t j, size_t k){
        if(i*j*k != _data.size()){
            std::string err_str = "Cannot reshape field of shape (";
            err_str += std::to_string(_strides[0]);
            err_str += ", ";
            err_str += std::to_string(_strides[1]);
            err_str += ", ";
            err_str += std::to_string(_strides[2]);

            err_str += ") to shape (";
            err_str += std::to_string(i);
            err_str += ", ";
            err_str += std::to_string(j);
            err_str += ", ";
            err_str += std::to_string(k);
            err_str += ").\n";

            throw std::length_error(err_str);
        }

        _strides[0] = i;
        _strides[1] = j;
        _strides[2] = k;
    }

    /**
     * Fill an array with zero-valued elements.
    */
    inline void zeros(){std::fill(_data.begin(), _data.end(), 0);}

    /**
     * Fill an array with one-valued elements.
    */
    inline void ones(){std::fill(_data.begin(), _data.end(), 1);}

    /**
     * Fill an array with elements of a specified value.
    */
    inline void fill(const T& val){
        std::fill(_data.begin(), _data.end(), val);
    }

    /**
     * Convert a 3D (i,j,k) index to a 1D index into the master data array.
    */
    inline size_t get_index(size_t i, size_t j, size_t k){
        return i + j*_strides[0] + k*_strides[1]*_strides[0];
    }

    /**
     * Convert a 3D (i,j,k) index in STL vector format to a 1D index into the 
     * master data array.
    */
    inline size_t get_index(const std::vector<size_t>& vals){
        return vals[0] + vals[1]*_strides[0] + vals[2]*_strides[1]*_strides[0];
    }

    /**
     * Convert a 3D (i,j,k) index in STL array format to a 1D index into the 
     * master data array.
    */
    inline size_t get_index(const std::array<size_t,3>& vals){
        return vals[0] + vals[1]*_strides[0] + vals[2]*_strides[1]*_strides[0];
    }

    /**
     * Convert a 3D (i,j,k) index in C-array format to a 1D index into the 
     * master data array.
    */
    inline size_t get_index(const size_t* vals){
        return vals[0] + vals[1]*_strides[0] + vals[2]*_strides[1]*_strides[0];
    }

    /**
     * Get the element at (i,j,k) with bounds checking. Indices are specified
     * in STL vector format. 
    */
    inline T& at(const std::vector<size_t>& vals){
        auto idx = get_index(vals);
        return _data.at(idx);
    }

    /**
     * Get the element at (i,j,k) with bounds checking. Indices are specified
     * in STL array format. 
    */
    inline T& at(const std::array<size_t,3>& vals){
        auto idx = get_index(vals);
        return _data.at(idx);
    }

    /**
     * Get the element at (i,j,k) with bounds checking. Indices are specified
     * as individual arguments. 
    */
    inline T& at(size_t i, size_t j, size_t k){
        auto idx = get_index(i,j,k);
        return _data.at(idx);
    }

    /**
     * Get the element at index i in the raw array with bounds checking.
    */
    inline T& at(size_t i){
        return _data.at(i);
    }

    /**
     * Get the element at (i,j,k). Indices are specified as individual
     * arguments. 
    */
    inline T& operator()(size_t i, size_t j, size_t k){
        auto idx = get_index(i,j,k);
        return _data[idx];
    }

    /**
     * Get the element at (i,j,k) with bounds checking. Indices are specified
     * in STL vector format. 
    */
    inline T& operator()(const std::vector<size_t>& vals){
        auto idx = get_index(vals);
        return _data[idx];
    }

    /**
     * Get the element at (i,j,k) with bounds checking. Indices are specified
     * in STL array format. 
    */
    inline T& operator()(const std::array<size_t,3>& vals){
        auto idx = get_index(vals);
        return _data[idx];
    }

    /**
     * Get the element at (i,j,k) with bounds checking. Indices are specified
     * in C-array format. 
    */
    inline T& operator()(const size_t* vals){
        auto idx = get_index(vals);
        return _data[idx];
    }

    /**
     * Print a single row of the field to stdout. The row is specified at
     * (i, j, k).  
    */
    void print_row(size_t i=0, size_t j=0, size_t k=0){

        auto ncols = this->shape()[0];

        auto sep = "    ";
        auto prec = 10;
        auto width = 12;

        std::cout << "[";
        for (size_t i=0; i<ncols; i++){
            std::cout << 
            std::setw(width) <<
            std::setprecision(prec) <<
            this->at(get_index(i, j, k)) <<
            sep;
        }

        std::cout << "]\n";
    }

    /**
     * Print a single panel of the field to stdout. The panel is specified at
     * (i, j, k).  
    */
    void print_panel(size_t i=0, size_t j=0, size_t k=0){

        auto nrows = this->shape()[1];

        for (size_t i=0; i<nrows; i++){
            print_row(i, j, k);
        }
        std::cout << "\n";
    }

    /**
     * Print a single panel of the field to stdout. The panel is specified at
     * (i, j, k) and is printed with the maximum y-index oriented at the TOP of
     * the screen.  
    */
    void print_field2d(size_t i=0, size_t j=0, size_t k=0){

        auto nrows = this->shape()[1];

        for (size_t i=nrows-1; i>=0; i--){
            print_row(i, j, k);
        }
        std::cout << "\n";
    }

    /**
     * Print a block of the field to stdout. The panel is specified at
     * (i, j, k) and is printed with the maximum y-index oriented at the
     * BOTTOM of the screen.  
    */
    void print_block(size_t i=0, size_t j=0, size_t k=0){

        auto nlayers = this->shape()[2];

        for (size_t i=0; i<nlayers; i++){
            print_panel(i, j, k);
        }
        std::cout << "\n";
    }

    /**
     * Print a block of the field to stdout. The panel is specified at
     * (i, j, k) and is printed with the maximum y-index oriented at the
     * TOP of the screen.  
    */
    void print_field3d(size_t i=0, size_t j=0, size_t k=0){

        auto nlayers = this->shape()[2];

        for (size_t i=nlayers-1; i>=0; i--){
            print_field2d(i, j, k);
        }
        std::cout << "\n";
    }

    /**
     * Return the shape of the field in 3-dimension.
    */
    std::array<size_t, 3> shape(){
        return _strides;
    }

    /**
     * Return the total number of elements in the field.
    */
    size_t size(){
        return _data.size();
    }

    /**
     * Return a reference to the master data array of the field.
    */
    std::vector<T>& data(){
        return _data;
    }
};

/**
 * Throw an error indicating a shape mismatch between two fields.
*/
template<typename T>
void shape_mismatch(const Field3<T>& a, const Field3<T>& b){

    std::string asize_str = "(";
    std::string bsize_str = "(";

    for(int i=0; i<3; i++){
        asize_str += std::to_string(a.shape()[i]);
        bsize_str += std::to_string(a.shape()[i]);
        if(i<2){
            asize_str += ", ";
            bsize_str += ", ";
        }
        asize_str += ")";
        bsize_str += ")";
    }

    auto err_str = "Shape Mismatch: field A shape " + 
                    asize_str +
                    "does not match field B shape " +
                    bsize_str +
                    "\n";

    throw std::length_error(err_str);
}

/**
 * Perform an element-wise sum of two 3D fields.
*/
template<typename T>
Field3<T> operator+(const Field3<T>& a, const Field3<T>& b){

    // Shapes of fields need to match for this operation to make sense
    if(a.shape() != b.shape()){shape_mismatch(a, b);}

    // Create a new field with the same shape as the operands
    auto newshape = a.shape();
    Field3<T> result(newshape[0], newshape[1], newshape[2]);

    for(int i=0; i<result.size(); i++){
        result.data()[i] = a.data()[i] + b.data()[i];
    }
}

/**
 * Perform an element-wise subtraction of two 3D fields.
*/
template<typename T>
Field3<T> operator-(const Field3<T>& a, const Field3<T>& b){

    // Shapes of fields need to match for this operation to make sense
    if(a.shape() != b.shape()){shape_mismatch(a, b);}

    // Create a new field with the same shape as the operands
    auto newshape = a.shape();
    Field3<T> result(newshape[0], newshape[1], newshape[2]);

    for(int i=0; i<result.size(); i++){
        result.data()[i] = a.data()[i] - b.data()[i];
    }
}

/**
 * Perform an element-wise multiplcation of two 3D fields.
*/
template<typename T>
Field3<T> operator*(const Field3<T>& a, const Field3<T>& b){

    // Shapes of fields need to match for this operation to make sense
    if(a.shape() != b.shape()){shape_mismatch(a, b);}

    // Create a new field with the same shape as the operands
    auto newshape = a.shape();
    Field3<T> result(newshape[0], newshape[1], newshape[2]);

    for(int i=0; i<result.size(); i++){
        result.data()[i] = a.data()[i] * b.data()[i];
    }
}

/**
 * Perform an element-wise division of two 3D fields.
*/
template<typename T>
Field3<T> operator/(const Field3<T>& a, const Field3<T>& b){

    // Shapes of fields need to match for this operation to make sense
    if(a.shape() != b.shape()){shape_mismatch(a, b);}

    // Create a new field with the same shape as the operands
    auto newshape = a.shape();
    Field3<T> result(newshape[0], newshape[1], newshape[2]);

    for(int i=0; i<result.size(); i++){
        result.data()[i] = a.data()[i] * b.data()[i];
    }
}

#endif