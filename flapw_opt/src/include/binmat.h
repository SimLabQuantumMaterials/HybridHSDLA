#pragma once
/** \file 
 * Matrix dump reader
 */
#include <string>
#include <complex>
#include <vector>
#include <fstream>
#include <exception>
#include <utility> // for std::pair
#include <armadillo>

#include "easylogging++.h"

namespace flapw {

/** \brief Read Fleur matrix dumps
 *
 * Template class to read simple binary dumps from Fleur (by the binmat.f90 module)
 * and parse them into an Armadillo matrix/vector object.
 */
template <typename T>
class Binmat {
public:
    Binmat() = default;
    Binmat(const std::string& path);
    void read(const std::string& path);
    const std::pair<int, int> get_size() const;
    const arma::Mat<T>& get_data() const { return data; };
    static char file_test_type(const std::string& path);
protected:   
    void readheader(std::ifstream& f);
    const int bytes_header = 9;
    /** 
     * Header comprises the first 9 bytes of a file:
    *      * char: 'R' or 'C' or 'I', real or complex or integer
    *      * int: number of rows
    *      * int: number of columns
    */
    struct {
        char t;
        int  r, c;
    } header;
    arma::Mat<T> data;
    void parse_data(std::ifstream&);
    bool check_type(const char c);
};

typedef Binmat<double> RealBinmat;
typedef Binmat<std::complex<double>> ComplexBinmat;
typedef Binmat<int> IntBinmat;
/* Constructors */

/** \brief Open-and-read constructor, thin wrapper for read() */
template <class T>
Binmat<T>::Binmat(const std::string& path)
{
        read(path);
}

/* Member functions */

/** \brief Peek into the file to check the data type
 */
template <class T>
char Binmat<T>::file_test_type(const std::string& path)
{
    char c;
    std::ifstream file(path, std::ifstream::binary);
    file.read(&c, 1);
    return c;
}

/** Load data from file
*
* Opens file at `path` and tries to read it. If successful,
* use get_data() to access it.
* \throws std::invalid_argument if datatype from file does not match own type
*/
template <typename T>
void Binmat<T>::read(const std::string& path)
{
    LOG(INFO) << "Opening " << path;
    std::ifstream file;
    // make caller handle invalid file that cannot be opened
    file.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    file.open(path, std::ifstream::binary);
    readheader(file);
    if (!check_type(header.t)) {
        LOG(ERROR) << "Wrong dtype, got " << header.t;
        throw std::invalid_argument("File has wrong datatype");
    }
    parse_data(file);
}

template <typename T>
const std::pair<int, int> Binmat<T>::get_size() const
{
    return std::make_pair(header.r, header.c);
}
/**
 * Read the file header 
 *
 * See \ref header for a description of the format.
 */
template <class T>
void Binmat<T>::readheader(std::ifstream& f)
{
    f.read(&header.t, 1);
    f.read(reinterpret_cast<char*>(&header.r), 4);
    f.read(reinterpret_cast<char*>(&header.c), 4);
    LOG(DEBUG) << "Header:" << header.t << "|" 
                            << header.r << "|" 
                            << header.c; 
}


/**
 * Read the data from file into memory
 *
 * Data is stored in column-major (Fortran) order.
 */
template <typename T>
void Binmat<T>::parse_data(std::ifstream& f)
{
    //size_t num_el = header.r*header.c;
    //data.resize(num_el, 0);
    data.set_size(header.r, header.c);
    f.read(reinterpret_cast<char*>(data.memptr()), data.n_elem*sizeof(T));
}

} // namespace flapw
