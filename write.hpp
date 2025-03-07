#ifndef write_hpp
#define write_hpp

#include "initialize.hpp"
#include "initialize.cpp"

namespace FileUtils
{
    int demintions = SIZE + 1;
    void write(std::ofstream &output, double tstar, int spins[SIZE + 1]);
    void w_output(std::ofstream &file, double tstar, double *m);
    void w_error(std::ofstream &outputBining, double tstar, double *bins, double *m);
}
#endif // WRITE_H
