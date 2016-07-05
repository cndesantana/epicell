#ifndef CELL_IO_H
#define CELL_IO_H

#include <fstream>
#include <string>
#include <vector>
#include "io.h"
#include "core/cell.h"
#include "core/epcDebug.h"

namespace epc {

template<int d>
struct cellIO {

/// writes one cell in file fname
static void writeCell(const Cell<d> &cell, std::string fname) {
    std::ofstream fout(fname.c_str());
    fout << cell;
    fout.close();
}

/// writes a list of cells in a file named fname
static void writeCells(const std::vector<std::shared_ptr<Cell<d>>> &cells, std::string fname) {
    std::ofstream fout(fname.c_str());
    for (auto c:cells) fout << *c << std::endl;
    fout.close();
}

};

} // end namespace epc

#endif

