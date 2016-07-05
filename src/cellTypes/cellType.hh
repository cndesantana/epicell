#ifndef CELLTYPE_HH
#define CELLTYPE_HH

#include "cellType.h"

namespace epc {

// =============== Constuctors ========================= //
CellType::CellType()
{
	createNewId();
}



void CellType::createNewId() {
    //static unsigned long cellTypesIds = 0; 
    std::cout << " --> start CellType::createNewId" << std::endl;
    id = cellTypesIds++;
    std::cout << " --> end CellType::createNewId" << std::endl;
}

// ========== Interface functions ================ //
unsigned long const& CellType::getId() const {
    return id;
}

}  // end namespace epc

#endif
