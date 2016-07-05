#ifndef CELLTYPE_H
#define CELLTYPE_H

namespace epc {

class CellType {
public:
    CellType();
private:
	/// increments the id of a cell type
    void createNewId();
public:
	unsigned long const& getId() const;
	
private:
	unsigned long id; // unique identifier for each CellType

	static unsigned long cellTypesIds;
};

unsigned long CellType::cellTypesIds = 0;
}  // end namespace epc

#endif
