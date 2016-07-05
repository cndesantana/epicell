#ifndef CELL_GROWTH_H
#define CELL_GROWTH_H

#include "core/epcDebug.h"
#include "core/cell.h"


namespace epc {

template<int d>
class CellGrowth {
public:
	virtual ~CellGrowth(){};

	virtual void grow(std::shared_ptr<Cell<d>> &cell, double incr) = 0;
};


template<int d>
class IncrementedCellGrowth final : public CellGrowth<d> {
public:
	virtual void grow(std::shared_ptr<Cell<d>> &cell, double incr) override;
};


template<int d>
void IncrementedCellGrowth<d>::grow(std::shared_ptr<Cell<d>> &cell, double incr) {
	// implement the growth type of the cell.
	double actualA0 = cell->getA0(); 
//	cell->setA0(incr*actualA0);//increases the preferred area (A0) of the cell
	cell->setA0(incr+actualA0);//increases the preferred area (A0) of the cell
	return;
}

}  // end namespace epc

#endif
