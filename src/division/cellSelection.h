#ifndef CELL_SELECTION_H
#define CELL_SELECTION_H

#include "mathematics/rng.h"
#include "core/epcDebug.h"
#include "core/cell.h"
#include <vector> 

namespace epc {

template<int d>
class 	CellSelection {
public:
	virtual ~CellSelection(){}
	virtual std::shared_ptr<Cell<d>> select(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) = 0;
	std::shared_ptr<Cell<d> > select(const CellCluster<d> &tissue,unsigned long iT) {
		return select(tissue.getCells(),iT);
	}
	virtual std::vector<std::shared_ptr<Cell<d>>> selectCells(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) = 0;
	std::vector<std::shared_ptr<Cell<d>>> selectCells(const CellCluster<d> &tissue,unsigned long iT) {
		return selectCells(tissue.getCells(),iT);
	}
};

// ======================================================================== //
// ========================RandomCellSelection============================ //
// ======================================================================== //

template<int d>
class RandomCellSelection final : public CellSelection<d> {
//class RandomCellSelection final {
public:
	virtual std::shared_ptr<Cell<d> > select(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) override;
	virtual std::vector<std::shared_ptr<Cell<d>>> selectCells(const std::vector<std::shared_ptr<Cell<d>>> &cells, unsigned long iT) override;
};

template<int d>
std::shared_ptr<Cell<d> > RandomCellSelection<d>::select(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) {
	// randomly choose a cell and tell it to start the division...
	
	int ncells=cells.size();
    RNG myirn(0,(ncells-1));
	int cellSelected = myirn.igenerate();

	// now to increase the size of "cellSelected".	
	return(cells[cellSelected]);
}

template<int d>
std::vector<std::shared_ptr<Cell<d>>> RandomCellSelection<d>::selectCells(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) {
	std::vector<std::shared_ptr<Cell<d>>> selectdCells;

	/// TO DO

	return selectdCells;
}


// ======================================================================== //
// ========================EldestCellSelection============================ //
// ======================================================================== //

template<int d>
class EldestCellSelection final : public CellSelection<d> {
public:
	virtual std::shared_ptr<Cell<d>> select(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) override;
	virtual std::vector<std::shared_ptr<Cell<d>>> selectCells(const std::vector<std::shared_ptr<Cell<d>>> &cells, unsigned long iT) override;
};

template<int d>
std::shared_ptr<Cell<d>> EldestCellSelection<d>::select(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) {
	CellSelectionHelpers<2> helper;

	auto ageAllCells = helper.getAgeAllCells(cells,iT);
	auto ageCum = helper.CumulativeSum(ageAllCells);

        RNG myfrn(float(0.0),float(1.0));
	auto myrand = myfrn.fgenerate();//generate a random integer number

	auto cellSelected = helper.Multinomial(ageCum,myrand);

	// chose a cell according to its age (eldest cell first) and tell it to start the division...
	return(cells[cellSelected]);
}

template<int d>
std::vector<std::shared_ptr<Cell<d>>> EldestCellSelection<d>::selectCells(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) {
	std::vector<std::shared_ptr<Cell<d>>> selectdCells;

	/// TO DO

	return selectdCells;
}

// ======================================================================== //
// ========================HigherSignalCellSelection============================ //
// ======================================================================== //

template<int d>
class HigherSignalCellSelection final : public CellSelection<d> {
public:
	HigherSignalCellSelection(){};
	HigherSignalCellSelection(double signalThresh_);
	virtual std::shared_ptr<Cell<d>> select(const std::vector<std::shared_ptr<Cell<d>>> &cells, unsigned long iT) override;
	virtual std::vector<std::shared_ptr<Cell<d>>> selectCells(const std::vector<std::shared_ptr<Cell<d>>> &cells, unsigned long iT) override;
private:
	double signalThresh;
};

template<int d>
HigherSignalCellSelection<d>::HigherSignalCellSelection(double signalThresh_):signalThresh(signalThresh_){}

template<int d>
std::shared_ptr<Cell<d>> HigherSignalCellSelection<d>::select(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) {
	CellSelectionHelpers<2> helper;

	//std::vector<double> signalAllCells = helper.getSignalAllCells(cells);
	//std::vector<double> signalAllCells = helper.getSqrSignalAllCells(cells);
	std::vector<double> signalAllCells = helper.getCubicSignalAllCells(cells);
	std::vector<float> signalCum = helper.CumulativeSum(signalAllCells);

	if (signalCum[signalCum.size()-1] == 0){
		std::cout << "WARNING! Select cell based on Signal but all cells have signal = 0 ..." << std::endl;
		return std::shared_ptr<Cell<d>>(nullptr);
	}

    RNG myfrn(float(0.0),float(1.0));
	auto myrand = myfrn.fgenerate();//generate a random integer number

	auto cellSelected = helper.Multinomial(signalCum,myrand);

	// chose a cell according to its signal (eldest cell first) and tell it to start the division...
	return(cells[cellSelected]);
}

template<int d>
std::vector<std::shared_ptr<Cell<d>>> HigherSignalCellSelection<d>::selectCells(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT) {
	std::vector<std::shared_ptr<Cell<d>>> selectdCells;

	for (auto c: cells){
		// if (c->getSignal()>=signalThresh){
		// 	selectdCells.push_back(c);
		// }
		RNG myfrn(float(0.0),float(1.0));
		auto myrand = myfrn.fgenerate();//generate a random integer number
		double signalSigmoid = 1/(1+exp(-(c->getSignal()-signalThresh)));
		if (myrand < signalSigmoid){
			selectdCells.push_back(c);
		}

	}

	if (selectdCells.size() == 0){
		std::cout << "WARNING! Select cells based on Signal > " << signalThresh << " but no cell was selected ..." << std::endl;
	}

	return selectdCells;
}

}  // end namespace epc

#endif
