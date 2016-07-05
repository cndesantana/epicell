#ifndef CELLTYPESMANAGER_H
#define CELLTYPESMANAGER_H

#include "cellType.h"

namespace epc {


class CellTypesManager {
public:
    CellTypesManager();
    ~CellTypesManager();

    void addCellType(std::shared_ptr<CellType> &type);
    void setLambdaBetweenCellTypes(double lambda_, std::shared_ptr<CellType> &type1, std::shared_ptr<CellType> &type2);
    double getLambdaBetweenCellTypes(CellType type1, CellType type2);


private:
	std::vector<std::shared_ptr<CellType>> cellTypes;
	std::vector<std::vector<std::tuple<std::shared_ptr<CellType>, double>>> lambdas;
};


}  // end namespace epc

#endif
