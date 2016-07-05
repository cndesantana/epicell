#ifndef CELLTYPESMANAGER_HH
#define CELLTYPESMANAGER_HH


#include "cellTypesManager.h"

namespace epc {

// =============== Constuctors ========================= //

CellTypesManager::CellTypesManager()
{
	cellTypes.clear();
}

CellTypesManager::~CellTypesManager()
{
	/*cellTypes.clear();
	for (unsigned i = 0; i< lambdas.size(); i++){
		lambdas[i].clear();
	}
	lambdas.clear();*/
}

// Add a cell type.
void CellTypesManager::addCellType(std::shared_ptr<CellType> &type) {
	std::cout << "--> start addCellType()" << std::endl;
	cellTypes.push_back(type);
	lambdas.resize(cellTypes.size());
	std::cout << "<-- end addCellType()" << std::endl;
}

// Set the line tension along the edges separating a cell of type1 and a cell of type2.
void CellTypesManager::setLambdaBetweenCellTypes(double lambda_, std::shared_ptr<CellType> &type1, std::shared_ptr<CellType> &type2){
	// std::cout << "--> start setLambdaBetweenCellTypes("<< type1.getId() << ", " << type2.getId() <<" )" << std::endl;
	// Get the indices of the given types in the list of cellTypes.
	int i1 = -1;
	int i2 = -1;
	for (int iA = 0; iA < cellTypes.size(); iA ++){
		if (cellTypes[iA]->getId() == type1->getId()){
			i1 = iA;
		}
		if (cellTypes[iA]->getId() == type2->getId()){
			i2 = iA;
		}
		if ((i1 > -1) && (i2 >-1)){
			break;
		}
	}
	EPC_ASSERT(((i1>-1) && (i2>-1)) &&
                        "Error Setting Lambda between inexisting cell types");

	// std::cout << "index of type("<< type1.getId() << ") = " << i1 <<" & index of type("<< type2.getId() << ") = " << i2 << std::endl; 
	// Set lambda between type 1 and type 2.
	bool found12 = false;
	for (int iB = 0; iB < lambdas[i1].size();iB++){
		if (std::get<0>(lambdas[i1][iB])->getId()==type2->getId()){
			std::get<1>(lambdas[i1][iB]) = lambda_;
			found12 = true; 
			break;
		}
	}

	if (found12 == false){
		// std::cout << "Lambda between these types did not exist" << std::endl;
		lambdas[i1].push_back(std::tuple<std::shared_ptr<CellType>,double>(type2,lambda_));		
		// std::cout << " Set Lambda between type1->type2 done" << std::endl;
		// for (int iB=0; iB < lambdas[i1].size(); iB++){
		// 	std::cout << "lambda( " << cellTypes[i1]->getId() << ", " << std::get<0>(lambdas[i1][iB])->getId() << ") = " << std::get<1>(lambdas[i1][iB])<< std::endl;
		
		// }
	}

	// Set lambda between type 2 and type 1.
	bool found21 = false;
	for (int iB = 0; iB < lambdas[i2].size();iB++){
		if (std::get<0>(lambdas[i2][iB])->getId()==type1->getId()){
			std::get<1>(lambdas[i2][iB]) = lambda_;
			found21 = true; 
			break;
		}
	}
	if (found21 == false){
		// std::cout << "Lambda between these types did not exist" << std::endl;
		lambdas[i2].push_back(std::tuple<std::shared_ptr<CellType>,double>(type1,lambda_));
		// std::cout << " Set Lambda between type2->type1 did not exist" << std::endl;
	}
	// std::cout << "<-- end setLambdaBetweenCellTypes()" << std::endl;
}

// Get the line tension along the edges separating a cell of type1 and a cell of type2.
double CellTypesManager::getLambdaBetweenCellTypes(CellType type1, CellType type2){
	// std::cout << "--> start getLambdaBetweenCellTypes("<< type1.getId() << ", " << type2.getId() <<" )" << std::endl;
	// Get the indices of the given types in the list of cellTypes.
	int i1 = -1;
	int i2 = -1;
	for (int iA = 0; iA < cellTypes.size(); iA ++){
		if (cellTypes[iA]->getId() == type1.getId()){
			i1 = iA;
		}
		if (cellTypes[iA]->getId() == type2.getId()){
			i2 = iA;
		}
		if ((i1 > -1) && (i2 >-1)){
			break;
		}
	}
	EPC_ASSERT(((i1>-1) && (i2>-1)) &&
                        "Error Getting Lambda between inexisting cell types");
	// std::cout << "index of type("<< type1.getId() << ") = " << i1 <<" & index of type("<< type2.getId() << ") = " << i2 << std::endl; 
	
	double lambda12;
	bool found12 = false;
	// std::cout << "lambdas[i1].size() = " << lambdas[i1].size() << std::endl;
	for (int iB = 0; iB < lambdas[i1].size(); iB++){
		// std::cout << "lambda("<< cellTypes[i1]->getId()<< ", "<< std::get<0>(lambdas[i1][iB])->getId() << ") = " << std::get<1>(lambdas[i1][iB]) << std::endl;
		if (std::get<0>(lambdas[i1][iB])->getId()==type2.getId()){			
			lambda12 = std::get<1>(lambdas[i1][iB]);
			found12 = true; 
			break;
		}
	}
	EPC_ASSERT((found12 == true) &&
                        "Error Lambda between cell types 1->2 not found");
	double lambda21;
	bool found21 = false;
	for (int iB = 0; iB < lambdas[i2].size();iB++){
		if (std::get<0>(lambdas[i2][iB])->getId()==type1.getId()){
			lambda21 = std::get<1>(lambdas[i2][iB]);
			found21 = true; 
			break;
		}
	}
	EPC_ASSERT((found21 == true) &&
                        "Error Lambda between cell types 2->1 not found");

	EPC_ASSERT((lambda12 == lambda21) &&
                        "Error lambda12 different from lambda21");

	// std::cout << "<-- end getLambdaBetweenCellTypes()" << std::endl;

	return lambda21;

}

}  // end namespace epc

#endif
