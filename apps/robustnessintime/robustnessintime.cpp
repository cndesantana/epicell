//Rectangular tissue with hexagonal cells
//Tissue Relaxation
//Cell proliferation until nbDiv is reached
//Select a subset of cells (in circle or below line) and assign different type
#include "epicell.h"
#include "epicell.hh"
#include <iostream>
#include <memory>
#include <stdlib.h>
#include <time.h>
//to create a directory
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string>

#define MAXIT 10000000
using namespace epc;
using namespace std;

typedef double T;

std::string getStringByParameters(char** myargv, int myargc){ 
	std::ostringstream oss;

	oss << myargv[3];//seed
	std::string str = oss.str();
	for (int i=5;i<myargc;i++){//gamma, lambda, ncols, nrows, sdGamma, sdLambda
		str.append(std::string("_").append(std::string(myargv[i])));
	}
	return(str);
}

int main(int argc, char *argv[]){
	if(argc != 13){
		std::cout << std::endl
			<< "Incorrect number of parameters." << std::endl 
			<< "11 Parameters: modelId nRep myseed K A0 gamma lambda ncols nrows sdLambda sdGamma" << std::endl 
			<< std::endl
			<< "modelId - model-Id (1-random choice of cells; 2-eldest cells chosen before; 3- cells with higher signal first )" << std::endl
			<< "nRep - number of replicates" << std::endl
			<< "myseed - seed for random number generator" << std::endl
			<< "K - elasticity of the cell area" << std::endl
			<< "A0 -  cell preferred area" << std::endl
			<< "lambda - normalized line tension of the edges" << std::endl
			<< "gamma - normalized contractility of the cell perimeter" << std::endl
			<< "ncols - number of columns of rectangular epithelium" << std::endl
			<< "nrows - number of columns of rectangular epithelium" << std::endl
			<< "sdLambda - standard deviation of line tension of the edges" << std::endl
			<< "sdGamma - standard deviation of contractility of the cell perimeter" << std::endl
			<< "updateProb - probability to update cell mechanics" << std::endl
			<< std::endl;
		exit(1);
	}
	unsigned long modelId = atoi(argv[1]); 
	unsigned long modelId1 = modelId; 
	const unsigned long nRep = atoi(argv[2]); 
	unsigned long myseed = atoi(argv[3]);
	double K = atof(argv[4]);
	double A0 = atof(argv[5]);
	double lambdaNorm = atof(argv[6]);//give lambdaNorm as parameter
	double gammaNorm = atof(argv[7]);//give gammaNorm as parameter
	unsigned long ncols = atoi(argv[8]);
	unsigned long nrows = atoi(argv[9]);
	double sdLambda = atof(argv[10]);//give lambdaNorm as parameter
	double sdGamma = atof(argv[11]);//give gammaNorm as parameter
	double updateCellMech = atof(argv[12]);//give probability to update cell mechanics 

	unsigned long myRep = 1;

	double r0 = sqrt(2*A0/(3*sqrt(3)));//0.5;
	double dt = 0.1;
	double damping = 1.0;	
	double gamma = gammaNorm*(K*A0);//give gammaNorm as parameter
	double lambda = lambdaNorm*(K*pow(A0,1.5));//give lambdaNorm as parameter
	//double boundaryLambda = lambda/2;
	double boundaryLambda = 0.0;

	std::ostringstream oss0,oss1,oss2,oss3;

	oss0 << modelId;
	std::string oss0Str= "./tmp/model_"; 
	oss0Str.append(oss0.str());
	std::string workdir0 = oss0Str;

	struct stat st = {0};
	//Check if the upper folder exists	
	if (stat(workdir0.c_str(), &st) == -1) {
	    mkdir(workdir0.c_str(), 0700);
	}

	oss1 << gammaNorm;
	oss2 << lambdaNorm;
	std::string oss1Str = "gammaN_";
	oss1Str.append(oss1.str());
	oss1Str.append("_lambdaN_");
	oss1Str.append(oss2.str());

	std::string workdir1 = workdir0.append(std::string("/").append(oss1Str));
	
	//Check if the upper folder exists	
	if (stat(workdir1.c_str(), &st) == -1) {
	    mkdir(workdir1.c_str(), 0700);
	}

	oss3 << myseed;
	std::string oss2Str = oss3.str();
	for (int i=4;i<argc;i++){
		oss2Str.append(std::string("_").append(std::string(argv[i])));
	}
	std::string workdir2 = workdir1.append(std::string("/").append(oss2Str));

	//Check if the upper folder exists	
	if (stat(workdir2.c_str(), &st) == -1) {
	    mkdir(workdir2.c_str(), 0700);
	}

	std::cout << "nRep= "<<nRep<<std::endl; 
	while(myRep <= nRep) {
	
		unsigned long iT = 0; 
		auto dyn = std::shared_ptr<Dynamics<2>>(new Verlet<2>(damping));
	
		RNG::seed(myseed);
		NRNG::seed(myseed);

	    // Create epithelium tissue.
		Epithelium<2> epithelium(HexagonalCellsRectangularCluster<2>(nrows,ncols,r0,K,gamma,lambda,sdGamma,sdLambda,dyn, false,iT));
		double regularity = dataAnalysis::computeRegularness(epithelium.getCellCluster());
		// Create a default cell type.
	        std::shared_ptr<CellType> type0 = std::shared_ptr<CellType>(new CellType());
		std::cout << "type0 : id = " << type0->getId() << std::endl;
		epithelium.getCellCluster().getCellTypesManager().addCellType(type0);
		epithelium.getCellCluster().getCellTypesManager().setLambdaBetweenCellTypes(lambda,type0,type0);
		// Assign it to all cells.
		for (auto c : epithelium.getCellCluster().getCells()){
			c->setType(*type0);
			std::cout << "cell type : id = " << c->getType().getId()<< std::endl;
		}
		// Set the lambda of boundary edges.
		epithelium.getCellCluster().setBoundaryLambda(boundaryLambda);
		// Set boundary vertices.
		//epithelium.getCellCluster().setBoundaryVertices();
		
		double energy = dataAnalysis::computeNetEnergy(epithelium.getCellCluster());
		//Convergence (threshold, maximum-window-data)	
		double conv_relax_thresh = double(mathConstants::ForceEpsilon);
		statistics::Convergence conv_relax(conv_relax_thresh, 1000/dt);
		conv_relax.takeValue(energy);//using the energy of the tissue as a parameter to check the convergence
		double rndFloat;//random number between 0 and 1
		RNG myfrn(float(0.0),float(1.0));
		while((conv_relax.hasConverged() == false)&&(iT < MAXIT)){//the model will run until there is no significative change in the regularity of the tissue (by significative we mean that the stdev/mean of the retularity is not higher than an epsilon
			rndFloat = myfrn.fgenerate();
			if (rndFloat < updateCellMech){
                                std::cout <<  "-------------------------------" << std::endl;
                                std::cout << "BEFORE UPDATE CELL MECHANICS" << std::endl;
                                for(auto ic: epithelium.getCellCluster().getCells()){
                                   std::cout << "Gamma(c" << ic->getId() << ") = " << ic->getGamma() << std::endl;
                                   for(auto ie: ic->getEdges()){
                                       std::cout << "Lambda(e" << ie->getId() << ") = " << ie->getLambda() << std::endl;
                                   } 
                                }
				std::cout << rndFloat << " < " << updateCellMech 
                                          << std::endl << "Update cell mechanics at time step " << iT 
                                          << std::endl; 
				epithelium.updateCellMechanics(gamma,lambda,sdGamma,sdLambda);
				epithelium.getCellCluster().setBoundaryLambda(boundaryLambda);
                                std::cout << "AFTER UPDATE CELL MECHANICS" << std::endl;
                                for(auto ic: epithelium.getCellCluster().getCells()){
                                   std::cout << "Gamma(c" << ic->getId() << ") = " << ic->getGamma() << std::endl;
                                   for(auto ie: ic->getEdges()){
                                       std::cout << "Lambda(e" << ie->getId() << ") = " << ie->getLambda() << std::endl;
                                   } 
                                }

                                std::cout <<  "-------------------------------" << std::endl;
			} 

                        if (iT % 1000 == 0){
            	               std::cout << "Relaxation: iT = " << iT << ", converged = " << conv_relax.hasConverged()<< "(" << conv_relax.getNormalizedStd () << " ?< "<< conv_relax_thresh << ")" <<std::endl;
			}
			epithelium.performTimeStep(dt);
			iT++;
			regularity = dataAnalysis::computeRegularness(epithelium.getCellCluster());
			energy = dataAnalysis::computeNetEnergy(epithelium.getCellCluster());
			conv_relax.takeValue(energy);//using the energy of the tissue as a parameter to check the convergence
                        if (iT % 500 == 0){
	        		VTKwriter<2> vtk(epithelium.getCellCluster(), workdir2+"/tissue"+std::to_string(iT)+".vtk");
	        		vtk.addScalarToCell("degree",dataFunctionals2D::computeDegree);
	        		vtk.addScalarToCell("onMitosis",dataFunctionals2D::isOnMitosis);
	        		vtk.addScalarToCell("type",dataFunctionals2D::getType);
	        		vtk.addScalarToVertex("area",dataFunctionals2D::computeArea);
	        		vtk.addVectorToVertex("force",dataFunctionals2D::computeForce);
				cellClusterIO<2>::writeREG(workdir2+"/regularness_celldiv.reg",iT,regularity,gammaNorm,lambdaNorm,sdGamma,sdLambda);
				cellClusterIO<2>::writeSTA(epithelium.getCellCluster(),workdir2+"/tissue_celldiv.sta",iT,gammaNorm,lambdaNorm,sdGamma,sdLambda);
			}//end-if
		}//end-while of Relaxation
        
        	std::cout << "Relaxation: iT = " << iT << ", converged = " << conv_relax.hasConverged()<< "(" << conv_relax.getNormalizedStd () << " ?< "<< conv_relax_thresh << ")" <<std::endl;
				
		//Start cell divisions.
		std::vector<double> vecOfProportions = dataAnalysis::computeDistributionOfSides(epithelium.getCellCluster());
		vecOfProportions = dataAnalysis::computeDistributionOfSides(epithelium.getCellCluster());
		int nbDivisions = 0;
		int nbDivMax = 300;
		T1Swap<2> t1(damping);
		while((nbDivisions < nbDivMax) && (iT < 2*MAXIT)){//the model will run until there is no significant change in the proportion of cells with 3,4,5,6,7,8,9 or more sides (by significant we mean that the stdev/mean of the regularity is not higher than an epsilon)
			rndFloat = myfrn.fgenerate();
			if (rndFloat < updateCellMech){
                                std::cout <<  "-------------------------------" << std::endl;
                                std::cout << "BEFORE UPDATE CELL MECHANICS" << std::endl;
                                for(auto ic: epithelium.getCellCluster().getCells()){
                                   std::cout << "Gamma(c" << ic->getId() << ") = " << ic->getGamma() << std::endl;
                                   for(auto ie: ic->getEdges()){
                                       std::cout << "Lambda(e" << ie->getId() << ") = " << ie->getLambda() << std::endl;
                                   } 
                                }
				std::cout << rndFloat << " < " << updateCellMech 
                                          << std::endl << "Update cell mechanics at time step " << iT 
                                          << std::endl; 
				epithelium.updateCellMechanics(gamma,lambda,sdGamma,sdLambda);
				epithelium.getCellCluster().setBoundaryLambda(boundaryLambda);
                                std::cout << "AFTER UPDATE CELL MECHANICS" << std::endl;
                                for(auto ic: epithelium.getCellCluster().getCells()){
                                   std::cout << "Gamma(c" << ic->getId() << ") = " << ic->getGamma() << std::endl;
                                   for(auto ie: ic->getEdges()){
                                       std::cout << "Lambda(e" << ie->getId() << ") = " << ie->getLambda() << std::endl;
                                   } 
                                }

                                std::cout <<  "-------------------------------" << std::endl;
			} 

			epithelium.performTimeStep(dt);

			IncrementedCellGrowth<2> cellgr;//strategy of cell growth
			
			ShortestAxisCellDivision<2> celldiv;//strategy of cell division (along shortex cell axis)
//			int modelId1 = 2;
			std::shared_ptr<Cell<2>> rndCell;
			if((iT % 100)==0){//The growth and division of cells occur only each 100 time steps 
				// Select the cell that will enter mitosis (start growing).
				if(modelId1 == 1){//random choice of cells
					std::shared_ptr<CellSelection<2>> RandomCellSelect(new RandomCellSelection<2>());
					rndCell = RandomCellSelect->select(epithelium.getCellCluster(),iT);
				}
				else if(modelId1 == 2){//eldest cells chosen most probabaly first
					std::shared_ptr<CellSelection<2>> EldestCellSelect(new EldestCellSelection<2>());
					rndCell = EldestCellSelect->select(epithelium.getCellCluster(),iT);						
				}
				if(rndCell){		
					epithelium.getCellCluster().addCellOnMitosis(rndCell);	
				}

				//Once the cells are added to the cellOnMitosis list
				auto cellsOnMitosis = epithelium.getCellCluster().getCellsOnMitosis();
				if(cellsOnMitosis.size() > 0){//If there is any cell in CellsOnMitosis list
					//std::cout << " -------- Cells on Mitosis -------- " << std::endl;
					for(auto incCell:cellsOnMitosis){
						//std::cout << incCell->getId() << std::endl;
						cellgr.grow(incCell,0.1*incCell->getAreaBeforeMitosis());//increase the preferred area of the chosen cell by a factor of 20%
					}	
				}

				//After increasing the area of the cells on Mitosis, we will divide the cells that exceed the threshold area (2 * A0)
				auto cellsToDivide = epithelium.getCellCluster().getCellsToDivide();//get the cells wich current area is the double of their initial area
				//std::cout << " -------- Divided Cells -------- " << std::endl;
				if (cellsToDivide.size() > 0){
					for(auto cellToDivide:cellsToDivide){
						//std::cout << cellToDivide->getId() << std::endl;
						celldiv.divide(epithelium.getCellCluster(),cellToDivide,lambda, sdLambda, A0,damping,iT);
						nbDivisions+=1;
					}
					if (nbDivisions %10 == 0){
						std::cout << "nbDivisions = " << nbDivisions << std::endl;
					}
				}
				if ((iT % 2000) == 0){
					bool done_t1 = true;
					while(done_t1){
						done_t1 = t1.check(epithelium.getCellCluster(),0.1*r0);
					}	
				}
			}//end - if ((iT % 1000) ==0)
			iT++;
			regularity = dataAnalysis::computeRegularness(epithelium.getCellCluster());
			vecOfProportions = dataAnalysis::computeDistributionOfSides(epithelium.getCellCluster());
			if (iT % 1000 == 0){
//	        		VTKwriter<2> vtk(epithelium.getCellCluster(), workdir2+"/tissue"+std::to_string(iT)+".vtk");
//	        		vtk.addScalarToCell("degree",dataFunctionals2D::computeDegree);
//	        		vtk.addScalarToCell("onMitosis",dataFunctionals2D::isOnMitosis);
//	        		vtk.addScalarToCell("type",dataFunctionals2D::getType);
// 	        		vtk.addScalarToCell("prolifSignal",dataFunctionals2D::getProlifSignal);
// 	        		vtk.addScalarToCell("sqrProlifSignal",dataFunctionals2D::getSqrProlifSignal);
// 	        		vtk.addScalarToCell("cubicProlifSignal",dataFunctionals2D::getCubicProlifSignal);
//	        		vtk.addScalarToVertex("area",dataFunctionals2D::computeArea);
//	        		vtk.addVectorToVertex("force",dataFunctionals2D::computeForce);
				
				cellClusterIO<2>::writeREG(workdir2+"/regularness_celldiv.reg",iT,regularity,gammaNorm,lambdaNorm,sdGamma,sdLambda);
				cellClusterIO<2>::writeSTA(epithelium.getCellCluster(),workdir2+"/tissue_celldiv.sta",iT,gammaNorm,lambdaNorm,sdGamma,sdLambda);
				cellClusterIO<2>::writeDIST(workdir2+"/distributionofshapes.dist",iT,vecOfProportions);
			}
		}//end-while modelId=2
//		regularity = dataAnalysis::computeRegularness(epithelium.getCellCluster());
//		vecOfProportions = dataAnalysis::computeDistributionOfSides(epithelium.getCellCluster());
//		cellClusterIO<2>::writeREG(workdir2+"/regularness_celldiv.reg",iT,regularity,gammaNorm,lambdaNorm,sdGamma,sdLambda);
//		cellClusterIO<2>::writeSTA(epithelium.getCellCluster(),workdir2+"/tissue_celldiv.sta",iT,gammaNorm,lambdaNorm,sdGamma,sdLambda);
//		cellClusterIO<2>::writeDIST(workdir2+"/distributionofshapes.dist",iT,vecOfProportions);

		myseed++;
		myRep= myRep+1;
	}//end-while-myRep
	return 0;
}
