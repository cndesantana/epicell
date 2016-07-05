#ifndef TRANSITION_H
#define TRANSITION_H

#include "transitionHelpers.h"
#include "core/epcDebug.h"
#include "core/cell.h"
#include "core/vertex.h"
#include "core/edge.h"
#include "mathematics/array.h"
#include "mathematics/rng.h"
#include <limits>

namespace epc {


///////////////////////// T1Swap
///////////////////////////////// Merging Boundary Vertices 
///////////////////////////////// Edge Rearrangement 

/////////////////////Prototypes of the classes and subclasses

template<int d>
class T1Swap{
public:
	T1Swap(double damping_){
		this->damping = damping_;
	}
	~T1Swap(){};
	bool check(CellCluster<d> &myCellCluster, double dmin= 0.1);
	void updateMBV(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> &e);//Merging Boundary Vertices
	void updateERE(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> &e, 
					std::vector<std::shared_ptr<Cell<d>>> &cellsSharingEdge , 
					std::vector<std::shared_ptr<Cell<d>>> &cellsContainingOnlyV1, 
					std::vector<std::shared_ptr<Cell<d>>> &cellsContainingOnlyV2, 
					double dmin);//Edges Rearrangement
private:
	double damping;
};

///////////////////////// Merging Boundary Vertices 

template<int d>
bool T1Swap<d>::check(CellCluster<d> &myCellCluster, double dmin) {
	transitionHelpers<d> helper;
	bool resp = false;
	for (auto e : myCellCluster.getEdges()) {
		auto elength = e->computeLength();
//		std::cout << "elength = " << elength << std::endl;
		if(elength <= dmin){
			std::vector<std::shared_ptr<Cell<d>>> cellsSharingEdge = helper.getCellsThatShareEdge(myCellCluster,e);
			std::vector<std::shared_ptr<Cell<d>>> cellsSharingV1 = helper.getCellsThatShareVertex(myCellCluster,e->getVertexOne());
			std::vector<std::shared_ptr<Cell<d>>> cellsSharingV2 = helper.getCellsThatShareVertex(myCellCluster,e->getVertexTwo());
			int nbCellsSharingEdge = cellsSharingEdge.size();	
			int nbCellsSharingV1 = cellsSharingV1.size();
			int nbCellsSharingV2 = cellsSharingV2.size();
			if((nbCellsSharingEdge == 1)&&((nbCellsSharingV1 == 1) || (nbCellsSharingV2 == 1))){
				std::cout << "T1 Merge Boundary Vertices " << e->getVertexOne()->getId() << " & " << e->getVertexTwo()->getId() << " of edge " << e->getId() << std::endl;
				this->updateMBV(myCellCluster,e);
			}else{
				std::cout << "T1 Exchange neighbours around edge " << e->getId() << " between cells ";
				for (auto c : cellsSharingEdge){
					std::cout << c->getId() << ", ";
				}
				std::cout << " of types ";
				for (auto c : cellsSharingEdge){
					std::cout << c->getType().getId() << ", ";
				}
				std::cout << std::endl;
				//std::cout << "T1 on edge ("<< e->getId()<<") with nbCellsSharingEdge = " << nbCellsSharingEdge <<", nbCellsSharingV1 = " << nbCellsSharingV1 << ", nbCellsSharingV2 = " << nbCellsSharingV2 << std::endl;
				// Get the cells containing v1 but not v2.
				std::vector<std::shared_ptr<Cell<d>>> cellsContainingOnlyV1(cellsSharingV1);
				for (int iC1 = 0; iC1<cellsSharingV1.size();iC1++){
					std::shared_ptr<Cell<d>> c1 = cellsSharingV1[iC1];
					for (int iC12 = 0; iC12<cellsSharingEdge.size(); iC12++){
						std::shared_ptr<Cell<d>> c12 = cellsSharingEdge[iC12];
						if (c12->getId() == c1->getId()){
							for (int i = 0; i<cellsContainingOnlyV1.size();i++){
								if (cellsContainingOnlyV1[i]->getId() == c1->getId()){
									cellsContainingOnlyV1.erase(cellsContainingOnlyV1.begin()+i);
									break;
								}
							}	
							break;	
						}
					}
				}
				// Get the cells containing v2 but not v1.
				std::vector<std::shared_ptr<Cell<d>>> cellsContainingOnlyV2(cellsSharingV2);
				for (int iC2 = 0; iC2<cellsSharingV2.size();iC2++){
					std::shared_ptr<Cell<d>> c2 = cellsSharingV2[iC2];
					for (int iC12 = 0; iC12<cellsSharingEdge.size(); iC12++){
						std::shared_ptr<Cell<d>> c12 = cellsSharingEdge[iC12];
						if (c12->getId() == c2->getId()){
							for (int i = 0; i<cellsContainingOnlyV2.size();i++){
								if (cellsContainingOnlyV2[i]->getId() == c2->getId()){
									cellsContainingOnlyV2.erase(cellsContainingOnlyV2.begin()+i);
									break;
								}
							}	
							break;	
						}
					}
				}
				//std::cout << "=> nbcellsContainingOnlyV2 = " << cellsContainingOnlyV1.size() <<", nbcellsContainingOnlyV2 = " << cellsContainingOnlyV2.size()<< std::endl;
				this->updateERE(myCellCluster,e,cellsSharingEdge,cellsContainingOnlyV1,cellsContainingOnlyV2,dmin);
			}
			return(true);
		}
	}
	return resp;
}

//MBV: Merging Boundary Vertices
//If two vertices at boundary edges are at a distance smaller than dmin, and at least one of them belongs to one cell only, they are merged into only one vertex
template<int d>
void T1Swap<d>::updateMBV(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> &e) {
	transitionHelpers<2> helper;
	helper.mergeEdgeAndUpdateCellCluster(myCellCluster, e, this->damping);
	return;
}


///////////////////////// EdgeRearrangement

template<int d>
void T1Swap<d>::updateERE(	CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> &e, 
							std::vector<std::shared_ptr<Cell<d>>> &cellsSharingEdge,	
							std::vector<std::shared_ptr<Cell<d>>> &cellsContainingOnlyV1, 
							std::vector<std::shared_ptr<Cell<d>>> &cellsContainingOnlyV2,
							double dmin) {
 	transitionHelpers<2> helper;
	helper.t1OnEdgeAndUpdateCellCluster(myCellCluster, e, cellsSharingEdge,cellsContainingOnlyV1,cellsContainingOnlyV2,dmin);
	return;
}


}  // end namespace epc

#endif
