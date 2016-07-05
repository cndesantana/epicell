#ifndef _CELLDIVISIONHELPERS_
#define _CELLDIVISIONHELPERS_

#include "mathematics/array.h"
#include <math.h>
#include "core/cell.h"
#include "core/vertex.h"
#include "core/edge.h"
#include "core/cellHelpers.h"
#include "core/edgeHelpers.h"
#include "core/vertexHelpers.h"
#include "printHelpers.h"
#include <functional>   // std::modulus, std::bind2nd

namespace epc {

template<int d>
class CellDivisionHelpers {
public:
	~CellDivisionHelpers(){};
	CellDivisionHelpers(){};
	int getEdgePosition(std::shared_ptr<Cell<d>> &mycell,int edgeId);
	int getOppositeEdgePosition(std::shared_ptr<Cell<d>> &mycell,int posEdgeSelected);
	void splitCell(CellCluster<d> &tissue,std::shared_ptr<Cell<d>> &cell,std::shared_ptr<Vertex<d>> &u1,std::shared_ptr<Vertex<d>> &u2,
					double defaultLambda, double defaultSdLambda, int posE1,int posE2,double A0,unsigned long iT);
	Array<double,d> randomPointInRange(double x1, double y1, double x2, double y2);
	Array<double,d> getCenterOfEdge(double x1, double y1, double x2, double y2);
	Array<double,d> getCenterOfEdge(std::shared_ptr<Edge<d>> e);
	Array<double,d> intersectionLineEdge(double x1, double y1, double c1, double c2);
	std::shared_ptr<Vertex<d>> getCentralVertexInEdge(CellCluster<d> &mytissue,std::shared_ptr<Vertex<d>> &v1,std::shared_ptr<Vertex<d>> &v2 ,double verletPar);
	std::vector<std::shared_ptr<Vertex<d>>> getListOfVertices(std::shared_ptr<Cell<d>> &mycell,int direction,std::shared_ptr<Vertex<d>> &u1,
															std::shared_ptr<Vertex<d>> &u2);
	std::vector<std::shared_ptr<Edge<d>>> updateListOfEdges(CellCluster<d> &mytissue,std::shared_ptr<Cell<d>> &mycell,int direction,
															std::shared_ptr<Vertex<d>> &u1,std::shared_ptr<Vertex<d>> &u2,
															double defaultLambda, double defaultSdLambda);
	void changeConnectivityTable(CellCluster<d> &mytissue, std::shared_ptr<Cell<d>> &oldcell,std::shared_ptr<Cell<d>> &newcell,
															std::shared_ptr<Cell<d>> &newcell2,std::shared_ptr<Vertex<d>> &u1, 
															std::shared_ptr<Vertex<d>> &u2);
	bool cellsShareEdge(CellCluster<d> &mytissue, std::shared_ptr<Cell<d>> &cell1, std::shared_ptr<Cell<d>> &cell2);
	std::shared_ptr<Edge<d>> getEdgesByVertices(std::vector<std::shared_ptr<Edge<d>>> &mydges,std::shared_ptr<Vertex<d>> &v1,std::shared_ptr<Vertex<d>> &v2);

private:
	cellHelpers<d> mycellhelper;
};

template<int d>
Array<double,d> CellDivisionHelpers<d>::randomPointInRange(double x1, double y1, double x2, double y2){
		Array<double,2> v3;
		std::random_device rd;
		std::default_random_engine generator(rd());
		std::uniform_real_distribution<double> dist(x1,x2);
		double m,b;
	
		m = (double) ((y2-y1)/(x2-x1));
		b = y1 - (double) (x1 * (double)m);
		
		double my_x = dist(generator);//Random x different from x1 and x2
		double my_y = m*my_x + b;//y associated to that random x, considering the equation of the line
		v3[0] = my_x;
		v3[1] = my_y;
		
		return(v3);
}

template<int d>
Array<double,d> CellDivisionHelpers<d>::getCenterOfEdge(double x1, double y1, double x2, double y2){
		Array<double,2> v3;
		double m,b;
	
		m = (double) ((y2-y1)/(x2-x1));
		b = y1 - (double) (x1 * (double)m);
	
		double my_x = (x2+x1)/2;//x referent to the center of the edge
		double my_y = m*my_x + b;//y associated to that random x, considering the equation of the line
		v3[0] = my_x;
		v3[1] = my_y;
		
		return(v3);
}

template<int d>
Array<double,d> CellDivisionHelpers<d>::getCenterOfEdge(std::shared_ptr<Edge<d>> e){
		std::shared_ptr<Vertex<d>> v1 = e->getVertexOne();
		std::shared_ptr<Vertex<d>> v2 = e->getVertexTwo();
		
		return 0.5*(v1->getPosition()+v2->getPosition());
}

template<int d>
Array<double,d> CellDivisionHelpers<d>::intersectionLineEdge(double x, double y, double cx, double cy){
//This is a "fake method" and needs to be implemented correctly
	Array<double,d> v4;
	double m,b;
	
	m = (double) ((cy-y)/(cx-x));
	b = y - (double) (x * (double)m);
	
	v4[0]=x;
	v4[1]=y;
		
	return(v4);
}

template<int d>
int CellDivisionHelpers<d>::getOppositeEdgePosition(std::shared_ptr<Cell<d>> &mycell,int posEdgeSelected){
	auto myedges = mycell->getEdges();
	int nEdges = myedges.size();
	int posOppositeEdge=0;


	if ((nEdges % 2) == 0){ //even number of edges
		posOppositeEdge = (posEdgeSelected + nEdges/2)%nEdges;
	}
	else{//odd number of edges
		posOppositeEdge = (posEdgeSelected + (nEdges-1)/2)%nEdges;
	}
	return(posOppositeEdge);
}

template<int d>
std::shared_ptr<Vertex<d>> CellDivisionHelpers<d>::getCentralVertexInEdge(CellCluster<d> &mytissue,std::shared_ptr<Vertex<d>> &v1,std::shared_ptr<Vertex<d>> &v2 , double verletPar){

	auto posv3 = getCenterOfEdge(v1->getPosition()[0],v1->getPosition()[1],v2->getPosition()[0],v2->getPosition()[1]);//choose the center of the selected edge to be one vertex of the new edge
	Vertex<2> v3(posv3,std::shared_ptr<Dynamics<2>>(new Verlet<2>(verletPar)));
	auto sharedv3 = std::make_shared<Vertex<2>>(v3);

	return(sharedv3);
}

template<int d>
void  CellDivisionHelpers<d>::splitCell(CellCluster<d> &mytissue,std::shared_ptr<Cell<d>> &mycell,std::shared_ptr<Vertex<d>> &u1,
									    std::shared_ptr<Vertex<d>> &u2, double defaultLambda, double defaultSdLambda, 
									    int posE1,int posE2,double A0_12,unsigned long iT){

	auto e1 = mycell->getEdges()[posE1];
	auto e2 = mycell->getEdges()[posE2];

	int posmycell = 0;
	posmycell = std::find(mytissue.getCells().begin(), mytissue.getCells().end(), mycell) - mytissue.getCells().begin();
	double K_ = mycell->getK();
	double Gamma_ = mycell->getGamma();
	double sdGamma_ = mycell->getSdGamma();

	// To create two new cells to replace "mycell", cell1 and cell2, from the list of vertices
 	// Cell 1.
 	auto listOfEdges1 = updateListOfEdges(mytissue,mycell,1,u1,u2,defaultLambda, defaultSdLambda);
	auto cell1 = std::shared_ptr<Cell<d>>(new Cell<d>(listOfEdges1, K_, A0_12, Gamma_, sdGamma_, mycell->getType(), iT));
	double signalProportion = 0.5;
	cell1->setSignal(signalProportion*mycell->getSignal());
	// Cell 2.
	auto listOfEdges2 = updateListOfEdges(mytissue,mycell,2,u1,u2,defaultLambda, defaultSdLambda);
	auto cell2 = std::shared_ptr<Cell<d>>(new Cell<d>(listOfEdges2, K_, A0_12, Gamma_, sdGamma_, mycell->getType(), iT));
	cell2->setSignal((1-signalProportion)*mycell->getSignal());

 	mytissue.addCell(cell1);
 	mytissue.addCell(cell2);
	mytissue.getVertices().push_back(u1);
	mytissue.getVertices().push_back(u2);

	posmycell = std::find(mytissue.getCells().begin(), mytissue.getCells().end(), mycell) - mytissue.getCells().begin();
	auto neighCellsOfMycell = mytissue.getConnectivity()[posmycell];
	bool e1_boundary = true;
	bool e2_boundary = true;
	for(auto neigh: neighCellsOfMycell){
		if(std::get<1>(neigh) == e1){
			auto auxcell = std::get<0>(neigh);
			auto posE1 = std::find(auxcell->getEdges().begin(), auxcell->getEdges().end(),e1) - auxcell->getEdges().begin();
			auxcell->getEdges().erase(auxcell->getEdges().begin() + posE1);
//to push back the edges: e1.One,u1 and e1.Two,u1. However, we need to use the existing edges (see the list of all edges in the tissue - tissue.getEdges())
			auto e3 = getEdgesByVertices(mytissue.getEdges(),e1->getVertexOne(),u1);	
			auto e4 = getEdgesByVertices(mytissue.getEdges(),e1->getVertexTwo(),u1);	
			auxcell->getEdges().push_back(e3);
			auxcell->getEdges().push_back(e4);
			e1_boundary = false;
		}
		else if(std::get<1>(neigh) == e2){
			auto auxcell = std::get<0>(neigh);
			auto posE2 = std::find(auxcell->getEdges().begin(), auxcell->getEdges().end(),e2) - auxcell->getEdges().begin();
			auxcell->getEdges().erase(auxcell->getEdges().begin() + posE2);
//to push back the edges: e2.One,u2 and e2.Two,u2. However, we need to use the existing edges (see the list of all edges in the tissue - tissue.getEdges())
			auto e3 = getEdgesByVertices(mytissue.getEdges(),e2->getVertexOne(),u2);
			auto e4 = getEdgesByVertices(mytissue.getEdges(),e2->getVertexTwo(),u2);
			auxcell->getEdges().push_back(e3);
			auxcell->getEdges().push_back(e4);	
			e2_boundary = false;							
		}
	}

	changeConnectivityTable(mytissue,mycell,cell1,cell2,u1,u2);

	if (mytissue.getConnectivity()[posmycell].size() > 0){
		mytissue.getConnectivity()[posmycell].erase(mytissue.getConnectivity()[posmycell].begin());
	}
	mytissue.getConnectivity().erase(mytissue.getConnectivity().begin() + posmycell);
	mytissue.removeEdge(e1);
	mytissue.removeEdge(e2);

////Removing the former edges and the former cell 	
	e1.reset();	e2.reset();
	mytissue.remCell(mycell);

	// Check if new vertices are on boundary
	if (e1_boundary){
		mytissue.addToBoundaryVertices(u1);
		std::cout << " After cell division, u1("<< u1->getId() <<") is added to boundary vertices" << std::endl;
	}
	if (e2_boundary){
		mytissue.addToBoundaryVertices(u2);
		std::cout << " After cell division, u2("<< u2->getId() <<") is added to boundary vertices" << std::endl;
	}
	return;
}

template<int d>
std::shared_ptr<Edge<d>> CellDivisionHelpers<d>::getEdgesByVertices(std::vector<std::shared_ptr<Edge<d>>> &myedges,std::shared_ptr<Vertex<d>> &v1,std::shared_ptr<Vertex<d>> &v2){
	auto lambda_ = myedges[0]->getLambda();
	auto sdLambda_ = myedges[0]->getSdLambda();
	for(auto e: myedges){
		auto ev1 = e->getVertexOne();
		auto ev2 = e->getVertexTwo();
		if( ((ev1 == v1)&(ev2==v2)) || ((ev1 == v2)&(ev2==v1)) ){
			return(e);
		}
	}

	return(std::shared_ptr<Edge<d>>(new Edge<d>(v1,v2,lambda_,sdLambda_)));
}

template<int d>
void CellDivisionHelpers<d>::changeConnectivityTable(CellCluster<d> &mytissue, std::shared_ptr<Cell<d>> &oldcell,std::shared_ptr<Cell<d>> &newcell1,std::shared_ptr<Cell<d>> &newcell2,std::shared_ptr<Vertex<d>> &u1, std::shared_ptr<Vertex<d>> &u2){
	auto allcells = mytissue.getCells();
	auto posOld = std::find(allcells.begin(), allcells.end(),oldcell) - allcells.begin();	
	auto posNew1 = std::find(allcells.begin(), allcells.end(),newcell1) - allcells.begin();
	auto posNew2 = std::find(allcells.begin(), allcells.end(),newcell2) - allcells.begin();
	auto edges1 = newcell1->getEdges();
	auto edges2 = newcell2->getEdges();	
	std::shared_ptr<Edge<2>> ce,ce2;

	mytissue.getConnectivity().resize(allcells.size());
	if(cellsShareEdge(mytissue,newcell1,newcell2)){
		ce = cellHelpers<2>::getCommonEdge(newcell1,newcell2);
		mytissue.getConnectivity()[posNew1].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(newcell2,ce));									
//		std::cout << "Adding connectivity between c" << newcell1->getId() << " and c" << newcell2->getId() << " through edge e" << ce->getId() << std::endl;
		mytissue.getConnectivity()[posNew2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(newcell1,ce));
//		std::cout << "Adding connectivity between c" << newcell2->getId() << " and c" << newcell1->getId() << " through edge e" << ce->getId() << std::endl;
	}

	auto nNeigh = mytissue.getConnectivity()[posOld].size();

	//  to remove the column associated to posOld from the rows of its Neighbour cells from the connectivity table
	for(auto iN=0;iN<(long int)nNeigh;iN++){
		auto nCell = std::get<0>(mytissue.getConnectivity()[posOld][iN]);
		auto nEdge = std::get<1>(mytissue.getConnectivity()[posOld][iN]);		
		auto posNCell = std::find(allcells.begin(), allcells.end(),nCell) - allcells.begin();
		
		if(cellsShareEdge(mytissue,nCell,newcell1)){//connecting the neighbour cell to newcell1
			ce = cellHelpers<2>::getCommonEdge(nCell,newcell1);
			mytissue.getConnectivity()[posNew1].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(nCell,ce));									
//			std::cout << "Adding connectivity between c" << newcell1->getId() << " and c" << nCell->getId() << " through edge e" << ce->getId() << std::endl;
			mytissue.getConnectivity()[posNCell].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(newcell1,ce));
//			std::cout << "Adding connectivity between c" << nCell->getId() << " and c" << newcell1->getId() << " through edge e" << ce->getId() << std::endl;
		}
		if(cellsShareEdge(mytissue,nCell,newcell2)){//connecting the neighbour cell to newcell2
			ce = cellHelpers<2>::getCommonEdge(nCell,newcell2);
			mytissue.getConnectivity()[posNew2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(nCell,ce));									
//			std::cout << "Adding connectivity between c" << newcell2->getId() << " and c" << nCell->getId() << " through edge e" << ce->getId() << std::endl;
			mytissue.getConnectivity()[posNCell].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(newcell2,ce));
//			std::cout << "Adding connectivity between c" << nCell->getId() << " and c" << newcell2->getId() << " through edge e" << ce->getId() << std::endl;
		}

		for(auto iCN=0; iCN<(long int)mytissue.getConnectivity()[posNCell].size(); iCN++){
			if(std::get<0>(mytissue.getConnectivity()[posNCell][iCN])==oldcell){
//				std::cout << "Removing Cell c" << std::get<0>(mytissue.getConnectivity()[posNCell][iCN])->getId() << " From the Neighbourhood of Cell c" << nCell->getId() << " Through Edge e" << std::get<1>(mytissue.getConnectivity()[posNCell][iCN])->getId() << std::endl;
				mytissue.getConnectivity()[posNCell].erase(mytissue.getConnectivity()[posNCell].begin() + iCN);
			}
		}
	}

	return;
}


//if cells "cell1" and "cell2" are connected by the edge "edge", this function returns the column of the table matrix. If they are not connected, the function returns -1
template<int d>
bool CellDivisionHelpers<d>::cellsShareEdge(CellCluster<d> &mytissue, std::shared_ptr<Cell<d>> &cell1, std::shared_ptr<Cell<d>> &cell2){
	auto allcells = mytissue.getCells();
	int ok=false;

	for(auto e1: cell1->getEdges()){
		for(auto e2: cell2->getEdges()){
			if(e1->getId() == e2->getId()){
				ok=true;
				return(ok);
			}
		}
	}
	return(ok);
}

template<int d>
std::vector<std::shared_ptr<Vertex<d>>> CellDivisionHelpers<d>::getListOfVertices(std::shared_ptr<Cell<d>> &mycell,int direction,std::shared_ptr<Vertex<d>> &u1,std::shared_ptr<Vertex<d>> &u2){
	std::vector<std::shared_ptr<Vertex<d>>> listOfVert;
	auto cellvertices = cellHelpers<2>::getVertices(*mycell);
	int subsetSize=0;

	cellvertices.push_back(u1);
	cellvertices.push_back(u2);	
	vertexHelpers<d>::sortVerticesCounterClockwiseInPlace(cellvertices);
	auto posv1 = std::find(cellvertices.begin(), cellvertices.end(), u1) - cellvertices.begin();

//	std::cout << "AllVerticesOfCell:";
//	for(auto iv: cellvertices){
//		std::cout << " v" << iv->getId();
//	}
//	std::cout << " --- ";
	auto nVerts = cellvertices.size();
	if ((nVerts % 2) == 0){ //even number of vertices (considering v3)
		subsetSize = nVerts/2;
	}
	else{//odd number of edges
		subsetSize = (nVerts+1)/2;
	}
//	std::cout << subsetSize+1 << " verts starting/ending by v" << u1->getId() << " that is at position " << posv1;
//	std::cout << std::endl;
//	std::cout << "AllVerticesOfCell:";
//create the list of vertices for each subset of the original cell
	if(direction == 1){//from u1 to u2:
    		for (int iv = 0; iv <= subsetSize; ++iv) {
			listOfVert.push_back(cellvertices[(posv1 + iv)%nVerts]);
	    //		std::cout << " v" << cellvertices[(posv1 + iv)%nVerts]->getId();
		}
	}
	else if(direction == 2){//from u2 to u1:
	    	for (unsigned iv = subsetSize; iv <= nVerts; ++iv) {
			listOfVert.push_back(cellvertices[(posv1 + iv)%nVerts]);
    	//	std::cout << " v" << cellvertices[(posv1 + iv)%nVerts]->getId();
		}
	}

//	std::cout << std::endl << "1st LIST OF VERTICES: ";
//	for(auto i: listOfVert){
//		std::cout << "v" << i->getId() << " ";
//	}
//	std::cout << std::endl << std::endl;

	return (listOfVert);
}

template<int d>
std::vector<std::shared_ptr<Edge<d>>> CellDivisionHelpers<d>::updateListOfEdges(
	CellCluster<d> &mytissue,std::shared_ptr<Cell<d>> &mycell,int direction,
	std::shared_ptr<Vertex<d>> &u1,std::shared_ptr<Vertex<d>> &u2,
	double defaultLambda, double defaultSdLambda)
{
//create the list of vertices for each subset of the original cell

	std::vector<std::shared_ptr<Edge<d>>> listOfEdges;
	auto myedges = mycell->getEdges();
//	unsigned long nEdges = myedges.size();
	bool found;

	auto listOfVert = getListOfVertices(mycell,direction,u1,u2);

	unsigned long nVerts = listOfVert.size();
    for (unsigned iV = 0; iV < nVerts; ++iV) {
    	found = false;
		auto v1 = listOfVert[iV];
		auto v2 = listOfVert[(iV+1)%nVerts];    	
    	for(auto ie: mytissue.getEdges()) { // we should correct this it's too long TODO
    		auto iev1 = ie->getVertexOne();
    		auto iev2 = ie->getVertexTwo();
 //			std::cout << "Comparing Edge e" << ie->getId() << " == (v" << iev1->getId() << " , v" << iev2->getId() << ")";
 //			std::cout << " With Vertices (v" << v1->getId() << " , v" << v2->getId() << ")" << std::endl; 			
    		if(
    			((iev1==v1) && (iev2==v2)) || 
    			((iev2==v1) && (iev1==v2)) ){
    				listOfEdges.push_back(ie);
//    				std::cout << "#1 - FOUND! Edge e" << ie->getId() << " == (v" << v1->getId() << " , v" << v2->getId() << ")" << std::endl;
    				found = true;
    				break;	
    		}
    	}

		if(!found){
			// we first try to find from which edge the new edge is derived
			double lambdaForNewEdge = defaultLambda;
			double sdLambdaForNewEdge = defaultSdLambda;
			for(auto ie: mytissue.getEdges()){ // we should correct this it's too long TODO
	    		auto iev1 = ie->getVertexOne();
	    		auto iev2 = ie->getVertexTwo();

	    		if (iev1 == v1 || iev2 == v2 || iev1 == v2 || iev2 == v1 ) {
	    			double edgeNorm = ie->computeLength();
	    			Array<double,2> newNormalizedVector = v1->getPosition()-v2->getPosition();
	    			newNormalizedVector /= newNormalizedVector.norm();

	    			if (arrayOperations<double,2>::isEqual(ie->computeVectorFromOneToTwo()/edgeNorm, newNormalizedVector) ||
	    			    arrayOperations<double,2>::isEqual(ie->computeVectorFromTwoToOne()/edgeNorm, newNormalizedVector) ) 
    			    {
			    		lambdaForNewEdge = ie->getLambda();
			    		sdLambdaForNewEdge = ie->getSdLambda();
			    		break;
	    			}
	    		}
    		}

			auto newEdge = std::shared_ptr<Edge<d>>(new Edge<d>(v1,v2,lambdaForNewEdge, sdLambdaForNewEdge));
			listOfEdges.push_back(newEdge);    
			myedges.push_back(newEdge);
			mytissue.getEdges().push_back(newEdge);
//			std::cout << "#2 - NOT-FOUND! Create new Edge called e" << newEdge->getId() << " == (v" << newEdge->getVertexOne()->getId() << " , v" << newEdge->getVertexTwo()->getId() << ")" << std::endl;
		}
   	}

 	return(listOfEdges);
}

template<int d>
int  CellDivisionHelpers<d>::getEdgePosition(std::shared_ptr<Cell<d>> &mycell,int edgeId){
	int pos = 0;
	for (auto e: mycell->getEdges()) {
		if (e->getId() == edgeId){
			return(pos);
		}
		else{
			pos++;
		}
	}
	return(pos);
}

}//end namespace epc

#endif
