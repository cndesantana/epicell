#ifndef _TRANSITIONHELPERS_H
#define _TRANSITIONHELPERS_H

#include "mathematics/array.h"
#include <math.h>
#include "core/cell.h"
#include "core/vertex.h"
#include "core/edge.h"
#include "core/cellCluster.h"
#include "core/cellHelpers.h"
#include "core/edgeHelpers.h"
#include "core/vertexHelpers.h"

namespace epc {

template<int d>
class transitionHelpers {
public:
	transitionHelpers(){};
	
	Array<double,d> getCenterOfEdge(const std::shared_ptr<Edge<d>> &e);
	std::vector<std::shared_ptr<Cell<d>>> getCellsThatShareEdge(const CellCluster<d> &myCellCluster,const std::shared_ptr<Edge<d>> &e);
	std::vector<std::shared_ptr<Cell<d>>> getCellsThatShareVertex(const CellCluster<d> &myCellCluster,const std::shared_ptr<Vertex<d>> &v);

    void mergeEdgeAndUpdateCellCluster(	CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> &e, double damping);
    void t1OnEdgeAndUpdateCellCluster(	CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> &e, 
    									std::vector<std::shared_ptr<Cell<d>>> &cellsSharingEdge , 
    									std::vector<std::shared_ptr<Cell<d>>> &cellsContainingOnlyV1, 
    									std::vector<std::shared_ptr<Cell<d>>> &cellsContainingOnlyV2,
    									double dmin);
    void replaceEdgeByVertexInCell(	std::shared_ptr<Edge<d>> &e, 
    								std::shared_ptr<Vertex<d>> &v3, 
    								std::shared_ptr<Cell<d>> &cell);
    
	
	Array<double,d> getCenterOfEdge(std::shared_ptr<Vertex<d>> v1, std::shared_ptr<Vertex<d>> v2);
    void addVertexToEdgeInCorrectPosition(std::shared_ptr<Edge<d>> e, std::shared_ptr<Vertex<d>> v, std::shared_ptr<Vertex<d>> v3);
    void removeOldVertices(CellCluster<d> &myCellCluster, std::shared_ptr<Vertex<d>> v1, std::shared_ptr<Vertex<d>> v2);
    void changeEdgeOfMergedVertices(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> e, std::shared_ptr<Vertex<d>> v1, std::shared_ptr<Vertex<d>> v2, std::shared_ptr<Vertex<d>> v3);
    void changeConnectivityAfterTransition(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> e);
	void printCellSharingVertices(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> e, std::shared_ptr<Vertex<d>> v1, std::shared_ptr<Vertex<d>> v2);
	void printEdgesAtCell(std::shared_ptr<Cell<d>> &cell);
	void printVerticesAtCell(std::shared_ptr<Cell<d>> &cell);
};

template<int d>
void transitionHelpers<d>::t1OnEdgeAndUpdateCellCluster(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> &e, 
    									std::vector<std::shared_ptr<Cell<d>>> &cellsSharingEdge , 
    									std::vector<std::shared_ptr<Cell<d>>> &cellsContainingOnlyV1, 
    									std::vector<std::shared_ptr<Cell<d>>> &cellsContainingOnlyV2,
    									double dmin){
	
	// DEBUG MESSAGES
	// std::cout << " -- cellsSharingEdge -- " << std::endl;
	// for (auto c : cellsSharingEdge){
	// 	std::cout << c->getId() << std::endl;
	// }
	// std::cout << " -- cellsContainingOnlyV1 -- " << std::endl;
	// for (auto c : cellsContainingOnlyV1){
	// 	std::cout << c->getId() << std::endl;
	// }
	// std::cout << " -- cellsContainingOnlyV2 -- " << std::endl;
	// for (auto c : cellsContainingOnlyV2){
	// 	std::cout << c->getId() << std::endl;
	// }
	// END DEBUG MESSAGES

	auto v1 = e->getVertexOne();
	auto v2 = e->getVertexTwo();
	// Replace shared edge by one of the two vertices.
	replaceEdgeByVertexInCell(e, v1,cellsSharingEdge[0]);

	if (cellsSharingEdge.size() > 1){
		replaceEdgeByVertexInCell(e, v2,cellsSharingEdge[1]);
		
	    // std::cout << "---Removing Connection between c12("<< cellsSharingEdge[0]->getId() <<") and c21("<< cellsSharingEdge[1]->getId()<< 
     //            ") through e("<<e->getId()<< ") "<< std::endl;

		// Remove the neighborhood connection between the cells sharing the edge e. 
		bool removedConnect = myCellCluster.removeConnectivityBetweenCellsThroughEdge(cellsSharingEdge[0],cellsSharingEdge[1],e);
	} else {
		bool foundEdge = false;
		auto c1Edges = cellsContainingOnlyV1[0]->getEdges();
		auto c12Edges = cellsSharingEdge[0]->getEdges();
		// Find the edge of cellsContainingOnlyV1[0] that contains vertex v1 and is not shared with cellsSharingEdge[0]. 
		// Replace its vertex v1 by v2.
		for (unsigned ie1 = 0; ie1 < c1Edges.size(); ie1++ ){
			if ((c1Edges[ie1]->getVertexOne() == v1) ||(c1Edges[ie1]->getVertexTwo() == v1)){
				bool shared = false;
				for (unsigned ie12 = 0; ie12 < c12Edges.size(); ie12++ ){
					if (c1Edges[ie1] == c12Edges[ie12]){
						shared = true;
						break;
					}
				}
				if (not shared){
					if (c1Edges[ie1]->getVertexOne() == v1){
						c1Edges[ie1]->getVertexOne() = v2;
					} else {
						c1Edges[ie1]->getVertexTwo() = v2;
					}
					foundEdge = true;
				}
			} 
		}
		EPC_ASSERT((foundEdge == true) && "ERROR Did not find edge in c1 that contains v1 but is not shared with c12");
		// We are in the case where edge e WAS a boundary edge. After T1, v2 stays a boundary vertex, but v1 BECOMES a bulk vertex.
		myCellCluster.removeFromBoundaryVertices(v1);
		std::cout << " After T1, v1("<< v1->getId() <<") is removed from boundary vertices" << std::endl;
	}

	// Add the edge e to the cells that contained previously only v1.
	for (auto c : cellsContainingOnlyV1){
		c->getEdges().push_back(e);
	}
	// Add the edge e to the cells that contained previously only v2.
	for (auto c : cellsContainingOnlyV2){
		c->getEdges().push_back(e);
	}
	
	if (cellsContainingOnlyV1.size()==1 && cellsContainingOnlyV2.size()==1){

	 // std::cout << "----adding Connection between c12("<< cellsContainingOnlyV1[0]->getId() <<") and c21("<< cellsContainingOnlyV2[0]->getId()<< 
     //            ") through e("<<e->getId()<< ") "<< std::endl;

		myCellCluster.addConnectivityBetweenCellsThroughEdge(cellsContainingOnlyV1[0],cellsContainingOnlyV2[0],e);
		e->setLambda(myCellCluster.getCellTypesManager().getLambdaBetweenCellTypes(cellsContainingOnlyV1[0]->getType(),cellsContainingOnlyV2[0]->getType()));
	} else {
		// std::cout << " e->getLambda() = " << e->getLambda() ;
		e->setLambda(myCellCluster.getBoundaryLambda());
		// Keep the list of boundary vertices up-to-date. Add v1 and v2 to the list of boundary vertices, because e BECOMES a boundary edge. They are added to the list only if they are not already in it.
		myCellCluster.addToBoundaryVertices(v1);
		myCellCluster.addToBoundaryVertices(v2);
		std::cout << " After T1, v1("<< v1->getId() <<") and v2(" << v2->getId() << ") are added to boundary vertices" << std::endl;
		// std::cout << " ==> e->getLambda() = " << e->getLambda() << std::endl;
	}
	// Move v1 apart from v2.
	Array<double, d> vectV1V2 = e->getVertexOne()->getPosition() - e->getVertexTwo()->getPosition();
	v1->setPosition(v1->getPosition()+(1.1*dmin-e->computeLength())*(vectV1V2/vectV1V2.norm()));
	// std::cout << "After T1, cells sharing edge("<< e->getId() << ") are:" ;
	// for (auto c : getCellsThatShareEdge(myCellCluster, e)){
	// 	std::cout << c->getId() << ",";
	// }
	// std::cout << std::endl;

}

// Merge the vertices of a given edge and update the cell cluster correspondingly.
template<int d>
void transitionHelpers<d>::mergeEdgeAndUpdateCellCluster(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> &e, double damping){
//  To create new vertex
	std::shared_ptr<Vertex<d>> v1 = e->getVertexOne();
	std::shared_ptr<Vertex<d>> v2 = e->getVertexTwo();
	Array<double,2> posv3 = getCenterOfEdge(e);
	std::shared_ptr<Vertex<2>> v3 = std::shared_ptr<Vertex<2>>(new Vertex<2>(posv3,std::shared_ptr<Dynamics<2>>(new Verlet<2>(damping))));

//  Update all edges containing one of the two merged vertices
	//std::cout << "Merge v1("<< v1->getId() << ") and v2("<< v2->getId() << ") of e("<< e->getId()<<") into v3("<<v3.getId()<<")" << std::endl;	
	auto myedges = myCellCluster.getEdges();
	for(auto ie: myedges){
		auto iev1 = ie->getVertexOne();
		auto iev2 = ie->getVertexTwo();
		if ( ((iev1 == v1)&&(iev2 != v2)) || ((iev1 == v2)&&(iev2 != v1)) ){//we will replace vertex v1 by new vertex v3 in edge ie
			ie->getVertexOne() = v3;		
		}else if (((iev1 != v1)&&(iev2 == v2)) || ((iev1 != v2)&&(iev2 == v1)) ){//we will replace vertex v2 by new vertex v3 in edge ie
			ie->getVertexTwo() = v3;
		}
	}

	//	Remove old edge from cell edges
	auto mycells = myCellCluster.getCells();
	for(auto ic: mycells){
		for(auto cellie: ic->getEdges()){
			if(cellie->getId() == e->getId()){
				ic->removeEdge(e);
				break;
			}
		}				
	}

	//  Remove the edge and its vertices from the tissue
	myCellCluster.removeVertex(e->getVertexOne());
	myCellCluster.removeVertex(e->getVertexTwo());
	myCellCluster.removeEdge(e);
	// Add the new vertex.
	myCellCluster.getVertices().push_back(v3);

	// Update the list of boundary vertices.
	myCellCluster.removeFromBoundaryVertices(e->getVertexOne());
	myCellCluster.removeFromBoundaryVertices(e->getVertexTwo());
	myCellCluster.addToBoundaryVertices(v3);
}


// Get the position of the point at the center of a given edge.
template<int d>
Array<double,d> transitionHelpers<d>::getCenterOfEdge(const std::shared_ptr<Edge<d>> &e){
		std::shared_ptr<Vertex<d>> v1 = e->getVertexOne();
		std::shared_ptr<Vertex<d>> v2 = e->getVertexTwo();
		
		return 0.5*(v1->getPosition()+v2->getPosition());
}


template<int d>
Array<double,d> transitionHelpers<d>::getCenterOfEdge(std::shared_ptr<Vertex<d>> v1, std::shared_ptr<Vertex<d>> v2){
	double x1,x2,y1,y2;
	x1 = v1->getPosition()[0];x2 = v2->getPosition()[0];
	y1 = v1->getPosition()[1];y2 = v2->getPosition()[1];

	double m,b;
	m = (double) ((y2-y1)/(x2-x1));
	b = y1 - (double) (x1 * (double)m);
		
	double my_x = abs(x2+x1)/2.0; 
	double my_y = m*my_x + b;//y associated to that random x, considering the equation of the line

	Array<double,2> v3; 
	v3[0] = my_x;
	v3[1] = my_y;

	return(v3);
}

template<int d>
std::vector<std::shared_ptr<Cell<d>>> transitionHelpers<d>::getCellsThatShareEdge(const CellCluster<d> &myCellCluster, const std::shared_ptr<Edge<d>> &e){
	std::vector<std::shared_ptr<Cell<d>>> cells;
	for(auto c : myCellCluster.getCells()){
		for(auto ie : c->getEdges()){
			if(ie == e){
				cells.push_back(c);
				break;//leave the edges-loop	
			}	
		}
	}
	if ((cells.size()<1) || (cells.size()>2)){
		std::cout << "Number of cells sharing edge(" << e->getId() << ") coonecting vertices (" << e->getVertexOne()->getId() << " & " << e->getVertexTwo()->getId() << ") = "  << cells.size() << std::endl;
		for (int i = 0; i < cells.size() ; i++){
			std::cout << cells[i]->getId() << ", ";
		}
		std::cout << std::endl;
	}
	EPC_ASSERT(((cells.size()>=1) && (cells.size()<=2)) &&
                        "Error Cells sharing an edge should be 1 or 2.");
	return cells;
}

template<int d>
std::vector<std::shared_ptr<Cell<d>>> transitionHelpers<d>::getCellsThatShareVertex(const CellCluster<d> &myCellCluster,const std::shared_ptr<Vertex<d>> &v){
	std::vector<std::shared_ptr<Cell<d>>> cells;
	for(auto c : myCellCluster.getCells()){
		for(auto iv : cellHelpers<d>::getVertices(*c)){
			if(iv->getId() == v->getId()){
				cells.push_back(c);
				break;//leave the edges-loop	
			}	
		}
	}
	if ((cells.size()<1) || (cells.size()>3)){
		std::cout << "Number of cells sharing vertex(" << v->getId() << ") = "  << cells.size() << std::endl;
		for (int i = 0; i < cells.size() ; i++){
			std::cout << cells[i]->getId() << ", ";
		}
		std::cout << std::endl;
	}
	EPC_ASSERT(((cells.size()>=1) && (cells.size()<=3)) &&
                        "Error Cells sharing a vertex should be 1, 2 or 3.");
	return cells;
}

template<int d>
void transitionHelpers<d>::addVertexToEdgeInCorrectPosition(std::shared_ptr<Edge<d>> e, std::shared_ptr<Vertex<d>> v, std::shared_ptr<Vertex<d>> v3){
	if( (e->getVertexOne() == v) && (e->getVertexTwo() != v3) ){
		e->getVertexTwo() = v3;	
	}else if((e->getVertexTwo() == v) && (e->getVertexOne() != v3)) {
		e->getVertexOne() = v3;	
	}
	return;
}

template<int d>
void transitionHelpers<d>::replaceEdgeByVertexInCell(std::shared_ptr<Edge<d>> &e, std::shared_ptr<Vertex<d>> &v3, std::shared_ptr<Cell<d>> &cell){
	auto ev1 = e->getVertexOne();
	auto ev2 = e->getVertexTwo();
	auto edges = cell->getEdges();
	for(auto ce: edges){
		auto cev1 = ce->getVertexOne();
		auto cev2 = ce->getVertexTwo();
		if ( ((ev1 == cev1)&&(ev2 != cev2)) || ((ev2 == cev1)&&(ev1 != cev2)) ){//we will replace vertex v1 by new vertex v3 in edge ie
			ce->getVertexOne() = v3;
		}else if ( ((ev1 == cev2)&&(ev2 != cev1)) || ((ev2 == cev2)&&(ev1 != cev1)) ){//we will replace vertex v2 by new vertex v3 in edge ie
			ce->getVertexTwo() = v3;
		}
	}
	cell->removeEdge(e);
	return;	
}

template<int d>
void transitionHelpers<d>::changeEdgeOfMergedVertices(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> e, 
													std::shared_ptr<Vertex<d>> v1, std::shared_ptr<Vertex<d>> v2, std::shared_ptr<Vertex<d>> v3){
	auto myedges = myCellCluster.getEdges();
	for(auto ie: myedges){
		auto iev1 = ie->getVertexOne();
		auto iev2 = ie->getVertexTwo();
		if ( ((v1 == iev1)&&(v2 != iev2)) || ((v1 == iev2)&&(v2 != iev1)) ){//we will replace vertex v1 by new vertex v3 in edge ie
			addVertexToEdgeInCorrectPosition(ie,v1,v3);
		}else if ( ((v2 == iev1)&&(v1 != iev2)) || ((v2 == iev2)&&(v1 != iev1)) ){//we will replace vertex v2 by new vertex v3 in edge ie
			addVertexToEdgeInCorrectPosition(ie,v2,v3);
		}
	}
	return;	
}


template<int d>
void transitionHelpers<d>::changeConnectivityAfterTransition(CellCluster<d> &myCellCluster,std::shared_ptr<Edge<d>> e){

	return;
}
template<int d>
void transitionHelpers<d>::removeOldVertices(CellCluster<d> &myCellCluster, std::shared_ptr<Vertex<d>> v1, std::shared_ptr<Vertex<d>> v2){
	auto allvertices = myCellCluster.getVertices();
	int nremoved = 0;
	for(auto iv : allvertices){
		if(iv == v1){
			myCellCluster.removeVertex(v1);
			nremoved++;
		}else if(iv == v2){
			myCellCluster.removeVertex(v2);
			nremoved++;
		} 
		if(nremoved == 2){
			break;
		}
	}	
	return;
}

template<int d>
void transitionHelpers<d>::printEdgesAtCell(std::shared_ptr<Cell<d>> &cell){
	auto edges = cell->getEdges();
	std::cout << "List of Edges: " << edges.size() << " edges!" << std::endl;
	for(auto ie : edges){
		std::cout << ie->getVertexOne()->getPosition() << " " << ie->getVertexTwo()->getPosition() << std::endl; 
	}
	return;
}

template<int d>
void transitionHelpers<d>::printVerticesAtCell(std::shared_ptr<Cell<d>> &cell){
	auto vertices = cellHelpers<d>::getVertices(*cell);
	std::cout << "List of Vertices: " << vertices.size() << " vertices!" << std::endl;
	for(auto iv : vertices){
		std::cout << iv->getPosition() << " " << iv->getPosition() << std::endl; 
	}
	return;
}

template<int d>
void transitionHelpers<d>::printCellSharingVertices(CellCluster<d> &myCellCluster, std::shared_ptr<Edge<d>> e, std::shared_ptr<Vertex<d>> v1, std::shared_ptr<Vertex<d>> v2){
	auto allcells = myCellCluster.getCells();
	for(auto ic : allcells){
		auto poscell = std::find(allcells.begin(), allcells.end(),ic) - allcells.begin();
		auto nNeighCells = myCellCluster.getConnectivity()[poscell].size();
		for(auto ineigh = 0; ineigh < nNeighCells; ineigh++){
			if (std::get<1>(myCellCluster.getConnectivity()[poscell][ineigh])->getId() == e->getId()){
				this->printEdgesAtCell(ic);
				this->printVerticesAtCell(ic);
				return;
			}
		}
	}
	return;
}


}
#endif
