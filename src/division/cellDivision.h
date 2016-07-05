#ifndef CELL_DIVISION_H
#define CELL_DIVISION_H

#include "core/epcDebug.h"
#include "core/cell.h"
#include "core/vertex.h"
#include "core/edge.h"
#include "edgeSelection.h"
#include "cellDivisionHelpers.h"
#include "mathematics/array.h"
#include "mathematics/rng.h"

namespace epc {

template<int d>
class CellDivision {
public:
	virtual ~CellDivision(){};

	virtual void divide(CellCluster<d> &tissue,std::shared_ptr<Cell<d>> &cell,double defaultLambda,double defaultSdLambda,double A0,double verletPar,unsigned long iT) = 0;
};

// =========== ShortestAxisCellDivision ====================== //
template<int d>
class ShortestAxisCellDivision final : public CellDivision<d> {
public:
	virtual void divide(CellCluster<d> &tissue,std::shared_ptr<Cell<d>> &cell,double defaultLambda,double defaultSdLambda,double A0,double verletPar,unsigned long iT) override;
};

template<int d>
void ShortestAxisCellDivision<d>::divide(CellCluster<d> &tissue,std::shared_ptr<Cell<d>> &cell,double defaultLambda,double defaultSdLambda,double A0,double verletPar,unsigned long iT) {
	double minDist;
	Array<double,d> final_v3_coord;
	Array<double,d> final_v4_coord;
	int final_posFirstEdge;
	int final_posOppositeEdge;

	std::vector<std::shared_ptr<Vertex<d>>> cellvertices = cellHelpers<d>::getVertices(*cell);
	int nVerts = cellvertices.size();
	CellDivisionHelpers<d> helper;
	for (int posFirstEdge = 0; posFirstEdge < cell->computeDegree(); posFirstEdge ++){
		
		std::shared_ptr<Edge<d>> firstEdge = cell->getEdges()[posFirstEdge];
		std::shared_ptr<Vertex<d>> firstEdgeV1 = firstEdge->getVertexOne();
		std::shared_ptr<Vertex<d>> firstEdgeV2 = firstEdge->getVertexTwo();

		// Sort Cell Vertices Counter Clockwise.
		vertexHelpers<d>::sortVerticesCounterClockwiseInPlace(cellvertices);
		int indexV1 = std::find(cellvertices.begin(), cellvertices.end(), firstEdgeV1) - cellvertices.begin();
		int indexV2 = std::find(cellvertices.begin(), cellvertices.end(), firstEdgeV2) - cellvertices.begin();
		// Minimum vertex (Counter Clockwise) of firstEdge.
		int minIndexV;
		if ((indexV1+1)%nVerts == indexV2){
			minIndexV = indexV1;
		} else {
			if ((indexV2+1)%nVerts == indexV1){
				minIndexV = indexV2;
			}
		}

		int indexOppositeVert;
		// Find minimum vertex of secondEdge.		
		if ((nVerts % 2) == 0){ //even number of vertices (considering v3)
			indexOppositeVert = (minIndexV + nVerts/2)%nVerts;
		}
		else{//odd number of edges
			indexOppositeVert = (minIndexV + (nVerts+1)/2)%nVerts;
		}

		//Find Opposite edge posOppositeEdge to posFirstEdge.
		std::shared_ptr<Edge<d>> secondEdge = helper.getEdgesByVertices(cell->getEdges(),cellvertices[indexOppositeVert],cellvertices[(indexOppositeVert+1)%nVerts]);
		int posOppositeEdge = std::find(cell->getEdges().begin(),cell->getEdges().end(),secondEdge) - cell->getEdges().begin();
		
		//Find coord of v3.
		Array<double,d> v3_coord = helper.getCenterOfEdge(firstEdge);
		
		//Find coord of v4.
		Array<double,d> v4_coord = helper.getCenterOfEdge(secondEdge);

		//If distance v3-v4 is smaller than the min distance found so far, keep v3, v4, posFirstEdge and posOppositeEdge.
		double v3_v4_dist = (v3_coord-v4_coord).norm();
		if ((posFirstEdge == 0) or (v3_v4_dist < minDist)){
			minDist = v3_v4_dist;
			final_v3_coord = v3_coord;
			final_v4_coord = v4_coord;
			final_posFirstEdge = posFirstEdge;
			final_posOppositeEdge = posOppositeEdge;
		}
	}

	// Add the vertices v3 and v4 to the cell vertices.
	std::shared_ptr<Vertex<d>> final_v3 = std::shared_ptr<Vertex<d>>(new Vertex<d>(final_v3_coord,std::shared_ptr<Dynamics<d>>(new Verlet<d>(verletPar))));
	std::shared_ptr<Vertex<d>> final_v4 = std::shared_ptr<Vertex<d>>(new Vertex<d>(final_v4_coord,std::shared_ptr<Dynamics<d>>(new Verlet<d>(verletPar))));
	cellvertices.push_back(final_v3);
	cellvertices.push_back(final_v4);

	// Split the cell.	
	vertexHelpers<d>::sortVerticesCounterClockwiseInPlace(cellvertices);
	helper.splitCell(tissue,cell,final_v3,final_v4,defaultLambda,defaultSdLambda,final_posFirstEdge,final_posOppositeEdge,A0,iT);
	return;
}

// =========== FarhadifarCellDivision ====================== //
template<int d>
class FarhadifarCellDivision final : public CellDivision<d> {
public:
	virtual void divide(CellCluster<d> &tissue,std::shared_ptr<Cell<d>> &cell,double defaultLambda,double defaultSdLambda,double A0,double verletPar,unsigned long iT) override;
};

template<int d>
void FarhadifarCellDivision<d>::divide(CellCluster<d> &tissue,std::shared_ptr<Cell<d>> &cell,double defaultLambda,double defaultSdLambda,double A0,double verletPar,unsigned long iT) {
//	RandomEdgeSelection<2> rndseledge;
//	Array<double,d> cent = cell->computeCentroid();
////	int firstEdge = rndseledge.select(cell);
////	auto e = cell->getEdges()[firstEdge];//randomly choose an edge in the selected cell, that will be crossed by the new edge that will divide the cell
////	auto v1 = e->getVertexOne();
////	auto v2 = e->getVertexTwo();
////	std::cout << v1->getPosition() << std::endl;
////	CellDivisionHelpers<2> helper;		
////	auto v3 = helper.randomPointInRange(v1->getPosition()[0],v1->getPosition()[1],v2->getPosition()[0],v2->getPosition()[1]);//randomly choose a point in the selected edge to be one vertex of the new edge
////	auto v4 = helper.intersectionLineEdge(v3[0],v3[1],cent[0],cent[1]);//v4 should be a vertex that is the intersection between one of the other edges of the cell and the line that crosses v3 and the centroid of the cell  -  TODO
////
	return;
// implement the division of the cell
}

// =========== RandomEqualCellDivision ====================== //
template<int d>
class RandomEqualCellDivision final : public CellDivision<d> {
public:
	virtual void divide(CellCluster<d> &tissue,std::shared_ptr<Cell<d>> &cell,double defaultLambda,double defaultSdLambda,double A0,double verletPar,unsigned long iT) override;
};

template<int d>
void RandomEqualCellDivision<d>::divide(CellCluster<d> &tissue,std::shared_ptr<Cell<d>> &cell,double defaultLambda,double defaultSdLambda,double A0,double verletPar,unsigned long iT) {
	RandomEdgeSelection<2> rndseledge;
	int posFirstEdge = rndseledge.select(cell);

	CellDivisionHelpers<2> helper;
	auto v3 = helper.getCentralVertexInEdge(tissue,cell->getEdges()[posFirstEdge]->getVertexOne(),cell->getEdges()[posFirstEdge]->getVertexTwo(),verletPar);

	int posOppositeVert;
	auto cellvertices = cellHelpers<2>::getVertices(*cell);
	cellvertices.push_back(v3);
	vertexHelpers<d>::sortVerticesCounterClockwiseInPlace(cellvertices);
	auto posv3 = std::find(cellvertices.begin(), cellvertices.end(), v3) - cellvertices.begin();
	auto nVerts = cellvertices.size();
	if ((nVerts % 2) == 0){ //even number of vertices (considering v3)
		posOppositeVert = (posv3 + nVerts/2)%nVerts;
	}
	else{//odd number of edges
		posOppositeVert = (posv3 + (nVerts-1)/2)%nVerts;
	}

	auto v4 = helper.getCentralVertexInEdge(tissue,cellvertices[posOppositeVert],cellvertices[(posOppositeVert+1)%nVerts],verletPar);
	auto e22 = helper.getEdgesByVertices(cell->getEdges(),cellvertices[posOppositeVert],cellvertices[(posOppositeVert+1)%nVerts]);
	int posOppositeEdge = std::find(cell->getEdges().begin(),cell->getEdges().end(),e22) - cell->getEdges().begin();
	cellvertices.push_back(v4);
	vertexHelpers<d>::sortVerticesCounterClockwiseInPlace(cellvertices);
	helper.splitCell(tissue,cell,v3,v4,defaultLambda,defaultSdLambda,posFirstEdge,posOppositeEdge,A0,iT);
	return;
}

}  // end namespace epc

#endif
