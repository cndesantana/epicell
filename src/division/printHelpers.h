#ifndef _PRINTHELPERS_
#define _PRINTHELPERS_

#include "mathematics/array.h"
#include "core/cell.h"
#include "core/vertex.h"
#include "core/edge.h"
#include "core/cellHelpers.h"
#include "core/vertexHelpers.h"
#include <vector> 
#include <iostream>
#include <algorithm>

namespace epc {

template<int d>
class PrintHelpers {
public:
	void printListOfVertices(std::vector<std::shared_ptr<Vertex<d>>> listOfVertices); 
	void printListOfEdges(std::vector<std::shared_ptr<Edge<d>>> listOfEdges); 
};

template<int d>
void PrintHelpers<d>::printListOfEdges(std::vector<std::shared_ptr<Edge<d>>> edges){

	std::cout << "************************************" << std::endl;
	std::cout << "LIST OF " << edges.size() << " EDGES" << std::endl;
	for (unsigned iE = 0; iE < edges.size(); ++iE) {
		std::cout << "e" << edges[iE]->getId() << " - " << edges[iE]->getVertexOne()->getPosition() << " ";
		std::cout << edges[iE]->getVertexTwo()->getPosition() << std::endl;
	}
	return;
}

template<int d>
void PrintHelpers<d>::printListOfVertices(std::vector<std::shared_ptr<Vertex<d>>> verts){
	unsigned long nVerts = verts.size();

	std::cout << "************************************" << std::endl;
	std::cout << "LIST OF VERTICES: EDGES FROM VERTICES: " << std::endl;
	for (unsigned iV = 0; iV < verts.size(); ++iV) {
		std::cout << verts[iV]->getPosition() << " " << verts[(iV+1)%nVerts]->getPosition() << std::endl;
	}
	return;
}
}//end namespace epd

#endif
