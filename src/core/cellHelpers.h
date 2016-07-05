#ifndef CELL_HELPERS_H
#define CELL_HELPERS_H

#include "vertex.h"
#include "vertexHelpers.h"
#include "cell.h"
#include "edge.h"
#include "edgeHelpers.h"
#include "mathematics/rng.h"

namespace epc {

template<int d> class Cell;    

template<int d>
inline bool operator==(Cell<d> const& c1, Cell<d> const& c2) {
    return (c1.getId() == c2.getId());
}

template<int d>
struct cellHelpers {
/// Checks if the edges composing a cell are forming a closed polygon
static bool checkCoherence(const Cell<d> &cell) {
    for ( auto edge : cell.getEdges() ) {
        const unsigned v1 = edge->getVertexOne()->getId();
        const unsigned v2 = edge->getVertexTwo()->getId();

        unsigned count1 = 0;
        unsigned count2 = 0;
        for ( auto const &edgeComp : cell.getEdges() ) {
            if (edge->getId() != edgeComp->getId()) {
                const unsigned v1comp = edgeComp->getVertexOne()->getId();
                const unsigned v2comp = edgeComp->getVertexTwo()->getId();

                if (v1 == v1comp) {
                    EPC_ASSERT(edge->getVertexOne()->getPosition() == edgeComp->getVertexOne()->getPosition() && 
                        "Error two vertex that must be on the same position are not");
                    count1++;
                }
                if (v2 == v2comp) {
                    count2++;
                    EPC_ASSERT(edge->getVertexTwo()->getPosition() == edgeComp->getVertexTwo()->getPosition() && 
                        "Error two vertex that must be on the same position are not");
                }
                if (v1 == v2comp) {
                    count1++;
                    EPC_ASSERT(edge->getVertexOne()->getPosition() == edgeComp->getVertexTwo()->getPosition() && 
                        "Error two vertex that must be on the same position are not");
                }
                if (v2 == v1comp) {
                    count2++;
                    EPC_ASSERT(edge->getVertexTwo()->getPosition() == edgeComp->getVertexOne()->getPosition() && 
                        "Error two vertex that must be on the same position are not:");
                }
            }
        }
        if (count1 != 1 || count2 != 1) {
            return false;
        }
    }
    return true;
}


/// Builds a cell from a list of vertices
static std::shared_ptr<Cell<d>> createCell(std::vector<std::shared_ptr<Vertex<d>>> &vertices, double K, double A0, 
    double gamma, double lambda, double sdGamma = 0.0, double sdLambda = 0.0, unsigned long iT = 0) 
{
    vertexHelpers<d>::sortVerticesCounterClockwiseInPlace(vertices); // sort edges counter clockwise
    std::vector<std::shared_ptr<Edge<d>>> edges;
    
    for (unsigned i = 0; i < vertices.size(); ++i) {
	
//        edges.push_back(std::shared_ptr<Edge<d>>(new Edge<d>(vertices[i],vertices[(i+1)%vertices.size()],NRNG(lambda,(lambda/100.0)).ngenerate()))); // create edges
        edges.push_back(std::shared_ptr<Edge<d>>(new Edge<d>(vertices[i],vertices[(i+1)%vertices.size()],lambda,sdLambda))); // create edges
    }

    
    auto cell = std::shared_ptr<Cell<d>>(new Cell<d>(edges, K, A0, gamma, sdGamma, iT));
    EPC_ASSERT(checkCoherence(*cell));
    return cell;
}

/// Removes duplicated edges in two cells
static bool removeDuplicatedEdge(std::shared_ptr<Cell<d>> const &cellOne, std::shared_ptr<Cell<d>> &cellTwo) {
    bool removed = false;
    for (auto &edgeOne: cellOne->getEdges()) {
        for (auto &edgeTwo : cellTwo->getEdges()) {
            if (*edgeOne == *edgeTwo) {
                edgeTwo.reset();
                edgeTwo = edgeOne;
                removed = true;
            }
        }
    }
    return removed;
}

/// Returns common edge between two cells
static std::shared_ptr<Edge<d>> getCommonEdge(std::shared_ptr<Cell<d>> const &cellOne, std::shared_ptr<Cell<d>> &cellTwo) {
    for (auto &edgeOne: cellOne->getEdges()) {
        for (auto &edgeTwo : cellTwo->getEdges()) {
            if (*edgeOne == *edgeTwo) {
                return edgeOne;
            }
        }
    }
    // if no common edge return 0 (for the moment raise an assert)
    EPC_ASSERT(false);
    return 0;
}

static std::vector<std::shared_ptr<Vertex<d>>> getVertices(Cell<d> const &cell) {
    std::vector<std::shared_ptr<Vertex<d>>> vertices;
    //TO SORT THE LIST OF EDGES MAY SOLVE THIS PROBLEM
    for (auto edge: cell.getEdges()) {
        vertexHelpers<d>::appendIfNotInVector(edge->getVertexOne(),vertices);
        vertexHelpers<d>::appendIfNotInVector(edge->getVertexTwo(),vertices);
    }
    if(vertices.size() != cell.getEdges().size()){
        auto nVert = vertices.size();
        std::cout << "Problem with c( " << cell.getId() << ")" << std::endl;
        std::cout << "Vertices: " << nVert << std::endl;
        for(auto iV=0;iV<(int)nVert;iV++){
            std::cout << vertices[iV]->getId() << " " << vertices[(iV+1)%nVert]->getId() << std::endl; 
        }
        auto nEdges = cell.getEdges().size();
        std::cout << "Edges = " << nEdges << std::endl; 
        for (auto edge: cell.getEdges()){
            std::cout << edge->getVertexOne()->getId() << " " << edge->getVertexTwo()->getId() << std::endl; 
        }
    }
    EPC_ASSERT(vertices.size() == cell.getEdges().size());
    return vertices;
}


};


}  // end namespace epc

#endif
