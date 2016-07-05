#ifndef EPITHELIUM_HELPERS_H
#define EPITHELIUM_HELPERS_H

#include "vertex.h"
#include "edge.h"
#include "cellCluster.h"
#include "mathematics/constants.h"

namespace epc {

template<unsigned d>
struct epitheliumHelpers {

// This function returns all the boundary cells (<=> cells with a number of edges higher than the number of neighbors)
static std::vector<std::shared_ptr<Cell<d>>> getBoundaryCells(const CellCluster<d> &tissue) 
{
    std::vector<std::shared_ptr<Cell<d>>> boundaryCells;
    const auto &cells = tissue.getCells();
    const auto &connectivity = tissue.getConnectivity();

    for (unsigned iC = 0; iC < cells.size(); ++iC) { // loop through all cells
        if (connectivity[iC].size() < cells[iC]->getEdges().size()) { // if a cell is connected to less cells than the number of its edges
            boundaryCells.push_back(cells[iC]);
        }
    }
    return boundaryCells;
}

/// get the cells inside a cercle
static std::vector<std::shared_ptr<Cell<d>>> getCellsInCircle(const CellCluster<d> &tissue, Array<double,d> center, double radius){
    std::vector<std::shared_ptr<Cell<d>>> cellsInCircle;
    for (auto c : tissue.getCells()){
        if ((c->computeCentroid()-center).norm() < radius){
            cellsInCircle.push_back(c);
        }
    }
    return cellsInCircle;
}

/// get the cells inside a rectangle
static std::vector<std::shared_ptr<Cell<d>>> getCellsInRectangle(const CellCluster<d> &tissue, Array<double,d> center, double length, double width){
    double nbCells = 0;
    double nbCellsIn = 0;
    std::vector<std::shared_ptr<Cell<d>>> cellsInRect;
    for (auto c : tissue.getCells()){
        nbCells += 1;
        if ((fabs(c->computeCentroid()[1]-center[1]) < 0.5*width) && (fabs(c->computeCentroid()[0]-center[0]) < 0.5*length)){
            cellsInRect.push_back(c);
            nbCellsIn += 1;
        }
    }
    std::cout << "nbCellTotal = " << nbCells << ", nbCellsInRect = " << nbCellsIn << std::endl;
    return cellsInRect;
}

/// get the cells below a line y = y_
static std::vector<std::shared_ptr<Cell<d>>> getCellsBelowLine(const CellCluster<d> &tissue, double y_){
    std::vector<std::shared_ptr<Cell<d>>> cellsBelowLine;
    for (auto c : tissue.getCells()){
        if (c->computeCentroid()[1] < y_){
            cellsBelowLine.push_back(c);
        }
    }
    return cellsBelowLine;
}


/// This function computes all the edges shared by cells of types 1 and type 2.
static std::vector<std::shared_ptr<Edge<d>>> getEdgesBetweenTypes(const CellCluster<d> &tissue, std::shared_ptr<CellType> type1, std::shared_ptr<CellType> type2){
    std::vector<std::shared_ptr<Edge<d>>> edges;
    const auto &cells = tissue.getCells();
    const auto &connectivity = tissue.getConnectivity();
    for (unsigned iC = 0; iC < cells.size(); ++iC) { // loop through all cells
        auto edgesIC(cells[iC]->getEdges()); // create a copy of the vector of pointers on the edges of the cell  
        for (unsigned iA = 0; iA < connectivity[iC].size(); ++iA) {
            if (((std::get<0>(connectivity[iC][iA])->getType().getId() == type1->getId()) && (cells[iC]->getType().getId() == type2->getId())) ||
                ((std::get<0>(connectivity[iC][iA])->getType().getId() == type2->getId()) && (cells[iC]->getType().getId() == type1->getId()))){
                for (unsigned iB = 0; iB < edgesIC.size(); ++iB) {
                    if (std::get<1>(connectivity[iC][iA])->getId() == edgesIC[iB]->getId()) { // remove the shared edge from the edge vector
                        edgeHelpers<d>::appendIfNotInVector(edgesIC[iB], edges); // add edges to the list of boundary edges (if not already in the vector?)        
                    }
                }
            }
        }
    }
    return edges;

}

/// This function computes the boundary edges of cells.
static std::vector<std::shared_ptr<Edge<d>>> getAllBoundaryEdges(
    const CellCluster<d> &tissue) 
{
    std::vector<std::shared_ptr<Edge<d>>> boundaryEdges;

    const auto &cells = tissue.getCells();
    const auto &connectivity = tissue.getConnectivity();

    for (unsigned iC = 0; iC < cells.size(); ++iC) { // loop through all cells
        if (connectivity[iC].size() < cells[iC]->getEdges().size()) { // if a cell is connected to less cells than the number of its edges
            auto edges(cells[iC]->getEdges()); // create a copy of the vector of pointers on the edges of the cell
            for (unsigned iA = 0; iA < connectivity[iC].size(); ++iA) {
                for (unsigned iB = 0; iB < edges.size(); ++iB) {
                    if (std::get<1>(connectivity[iC][iA])->getId() == edges[iB]->getId()) { // remove the shared edge from the edge vector
                        edges.erase(edges.begin()+iB);
                        break;
                    }
                }
            }
            // here we are left with an edge vector which only contains boundary edges
            EPC_ASSERT(edges.size()+connectivity[iC].size() == cells[iC]->getEdges().size() && "the number of boundary edges plus the ones connected to other cells must be the total number of the cell edges");
            EPC_ASSERT(edges.size() > 0 && "when connectivity is smaller than the number of cell edges there must be more boundary edges 0");
            for (auto edge: edges) {
                edgeHelpers<d>::appendIfNotInVector(edge, boundaryEdges); // add edges to the list of boundary edges (if not already in the vector?)
            }
        }
    }    
    return boundaryEdges;
}

/// This function returns the boundary vertices. A boundary vertex is a vertex which belongs to boundary edge.
static std::vector<std::shared_ptr<Vertex<d>>> getAllBoundaryVertices(
    const CellCluster<d> &tissue) 
{
    std::vector<std::shared_ptr<Edge<d>>> boundaryEdges = getAllBoundaryEdges(tissue);

    std::vector<std::shared_ptr<Vertex<d>>> bcVertices;

    for (auto edge: boundaryEdges) {
        vertexHelpers<d>::appendIfNotInVector(edge->getVertexOne(), bcVertices); // add vertex 1 to bc vertices if not already in the vector
        vertexHelpers<d>::appendIfNotInVector(edge->getVertexTwo(), bcVertices); // add vertex 2 to bc vertices if not already in the vector
    }
    return bcVertices;
}


/// This function computes the boundary vertices of cells with less than numBulkNeighbors neighbors.
/// Gets all the vertices on cells with less neighbors than the number edges
static bool isVertexOnBoundary(
    const CellCluster<d> &tissue, const std::shared_ptr<Vertex<d>> &v) 
{
    std::vector<std::shared_ptr<Vertex<d>>> bcVertices;

    const auto &cells = tissue.getCells();
    const auto &connectivity = tissue.getConnectivity();

    for (unsigned iC = 0; iC < cells.size(); ++iC) { // loop through all cells
        if (connectivity[iC].size() < cells[iC]->getEdges().size()) { // if a cell is connected to less than "num" other cells ... (it means it's a boundary cell)
            auto edges(cells[iC]->getEdges()); // create a copy of the vector of pointers on the edges of the cell
            for (unsigned iA = 0; iA < connectivity[iC].size(); ++iA) {
                for (unsigned iB = 0; iB < edges.size(); ++iB) {
                    if (std::get<1>(connectivity[iC][iA])->getId() == edges[iB]->getId()) { // remove the shared edge from the edge vector
                        edges.erase(edges.begin()+iB);
                        break;
                    }
                }
            }
            // here we are left with an edge vector which only contains boundary edges
            EPC_ASSERT(edges.size()+connectivity[iC].size() == cells[iC]->getEdges().size() && "the bc edge number plus the ones connected to other cells must be 'num'");
            EPC_ASSERT(edges.size() > 0 && "when connectivity is smaller than 'num' there must be more boundary edges 0");
            for (auto edge: edges) {
                if ((v->getId() == edge->getVertexOne()->getId()) || (v->getId() == edge->getVertexTwo()->getId())){
                    return true;
                }
            }
        }
    }    
    return false;
}

/// This function computes the boundary vertices of cells with less than numBulkNeighbors neighbors.
/// Gets all the vertices on cells with less neighbors than the number 
/// numBulkNeighbors: number of neighbors for bulk cells (hexagons -> 6, rectangle -> 4)
static std::vector<std::shared_ptr<Vertex<d>>> getAllBoundaryVertices(
    const CellCluster<d> &tissue, unsigned numBulkNeighbors) 
{
    std::vector<std::shared_ptr<Vertex<d>>> bcVertices;

    const auto &cells = tissue.getCells();
    const auto &connectivity = tissue.getConnectivity();

    for (unsigned iC = 0; iC < cells.size(); ++iC) { // loop through all cells
        if (connectivity[iC].size() < numBulkNeighbors) { // if a cell is connected to less than "num" other cells ... (it means it's a boundary cell)
            auto edges(cells[iC]->getEdges()); // create a copy of the vector of pointers on the edges of the cell
            for (unsigned iA = 0; iA < connectivity[iC].size(); ++iA) {
                for (unsigned iB = 0; iB < edges.size(); ++iB) {
                    if (std::get<1>(connectivity[iC][iA])->getId() == edges[iB]->getId()) { // remove the shared edge from the edge vector
                        edges.erase(edges.begin()+iB);
                        break;
                    }
                }
            }
            // here we are left with an edge vector which only contains boundary edges
            EPC_ASSERT(edges.size()+connectivity[iC].size() == numBulkNeighbors && "the bc edge number plus the ones connected to other cells must be 'num'");
            EPC_ASSERT(edges.size() > 0 && "when connectivity is smaller than 'num' there must be more boundary edges 0");
            for (auto edge: edges) {
                vertexHelpers<d>::appendIfNotInVector(edge->getVertexOne(), bcVertices); // add vertex 1 to bc vertices if not already in the vector
                vertexHelpers<d>::appendIfNotInVector(edge->getVertexTwo(), bcVertices); // add vertex 2 to bc vertices if not already in the vector
            }
        }
    }    
    return bcVertices;
}

/// This function computes the boundary vertices of cells with exactly numBoundaryNeighbors neighbors.
/// numBulkNeighbors: number of neighbors for bulk cells (hexagons -> 6, rectangle -> 4)
/// numBoundaryNeighbors: number of neighbors for "bulk" boundary cells (numBoundaryNeighbors = numBulkNeighbors-1 )
static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryVertices(
    const CellCluster<d> &tissue, unsigned numBulkNeighbors, unsigned numBoundaryNeighbors) 
{
    std::vector<std::shared_ptr<Vertex<d>>> bcVertices;

    const auto &cells = tissue.getCells();
    const auto &connectivity = tissue.getConnectivity();

    for (unsigned iC = 0; iC < cells.size(); ++iC) { // loop through all cells
        if (connectivity[iC].size() == numBoundaryNeighbors) { // if a cell is connected to less than "num" other cells ... (it means it's a boundary cell)
            auto edges(cells[iC]->getEdges()); // create a copy of the vector of pointers on the edges of the cell
            for (unsigned iA = 0; iA < connectivity[iC].size(); ++iA) {
                for (unsigned iB = 0; iB < edges.size(); ++iB) { // for all the edges
                    if (std::get<1>(connectivity[iC][iA])->getId() == edges[iB]->getId()) { // remove the shared edge from the edge vector
                        edges.erase(edges.begin()+iB);
                        break;
                    }
                }
            }
            // here we are left with an edge vector which only contains boundary edges
            EPC_ASSERT(edges.size()+connectivity[iC].size() == numBulkNeighbors && "the bc edge number plus the ones connected to other cells must be 'num'");
            EPC_ASSERT(edges.size() > 0 && "when connectivity is smaller than 'num' there must be more boundary edges 0");
            for (auto edge: edges) {
                vertexHelpers<d>::appendIfNotInVector(edge->getVertexOne(), bcVertices); // add vertex 1 to bc vertices if not already in the vector
                vertexHelpers<d>::appendIfNotInVector(edge->getVertexTwo(), bcVertices); // add vertex 2 to bc vertices if not already in the vector
            }
        }
    }    
    return bcVertices;
}

/// This function computes the boundary vertices of cells with exactly numBoundaryNeighbors neighbors.
/// numBulkNeighbors: number of neighbors for bulk cells (hexagons -> 6, rectangle -> 4)
/// numBoundaryNeighbors: number of neighbors for "bulk" boundary cells (numBoundaryNeighbors = numBulkNeighbors-1 )
static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryVertices(
    const CellCluster<d> &tissue, unsigned numBulkNeighbors, unsigned numBoundaryNeighbors, const std::vector<std::shared_ptr<Vertex<d>>> &otherVertices) 
{
    std::vector<std::shared_ptr<Vertex<d>>> bcVertices;

    const auto &cells = tissue.getCells();
    const auto &connectivity = tissue.getConnectivity();

    for (unsigned iC = 0; iC < cells.size(); ++iC) { // loop through all cells
        if (connectivity[iC].size() == numBoundaryNeighbors) { // if a cell is connected to less than "num" other cells ... (it means it's a boundary cell)
            auto edges(cells[iC]->getEdges()); // create a copy of the vector of pointers on the edges of the cell
            for (unsigned iA = 0; iA < connectivity[iC].size(); ++iA) {
                for (unsigned iB = 0; iB < edges.size(); ++iB) {
                    if (std::get<1>(connectivity[iC][iA])->getId() == edges[iB]->getId()) { // remove the shared edge from the edge vector
                        edges.erase(edges.begin()+iB);
                        break;
                    }
                }
            }
            // here we are left with an edge vector which only contains boundary edges
            EPC_ASSERT(edges.size()+connectivity[iC].size() == numBulkNeighbors && "the bc edge number plus the ones connected to other cells must be 'num'");
            EPC_ASSERT(edges.size() > 0 && "when connectivity is smaller than 'num' there must be more boundary edges 0");
            for (auto edge: edges) {
                if (!vertexHelpers<d>::isInVector(*edge->getVertexOne(),otherVertices)) {
                    vertexHelpers<d>::appendIfNotInVector(edge->getVertexOne(), bcVertices); // add vertex 1 to bc vertices if not already in the vector
                }
                if (!vertexHelpers<d>::isInVector(*edge->getVertexTwo(),otherVertices)) {
                    vertexHelpers<d>::appendIfNotInVector(edge->getVertexTwo(), bcVertices); // add vertex 2 to bc vertices if not already in the vector
                }
            }
        }
    }    
    return bcVertices;
}

/*static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMinimumPositionInDirection(const CellCluster<d> &tissue, 
    unsigned dir, double dist) 
{
    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getAllBoundaryVertices(tissue);
    EPC_ASSERT(allBcVertices.size() > 0 && "Error no bc vertices.");
    auto minV = *std::min_element(allBcVertices.begin(), allBcVertices.end(), 
        [=](const std::shared_ptr<Vertex<d>> &v1, const std::shared_ptr<Vertex<d>> &v2) {
                                 return v1->getPosition()[dir] < v2->getPosition()[dir]; // one wants the smaller value per direction "dir"
                             });

    double minVal = minV->getPosition()[dir];
    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (fabs(v->getPosition()[dir]-minVal) < dist)
                minVertices.push_back(v);
        });

    return minVertices;
}

static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMaximumPositionInDirection(const CellCluster<d> &tissue, 
    unsigned dir, double dist) 
{
    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getAllBoundaryVertices(tissue);
    EPC_ASSERT(allBcVertices.size() > 0 && "Error no bc vertices.");
    auto maxV = *std::max_element(allBcVertices.begin(), allBcVertices.end(), 
        [=](const std::shared_ptr<Vertex<d>> &v1, const std::shared_ptr<Vertex<d>> &v2) {
                                 return v1->getPosition()[dir] < v2->getPosition()[dir]; // one wants the smaller value per direction "dir"
                             });

    double maxVal = maxV->getPosition()[dir];
    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (fabs(v->getPosition()[dir]-maxVal) < dist)
                minVertices.push_back(v);
        });

    return minVertices;
}*/

static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMinimumPositionInDirection(const CellCluster<d> &tissue, 
    unsigned numBulkNeighbors, unsigned dir, double dist) 
{
    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getAllBoundaryVertices(tissue,numBulkNeighbors);
    EPC_ASSERT(allBcVertices.size() > 0 && "Error no bc vertices.");
    auto minV = *std::min_element(allBcVertices.begin(), allBcVertices.end(), 
        [=](const std::shared_ptr<Vertex<d>> &v1, const std::shared_ptr<Vertex<d>> &v2) {
                                 return v1->getPosition()[dir] < v2->getPosition()[dir]; // one wants the smaller value per direction "dir"
                             });

    double minVal = minV->getPosition()[dir];
    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (fabs(v->getPosition()[dir]-minVal) < dist)
                minVertices.push_back(v);
        });

    return minVertices;
}

// Get all the boundary vertices with a minimum position (smaller than the boundary cell centroid of minimum position) in a given direction. 
static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMinimumPositionInDirection(const CellCluster<d> &tissue, unsigned dir, double dist) 
{
    std::vector<std::shared_ptr<Cell<d>>> boundaryCells = getBoundaryCells(tissue);
    auto minC = *std::min_element(boundaryCells.begin(), boundaryCells.end(), 
        [=](const std::shared_ptr<Cell<d>> &c1, const std::shared_ptr<Cell<d>> &c2) {
                                 return c1->computeCentroid()[dir] < c2->computeCentroid()[dir]; // one wants the smaller value per direction "dir"
                             });

    double minCentr = minC->computeCentroid()[dir];

    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getAllBoundaryVertices(tissue);

    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (v->getPosition()[dir] < minCentr + dist)
                minVertices.push_back(v);
        });

    return minVertices;
}

// Get all the boundary vertices with a maximum position (smaller than the centroid with minimum position) in a given direction. 
static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMaximumPositionInDirection(const CellCluster<d> &tissue, unsigned dir, double dist) 
{
    std::vector<std::shared_ptr<Cell<d>>> boundaryCells = getBoundaryCells(tissue);
    auto maxC = *std::max_element(boundaryCells.begin(), boundaryCells.end(), 
        [=](const std::shared_ptr<Cell<d>> &c1, const std::shared_ptr<Cell<d>> &c2) {
                                 return c1->computeCentroid()[dir] < c2->computeCentroid()[dir]; // one wants the smaller value per direction "dir"
                             });

    double maxCentr = maxC->computeCentroid()[dir];

    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getAllBoundaryVertices(tissue);

    std::vector<std::shared_ptr<Vertex<d>>> maxVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (v->getPosition()[dir] > maxCentr - dist)
                maxVertices.push_back(v);
        });

    return maxVertices;
}

static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMaximumPositionInDirection(const CellCluster<d> &tissue, 
    unsigned numBulkNeighbors, unsigned dir, double dist) 
{
    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getAllBoundaryVertices(tissue,numBulkNeighbors);
    EPC_ASSERT(allBcVertices.size() > 0 && "Error no bc vertices.");
    auto maxV = *std::max_element(allBcVertices.begin(), allBcVertices.end(), 
        [=](const std::shared_ptr<Vertex<d>> &v1, const std::shared_ptr<Vertex<d>> &v2) {
                                 return v1->getPosition()[dir] < v2->getPosition()[dir]; // one wants the smaller value per direction "dir"
                             });

    double maxVal = maxV->getPosition()[dir];
    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (fabs(v->getPosition()[dir]-maxVal) < dist)
                minVertices.push_back(v);
        });

    return minVertices;
}

static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMinimumPositionInDirection(const CellCluster<d> &tissue, 
    unsigned numBulkNeighbors, unsigned numBoundaryNeighbors, unsigned dir, double dist) 
{
    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getBoundaryVertices(tissue,numBulkNeighbors,numBoundaryNeighbors);
    EPC_ASSERT(allBcVertices.size() > 0 && "Error no bc vertices.");
    auto minV = *std::min_element(allBcVertices.begin(), allBcVertices.end(), 
        [=](const std::shared_ptr<Vertex<d>> &v1, const std::shared_ptr<Vertex<d>> &v2) {
                                 return v1->getPosition()[dir] < v2->getPosition()[dir]; // one wants the smaller value per direction "dir"
                             });

    double minVal = minV->getPosition()[dir];
    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (fabs(v->getPosition()[dir]-minVal) < dist)
                minVertices.push_back(v);
        });

    return minVertices;
}

static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMaximumPositionInDirection(const CellCluster<d> &tissue, 
    unsigned numBulkNeighbors, unsigned numBoundaryNeighbors, unsigned dir, double dist) 
{
    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getBoundaryVertices(tissue,numBulkNeighbors,numBoundaryNeighbors);
    EPC_ASSERT(allBcVertices.size() > 0 && "Error no bc vertices.");
    auto maxV = *std::max_element(allBcVertices.begin(), allBcVertices.end(), 
        [=](const std::shared_ptr<Vertex<d>> &v1, const std::shared_ptr<Vertex<d>> &v2) {
                                 return v1->getPosition()[dir] < v2->getPosition()[dir]; // one wants the smaller value per direction "dir"
                             });

    double maxVal = maxV->getPosition()[dir];
    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (fabs(v->getPosition()[dir]-maxVal) < dist)
                minVertices.push_back(v);
        });

    return minVertices;
}


static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMinimumPositionInDirection(const CellCluster<d> &tissue, 
    unsigned numBulkNeighbors, unsigned numBoundaryNeighbors, 
    unsigned dir, std::vector<std::shared_ptr<Vertex<d>>> &vertices, double dist) 
{
    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getBoundaryVertices(tissue,numBulkNeighbors,numBoundaryNeighbors,vertices);
    EPC_ASSERT(allBcVertices.size() > 0 && "Error no bc vertices.");
    auto minV = *std::min_element(allBcVertices.begin(), allBcVertices.end(), 
        [=](const std::shared_ptr<Vertex<d>> &v1, const std::shared_ptr<Vertex<d>> &v2) {
                                 return v1->getPosition()[dir] < v2->getPosition()[dir]; // one wants the smaller value per direction "dir"
                             });

    double minVal = minV->getPosition()[dir];
    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (fabs(v->getPosition()[dir]-minVal) < dist)
                minVertices.push_back(v);
        });

    return minVertices;
}

static std::vector<std::shared_ptr<Vertex<d>>> getBoundaryWithMaximumPositionInDirection(const CellCluster<d> &tissue, 
    unsigned numBulkNeighbors, unsigned numBoundaryNeighbors, unsigned dir, std::vector<std::shared_ptr<Vertex<d>>> &vertices, double dist) 
{
    std::vector<std::shared_ptr<Vertex<d>>> allBcVertices = getBoundaryVertices(tissue,numBulkNeighbors,numBoundaryNeighbors,vertices);
    EPC_ASSERT(allBcVertices.size() > 0 && "Error no bc vertices.");
    auto maxV = *std::max_element(allBcVertices.begin(), allBcVertices.end(), 
        [=](const std::shared_ptr<Vertex<d>> &v1, const std::shared_ptr<Vertex<d>> &v2) {
                                 return v1->getPosition()[dir] < v2->getPosition()[dir]; // one wants the smaller value per direction "dir"
                             });

    double maxVal = maxV->getPosition()[dir];
    std::vector<std::shared_ptr<Vertex<d>>> minVertices;
    std::for_each(allBcVertices.begin(), allBcVertices.end(), 
        [&](const std::shared_ptr<Vertex<d>> &v){ 
            if (fabs(v->getPosition()[dir]-maxVal) < dist)
                minVertices.push_back(v);
        });

    return minVertices;
}


}; // end epithelium helper

}  // end namespace epc

#endif
