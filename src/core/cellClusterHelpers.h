#ifndef CELL_CLUSTER_HELPERS_H
#define CELL_CLUSTER_HELPERS_H

#include "mathematics/constants.h"
#include "cellHelpers.h"
#include "edgeHelpers.h"
#include "io/io.h"
#include <cmath>

namespace epc {

template<int d> struct cellClusterHelpersImpl;

template<int d>
struct cellClusterHelpers {

static void createRectangularPavementOfHexagonalCells(int nx, int ny, double circumRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, 
    std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
    std::vector<std::shared_ptr<Cell<d>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
    std::shared_ptr<Dynamics<d>> &dyn, bool isLarge,unsigned long iT) 
{
    cellClusterHelpersImpl<d>::createRectangularPavementOfHexagonalCells(nx, ny, circumRadius, K, gamma, lambda, sdGamma, sdLambda, vertices, edges, cells, connectivity, dyn, isLarge,iT);
}

static void createRectangularPavementOfHexagonalCells(int nx, int ny, double initCircumRadius, double K, double A0, double gamma, double lambda, 
    std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
    std::vector<std::shared_ptr<Cell<d>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
    std::shared_ptr<Dynamics<d>> &dyn, bool isLarge) 
{
    cellClusterHelpersImpl<d>::createRectangularPavementOfHexagonalCells(nx, ny, initCircumRadius, K, A0, gamma, lambda,vertices, edges, cells, connectivity, dyn, isLarge);
}

static void createCircularPavementOfHexagonalCells(int n, double circumRadius, double K, double gamma, double lambda,
    std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
    std::vector<std::shared_ptr<Cell<d>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
    std::shared_ptr<Dynamics<d>> &dyn) 
{
    cellClusterHelpersImpl<d>::createCircularPavementOfHexagonalCells(n, circumRadius, K, gamma, lambda,vertices, edges, cells, connectivity, dyn);
}

static void createRectangularPavementOfSquareCells(int nx, int ny, double initLength, double K, double A0, double gamma, double lambda,
    std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
    std::vector<std::shared_ptr<Cell<d>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
    std::shared_ptr<Dynamics<d>> &dyn) 
{
    cellClusterHelpersImpl<d>::createRectangularPavementOfSquareCells(nx, ny, initLength, K, A0, gamma, lambda,vertices, edges, cells, connectivity, dyn);
}

static void createRectangularPavementOfSquareCells(int nx, int ny, double length, double K, double gamma, double lambda, double sdGamma, double sdLambda, 
    std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
    std::vector<std::shared_ptr<Cell<d>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
    std::shared_ptr<Dynamics<d>> &dyn,unsigned long iT) 
{
    cellClusterHelpersImpl<d>::createRectangularPavementOfSquareCells(nx, ny, length, K, gamma, lambda, sdGamma, sdLambda, vertices, edges, cells, connectivity, dyn,iT);
}


static void createRingOfTrapezoides(unsigned nCells, double cellHeight, double ringRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, 
    std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
    std::vector<std::shared_ptr<Cell<d>>> &cells, std::shared_ptr<Dynamics<d>> &dyn,unsigned long iT) 
{
    cellClusterHelpersImpl<d>::createRingOfTrapezoides(nCells,cellHeight, ringRadius, K, gamma, lambda, sdGamma, sdLambda, vertices, edges, cells, dyn,iT);
}

static void createRingOfTrapezoides(int nCells, double cellHeight, double ringRadius, double K, double gamma, double lambda,
    std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
    std::vector<std::shared_ptr<Cell<d>>> &cells,
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
    std::shared_ptr<Dynamics<d>> &dyn) 
{
    cellClusterHelpersImpl<d>::createRingOfTrapezoides(nCells,cellHeight, ringRadius, K, gamma, lambda, vertices, edges, cells, connectivity, dyn);
}

static void createCircularPavementOfSquareCells(int n, double length, double K, double gamma, double lambda,
    std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
    std::vector<std::shared_ptr<Cell<d>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
    std::shared_ptr<Dynamics<d>> &dyn) 
{
    cellClusterHelpersImpl<d>::createCircularPavementOfSquareCells(n, length, K, gamma, lambda,vertices, edges, cells, connectivity, dyn);
}


};

template<int d>
struct cellClusterHelpersImpl {
    static void createRectangularPavementOfHexagonalCells(int nx, int ny, 
        double circumRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, 
        std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
        std::vector<std::shared_ptr<Cell<d>>> &cells, 
        std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
        std::shared_ptr<Dynamics<d>> &dyn,unsigned long iT)  
    {
        EPC_ASSERT(false && "No generic implementation for rectangular hexagonal pavement exists.");
    }

    static void createRectangularPavementOfHexagonalCells(int nx, int ny, 
        double initCircumRadius, double K, double A0, double gamma, double lambda,
        std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
        std::vector<std::shared_ptr<Cell<d>>> &cells, 
        std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
        std::shared_ptr<Dynamics<d>> &dyn)  
    {
        EPC_ASSERT(false && "No generic implementation for rectangular hexagonal pavement exists.");
    }

    static void createCircularPavementOfHexagonalCells(int n, double circumRadius, double K, double gamma, double lambda,
        std::vector<std::shared_ptr<Vertex<2>>> &vertices, std::vector<std::shared_ptr<Edge<2>>> &edges, 
        std::vector<std::shared_ptr<Cell<2>>> &cells, 
        std::vector<std::vector<std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>>> &connectivity,
        std::shared_ptr<Dynamics<2>> &dyn) 
    {
        EPC_ASSERT(false && "No generic implementation for circular hexagonal pavement exists.");
    }

    static void createRectangularPavementOfSquareCells(int nx, int ny, 
        double length, double K, double gamma, double lambda, double sdGamma, double sdLambda, 
        std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
        std::vector<std::shared_ptr<Cell<d>>> &cells, 
        std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
        std::shared_ptr<Dynamics<d>> &dyn,unsigned long iT)  
    {
        EPC_ASSERT(false && "No generic implementation for square pavement exists.");
    }

    static void createRingOfTrapezoides(unsigned nCells, double cellHeight, double ringRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, 
        std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
        std::vector<std::shared_ptr<Cell<d>>> &cells, std::shared_ptr<Dynamics<d>> &dyn,unsigned long iT)  
    {
        EPC_ASSERT(false && "No generic implementation for the creation of a ring of rectangles exists.");
    }

    static void createRectangularPavementOfSquareCells(int nx, int ny, 
        double initLength, double K, double A0, double gamma, double lambda,
        std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
        std::vector<std::shared_ptr<Cell<d>>> &cells, 
        std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
        std::shared_ptr<Dynamics<d>> &dyn)  
    {
        EPC_ASSERT(false && "No generic implementation for square pavement exists.");
    }

    static void createCircularPavementOfSquareCells(int n,
        double length, double K, double gamma, double lambda,
        std::vector<std::shared_ptr<Vertex<d>>> &vertices, std::vector<std::shared_ptr<Edge<d>>> &edges, 
        std::vector<std::shared_ptr<Cell<d>>> &cells, 
        std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &connectivity,
        std::shared_ptr<Dynamics<d>> &dyn)  
    {
        EPC_ASSERT(false && "No generic implementation for circular pavement of square cells exists.");
    }
    
};

template <>
struct cellClusterHelpersImpl<2> {
/* 

Exmaple with a 4 x 11 pavement
           _   _   _   _   _     
         _/ \_/ \_/ \_/ \_/ \_          
        / \_/ \_/ \_/ \_/ \_/ \
3:      \_/ \_/ \_/ \_/ \_/ \_/ 
        / \_/ \_/ \_/ \_/ \_/ \
2:      \_/ \_/ \_/ \_/ \_/ \_/ 
        / \_/ \_/ \_/ \_/ \_/ \
1:      \_/ \_/ \_/ \_/ \_/ \_/ 
        / \_/ \_/ \_/ \_/ \_/ \
0:      \_/ \_/ \_/ \_/ \_/ \_/ 

         0 1 2 3 4 5 6 7 8 9 10

*/

/// Creation of rectangular pavement of hexagons (nx * ny hexagons of radius circumRadius)
static void createRectangularPavementOfHexagonalCells(int nx, int ny, double circumRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, 
    std::vector<std::shared_ptr<Vertex<2>>> &vertices, std::vector<std::shared_ptr<Edge<2>>> &edges, 
    std::vector<std::shared_ptr<Cell<2>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>>> &connectivity,
    std::shared_ptr<Dynamics<2>> &dyn, bool isLarge,unsigned long iT) 
{
    using namespace mathConstants;

    vertices.clear();
    edges.clear();
    cells.clear();


    auto incircRadius = circumRadius * sqrt(3.0) / 2.0; // circumscrpit circle radius
    const double A0 = incircRadius*circumRadius/2.0 * 6.0; // area of an hexagon of side circumRadius

    // std::cout << incircRadius << " " << circumRadius << std::endl;

    for (int iX=0; iX<(int)nx; ++iX) {
        double x0 = (double)iX * 1.5 * circumRadius; // x component of the center of the cell
        for (int iY=-(iX % 2)*isLarge; iY<(int)(ny-(iX % 2)*(1-isLarge)); ++iY){ // iX % 2 is for removing the last cells on top of the square.
            double offset = (double)(iX % 2)*incircRadius; // offset for odd numbered cells
            double y0 = offset + 2.0 * incircRadius * (double)iY; // y component of the center of the cell 

            std::vector<std::shared_ptr<Vertex<2>>> localVertices; // tempoarary container vertices for (iX,iY) cell creation
            for (auto theta : HexagonAngles) { // for all angles in an hexagon ...
                Array<double,2> pos; //define the position of each vertex ...
                pos[0] = x0 + circumRadius * cos(theta);
                pos[1] = y0 + circumRadius * sin(theta); 

                // std::cout << iX << " , " << iY << " " << pos << std::endl;

                bool alreadyExisting = false; 
                for (auto &v : vertices) { // check if it hasn't already be created ...
                    if (arrayOperations<double,2>::isEqual(v->getPosition(),pos,0.1*sqrt(A0))) {
                        alreadyExisting = true;
                        localVertices.push_back(v); // if yes then add the existing vertex to the cell list of vertices ...
                        // std::cout << "========= ALREADY EXISTING ===============" << std::endl;
                        break;
                    }
                }
                if (!alreadyExisting) { // if not add the vertex to the global list and to the cell list of vertices
                    auto v = std::shared_ptr<Vertex<2>>(new Vertex<2>(pos,dyn));
                    vertices.push_back(v);
                    localVertices.push_back(v);
                }
            }
            // create the cell
            cells.push_back(cellHelpers<2>::createCell(localVertices, K, A0, gamma, lambda, sdGamma, sdLambda, iT)); 
        }
    }

    // vector conatining the connectivity between cells and the connecting edges
    connectivity.resize(cells.size());
    for (unsigned iC = 0; iC < cells.size(); ++iC) { // for each pair of cells
        for (unsigned iC2 = iC+1; iC2 < cells.size(); ++iC2) { 
            // if edge is duplicated (means that there is a common edge and therefore the cells are neighrbors)
            if (cellHelpers<2>::removeDuplicatedEdge(cells[iC],cells[iC2])) { 
                auto ce = cellHelpers<2>::getCommonEdge(cells[iC],cells[iC2]);
                // add cell and edge to connectivity
                connectivity[iC].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC2],ce)); 
                connectivity[iC2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC],ce));
            }
        }
    }

    // sanity check (max neighbors is 6 for hexagons). on boundaries can be less
    for (auto c: connectivity) {
        EPC_ASSERT(c.size() <= 6 && "max number of neighbors is 6. error here there are more.");
    }
    // add edges in the edge container
    for (auto c: cells) {
        for (auto &e: c->getEdges()) {
            edgeHelpers<2>::appendIfNotInVector(e, edges);
        }
    }
}

/// Creation of rectangular pavement of hexagons (nx * ny hexagons of radius circumRadius)
static void createRectangularPavementOfHexagonalCells(int nx, int ny, double initialCircumRadius, double K, double A0, double gamma, double lambda,
    std::vector<std::shared_ptr<Vertex<2>>> &vertices, std::vector<std::shared_ptr<Edge<2>>> &edges, 
    std::vector<std::shared_ptr<Cell<2>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>>> &connectivity,
    std::shared_ptr<Dynamics<2>> &dyn, bool isLarge) 
{
    using namespace mathConstants;

    vertices.clear();
    edges.clear();
    cells.clear();

    auto incircRadius = initialCircumRadius * sqrt(3.0) / 2.0; // circumscrpit circle radius
    //const double initA = incircRadius*initialCircumRadius/2.0 * 6.0; // area of an hexagon of side circumRadius

    for (int iX=0; iX<(int)nx; ++iX) {
        double x0 = iX * 1.5 * initialCircumRadius; // x component of the center of the cell
        for (int iY=-(iX % 2)*isLarge; iY<(int)(ny-(iX % 2)*(1-isLarge)); ++iY){ // iX % 2 is for removing the last cells on top of the square.
            double offset = (double)(iX % 2)*incircRadius; // offset for odd numbered cells
            double y0 = offset + 2.0 * incircRadius * iY; // y component of the center of the cell 

            std::vector<std::shared_ptr<Vertex<2>>> localVertices; // tempoarary container vertices for (iX,iY) cell creation
            for (auto theta : HexagonAngles) { // for all angles in an hexagon ...
                Array<double,2> pos; //define the position of each vertex ...
                pos[0] = x0 + initialCircumRadius * cos(theta);
                pos[1] = y0 + initialCircumRadius * sin(theta); 

                bool alreadyExisting = false; 
                for (auto &v : vertices) { // check if it hasn't already be created ...
                    if (arrayOperations<double,2>::isEqual(v->getPosition(),pos,0.1*sqrt(A0))) {
                        alreadyExisting = true;
                        localVertices.push_back(v); // if yes then add the existing vertex to the cell list of vertices ...
                        break;
                    }
                }
                if (!alreadyExisting) { // if not add the vertex to the global list and to the cell list of vertices
                    auto v = std::shared_ptr<Vertex<2>>(new Vertex<2>(pos,dyn));
                    vertices.push_back(v);
                    localVertices.push_back(v);
                }
            }
            // create the cell
            cells.push_back(cellHelpers<2>::createCell(localVertices, K, A0, gamma, lambda)); 
        }
    }

    // vector conatining the connectivity between cells and the connecting edges
    connectivity.resize(cells.size());
    for (unsigned iC = 0; iC < cells.size(); ++iC) { // for each pair of cells
        for (unsigned iC2 = iC+1; iC2 < cells.size(); ++iC2) { 
            // if edge is duplicated (means that there is a common edge and therefore the cells are neighrbors)
            if (cellHelpers<2>::removeDuplicatedEdge(cells[iC],cells[iC2])) { 
                auto ce = cellHelpers<2>::getCommonEdge(cells[iC],cells[iC2]);
                // add cell and edge to connectivity
                connectivity[iC].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC2],ce)); 
                connectivity[iC2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC],ce));
            }
        }
    }

    // sanity check (max neighbors is 6 for hexagons). on boundaries can be less
    for (auto c: connectivity) {
        EPC_ASSERT(c.size() <= 6 && "max number of neighbors is 6. error here there are more.");
    }
    // add edges in the edge container
    for (auto c: cells) {
        for (auto &e: c->getEdges()) {
            edgeHelpers<2>::appendIfNotInVector(e, edges);
        }
    }
}

/* 

Exmaple with a 4 x 11 pavement
           _   _   _   _   _     
         _/ \_/ \_/ \_/ \_/ \_          
        / \_/ \_/ \_/ \_/ \_/ \
3:      \_/ \_/ \_/ \_/ \_/ \_/ 
        / \_/ \_/ \_/ \_/ \_/ \
2:      \_/ \_/ \_/ \_/ \_/ \_/ 
        / \_/ \_/ \_/ \_/ \_/ \
1:      \_/ \_/ \_/ \_/ \_/ \_/ 
        / \_/ \_/ \_/ \_/ \_/ \
0:      \_/ \_/ \_/ \_/ \_/ \_/ 

         0 1 2 3 4 5 6 7 8 9 10

*/

/// Creation of rectangular pavement of hexagons (nx * ny hexagons of radius circumRadius)
static void createCircularPavementOfHexagonalCells(int n, double circumRadius, double K, double gamma, double lambda,
    std::vector<std::shared_ptr<Vertex<2>>> &vertices, std::vector<std::shared_ptr<Edge<2>>> &edges, 
    std::vector<std::shared_ptr<Cell<2>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>>> &connectivity,
    std::shared_ptr<Dynamics<2>> &dyn) 
{
    using namespace mathConstants;

    vertices.clear();
    edges.clear();
    cells.clear();

    auto incircRadius = circumRadius * sqrt(3.0) / 2.0; // circumscrpit circle radius
    const double A0 = incircRadius*circumRadius/2.0 * 6.0; // area of an hexagon of side circumRadius

    Array<double,2> center; 
    center[0] = n/2*1.5*circumRadius;
    center[1] = (double)((n/2) % 2)*incircRadius+2.0 * incircRadius * n/2;
    double totRadius = 1.5 * incircRadius * (n/2);


    for (int iX=0; iX<n; ++iX) {
        double x0 = iX * 1.5 * circumRadius; // x component of the center of the cell
        for (int iY=-(iX % 2); iY<n; ++iY){ // iX % 2 is for removing the last cells on top of the square.
            double offset = (double)(iX % 2)*incircRadius; // offset for odd numbered cells
            double y0 = offset + 2.0 * incircRadius * iY; // y component of the center of the cell 

            auto rad = (Array<double,2>(x0,y0)-center).norm();

            if (rad < totRadius) {

                std::vector<std::shared_ptr<Vertex<2>>> localVertices; // tempoarary container vertices for (iX,iY) cell creation
                for (auto theta : HexagonAngles) { // for all angles in an hexagon ...
                    Array<double,2> pos; //define the position of each vertex ...
                    pos[0] = x0 + circumRadius * cos(theta);
                    pos[1] = y0 + circumRadius * sin(theta); 

                    bool alreadyExisting = false; 
                    for (auto &v : vertices) { // check if it hasn't already be created ...
                        if (arrayOperations<double,2>::isEqual(v->getPosition(),pos,std::numeric_limits<double>::epsilon()*100.0*A0*n*n)) {
                            alreadyExisting = true;
                            localVertices.push_back(v); // if yes then add the existing vertex to the cell list of vertices ...
                            break;
                        }
                    }
                    if (!alreadyExisting) { // if not add the vertex to the global list and to the cell list of vertices
                        auto v = std::shared_ptr<Vertex<2>>(new Vertex<2>(pos,dyn));
                        vertices.push_back(v);
                        localVertices.push_back(v);
                    }
                }
                // create the cell
                cells.push_back(cellHelpers<2>::createCell(localVertices, K, A0, gamma, lambda)); 
            }
        }
    }

    // vector conatining the connectivity between cells and the connecting edges
    connectivity.resize(cells.size());
    for (unsigned iC = 0; iC < cells.size(); ++iC) { // for each pair of cells
        for (unsigned iC2 = iC+1; iC2 < cells.size(); ++iC2) { 
            // if edge is duplicated (means that there is a common edge and therefore the cells are neighrbors)
            if (cellHelpers<2>::removeDuplicatedEdge(cells[iC],cells[iC2])) { 
                auto ce = cellHelpers<2>::getCommonEdge(cells[iC],cells[iC2]);
                // add cell and edge to connectivity
                connectivity[iC].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC2],ce)); 
                connectivity[iC2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC],ce));
            }
        }
    }

    // sanity check (max neighbors is 6 for hexagons). on boundaries can be less
    for (auto c: connectivity) {
        EPC_ASSERT(c.size() <= 6 && "max number of neighbors is 6. error here there are more.");
    }
    // add edges in the edge container
    for (auto c: cells) {
        for (auto &e: c->getEdges()) {
            edgeHelpers<2>::appendIfNotInVector(e, edges);
        }
    }
}


/* 

Exmaple with a 4 x 11 pavement
         _ _ _ _ _ _ _ _ _ _ _        
6:      |_|_|_|_|_|_|_|_|_|_|_| 
5:      |_|_|_|_|_|_|_|_|_|_|_| 
4:      |_|_|_|_|_|_|_|_|_|_|_| 
3:      |_|_|_|_|_|_|_|_|_|_|_| 
2:      |_|_|_|_|_|_|_|_|_|_|_| 
1:      |_|_|_|_|_|_|_|_|_|_|_| 
0:      |_|_|_|_|_|_|_|_|_|_|_| 

         0 1 2 3 4 5 6 7 8 9 10

*/

/// Creation of rectangular pavement of hexagons (nx * ny hexagons of radius circumRadius)
static void createRectangularPavementOfSquareCells(int nx, int ny, double length, double K, double gamma, double lambda, double sdGamma, double sdLambda, 
    std::vector<std::shared_ptr<Vertex<2>>> &vertices, std::vector<std::shared_ptr<Edge<2>>> &edges, 
    std::vector<std::shared_ptr<Cell<2>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>>> &connectivity,
    std::shared_ptr<Dynamics<2>> &dyn,unsigned long iT) 
{
    using namespace mathConstants;

    vertices.clear();
    edges.clear();
    cells.clear();

    const auto halfDiag = length/sqrt(2.0);
    const auto A0 = length*length;

    for (int iX=0; iX<nx; iX++) {
        double x0 = iX * length; // x component of the center of the cell
        for (int iY=0; iY<ny; iY++){ // iX % 2 is for removing the last cells on top of the square.
            double y0 = length * iY; // y component of the center of the cell 

            std::vector<std::shared_ptr<Vertex<2>>> localVertices; // tempoarary container vertices for (iX,iY) cell creation
            for (auto theta : SquareAngles) { // for all angles in an hexagon ...
                Array<double,2> pos; //define the position of each vertex ...
                pos[0] = x0 + halfDiag * cos(theta);
                pos[1] = y0 + halfDiag * sin(theta); 

                // std::cout << pos << std::endl;

                bool alreadyExisting = false; 
                for (auto &v : vertices) { // check if it hasn't already be created ...
                    if (arrayOperations<double,2>::isEqual(v->getPosition(),pos,std::numeric_limits<double>::epsilon()*100.0*A0*nx*ny)) {
                        alreadyExisting = true;
                        localVertices.push_back(v); // if yes then add the existing vertex to the cell list of vertices ...
                        break;
                    }
                }
                if (!alreadyExisting) { // if not add the vertex to the global list and to the cell list of vertices
                    auto v = std::shared_ptr<Vertex<2>>(new Vertex<2>(pos,dyn));
                    vertices.push_back(v);
                    localVertices.push_back(v);
                }
            }
            // create the cell
            cells.push_back(cellHelpers<2>::createCell(localVertices, K, A0, gamma, lambda)); 
        }
    }

    // vector conatining the connectivity between cells and the connecting edges
    connectivity.resize(cells.size());
    for (unsigned iC = 0; iC < cells.size(); ++iC) { // for each pair of cells
        for (unsigned iC2 = iC+1; iC2 < cells.size(); ++iC2) { 
            // if edge is duplicated (means that there is a common edge and therefore the cells are neighrbors)
            if (cellHelpers<2>::removeDuplicatedEdge(cells[iC],cells[iC2])) { 
                auto ce = cellHelpers<2>::getCommonEdge(cells[iC],cells[iC2]);
                // add cell and edge to connectivity
                connectivity[iC].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC2],ce)); 
                connectivity[iC2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC],ce));
            }
        }
    }

    // sanity check (max neighbors is 6 for hexagons). on boundaries can be less
    for (auto c: connectivity) {
        EPC_ASSERT(c.size() <= 4 && "max number of neighbors is 4. error here there are more.");
    }
    // add edges in the edge container
    for (auto c: cells) {
        for (auto &e: c->getEdges()) {
            edgeHelpers<2>::appendIfNotInVector(e, edges);
        }
    }
    // exit(-1);
}

/// Creation of rectangular pavement of hexagons (nx * ny hexagons of radius circumRadius)
static void createRectangularPavementOfSquareCells(int nx, int ny, double initLength, double K, double A0, double gamma, double lambda,
    std::vector<std::shared_ptr<Vertex<2>>> &vertices, std::vector<std::shared_ptr<Edge<2>>> &edges, 
    std::vector<std::shared_ptr<Cell<2>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>>> &connectivity,
    std::shared_ptr<Dynamics<2>> &dyn) 
{
    using namespace mathConstants;

    vertices.clear();
    edges.clear();
    cells.clear();

    const auto halfDiag = initLength/sqrt(2.0);
    //const auto initA = initLength*initLength;

    for (int iX=0; iX<nx; iX++) {
        double x0 = iX * initLength; // x component of the center of the cell
        for (int iY=0; iY<ny; iY++){ // iX % 2 is for removing the last cells on top of the square.
            double y0 = initLength * iY; // y component of the center of the cell 

            std::vector<std::shared_ptr<Vertex<2>>> localVertices; // tempoarary container vertices for (iX,iY) cell creation
            for (auto theta : SquareAngles) { // for all angles in an hexagon ...
                Array<double,2> pos; //define the position of each vertex ...
                pos[0] = x0 + halfDiag * cos(theta);
                pos[1] = y0 + halfDiag * sin(theta); 

                // std::cout << pos << std::endl;

                bool alreadyExisting = false; 
                for (auto &v : vertices) { // check if it hasn't already be created ...
                    if (arrayOperations<double,2>::isEqual(v->getPosition(),pos,std::numeric_limits<double>::epsilon()*100.0*A0*nx*ny)) {
                        alreadyExisting = true;
                        localVertices.push_back(v); // if yes then add the existing vertex to the cell list of vertices ...
                        break;
                    }
                }
                if (!alreadyExisting) { // if not add the vertex to the global list and to the cell list of vertices
                    auto v = std::shared_ptr<Vertex<2>>(new Vertex<2>(pos,dyn));
                    vertices.push_back(v);
                    localVertices.push_back(v);
                }
            }
            // create the cell
            cells.push_back(cellHelpers<2>::createCell(localVertices, K, A0, gamma, lambda)); 
        }
    }

    // vector conatining the connectivity between cells and the connecting edges
    connectivity.resize(cells.size());
    for (unsigned iC = 0; iC < cells.size(); ++iC) { // for each pair of cells
        for (unsigned iC2 = iC+1; iC2 < cells.size(); ++iC2) { 
            // if edge is duplicated (means that there is a common edge and therefore the cells are neighrbors)
            if (cellHelpers<2>::removeDuplicatedEdge(cells[iC],cells[iC2])) { 
                auto ce = cellHelpers<2>::getCommonEdge(cells[iC],cells[iC2]);
                // add cell and edge to connectivity
                connectivity[iC].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC2],ce)); 
                connectivity[iC2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC],ce));
            }
        }
    }

    // sanity check (max neighbors is 6 for hexagons). on boundaries can be less
    for (auto c: connectivity) {
        EPC_ASSERT(c.size() <= 4 && "max number of neighbors is 4. error here there are more.");
    }
    // add edges in the edge container
    for (auto c: cells) {
        for (auto &e: c->getEdges()) {
            edgeHelpers<2>::appendIfNotInVector(e, edges);
        }
    }
    // exit(-1);
}


/* 

Exmaple with a 4 x 11 pavement
         _ _ _ _ _ _ _ _ _ _ _        
6:      |_|_|_|_|_|_|_|_|_|_|_| 
5:      |_|_|_|_|_|_|_|_|_|_|_| 
4:      |_|_|_|_|_|_|_|_|_|_|_| 
3:      |_|_|_|_|_|_|_|_|_|_|_| 
2:      |_|_|_|_|_|_|_|_|_|_|_| 
1:      |_|_|_|_|_|_|_|_|_|_|_| 
0:      |_|_|_|_|_|_|_|_|_|_|_| 

         0 1 2 3 4 5 6 7 8 9 10

*/

/// Creation of rectangular pavement of hexagons (nx * ny hexagons of radius circumRadius)
static void createCircularPavementOfSquareCells(int n, double length, double K, double gamma, double lambda,
    std::vector<std::shared_ptr<Vertex<2>>> &vertices, std::vector<std::shared_ptr<Edge<2>>> &edges, 
    std::vector<std::shared_ptr<Cell<2>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>>> &connectivity,
    std::shared_ptr<Dynamics<2>> &dyn) 
{
    using namespace mathConstants;

    vertices.clear();
    edges.clear();
    cells.clear();

    const auto halfDiag = length/sqrt(2.0);
    const auto A0 = length*length;

    Array<double,2> center; 
    center[0] = n/2*length;
    center[1] = n/2*length;
    double totRadius = n/2*length;

    for (int iX=-1; iX<n+1; iX++) {
        double x0 = iX * length; // x component of the center of the cell
        for (int iY=-1; iY<n+1; iY++){ // iX % 2 is for removing the last cells on top of the square.
            double y0 = length * iY; // y component of the center of the cell 

            auto rad = (Array<double,2>(x0,y0)-center).norm();

            if (rad < totRadius) {

                std::vector<std::shared_ptr<Vertex<2>>> localVertices; // tempoarary container vertices for (iX,iY) cell creation
                for (auto theta : SquareAngles) { // for all angles in an hexagon ...
                    Array<double,2> pos; //define the position of each vertex ...
                    pos[0] = x0 + halfDiag * cos(theta);
                    pos[1] = y0 + halfDiag * sin(theta); 

                    // std::cout << pos << std::endl;

                    bool alreadyExisting = false; 
                    for (auto &v : vertices) { // check if it hasn't already be created ...
                        if (arrayOperations<double,2>::isEqual(v->getPosition(),pos,std::numeric_limits<double>::epsilon()*100.0*A0*n*n)) {
                            alreadyExisting = true;
                            localVertices.push_back(v); // if yes then add the existing vertex to the cell list of vertices ...
                            break;
                        }
                    }
                    if (!alreadyExisting) { // if not add the vertex to the global list and to the cell list of vertices
                        auto v = std::shared_ptr<Vertex<2>>(new Vertex<2>(pos,dyn));
                        vertices.push_back(v);
                        localVertices.push_back(v);
                    }
                }
                // create the cell
                cells.push_back(cellHelpers<2>::createCell(localVertices, K, A0, gamma, lambda)); 
            }
        }
    }

    // vector conatining the connectivity between cells and the connecting edges
    connectivity.resize(cells.size());
    for (unsigned iC = 0; iC < cells.size(); ++iC) { // for each pair of cells
        for (unsigned iC2 = iC+1; iC2 < cells.size(); ++iC2) { 
            // if edge is duplicated (means that there is a common edge and therefore the cells are neighrbors)
            if (cellHelpers<2>::removeDuplicatedEdge(cells[iC],cells[iC2])) { 
                auto ce = cellHelpers<2>::getCommonEdge(cells[iC],cells[iC2]);
                // add cell and edge to connectivity
                connectivity[iC].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC2],ce)); 
                connectivity[iC2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC],ce));
            }
        }
    }

    // sanity check (max neighbors is 6 for hexagons). on boundaries can be less
    for (auto c: connectivity) {
        EPC_ASSERT(c.size() <= 4 && "max number of neighbors is 4. error here there are more.");
    }
    // add edges in the edge container
    for (auto c: cells) {
        for (auto &e: c->getEdges()) {
            edgeHelpers<2>::appendIfNotInVector(e, edges);
        }
    }
    // exit(-1);
}


/// Creation of ring of trapezoides (nCells trapezoides of height cellHeight in a ring of radius ringRadius)
static void createRingOfTrapezoides(int nCells, double cellHeight, double ringRadius, double K, double gamma, double lambda,
	double sdGamma, double sdLambda, 
    std::vector<std::shared_ptr<Vertex<2>>> &vertices, std::vector<std::shared_ptr<Edge<2>>> &edges, 
    std::vector<std::shared_ptr<Cell<2>>> &cells, 
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>>> &connectivity,
    std::shared_ptr<Dynamics<2>> &dyn,unsigned long iT) 
{
    using namespace mathConstants;

    vertices.clear();
    edges.clear();
    cells.clear();

    double trapezoidHeight = cos(pi/nCells)*cellHeight;
    double trapezoideBigBase = 2*sin(pi/nCells)*ringRadius;
    double trapezoideSmallBase = 2*sin(pi/nCells)*(ringRadius-cellHeight);
    double A0 = 0.5*trapezoidHeight*(trapezoideBigBase+trapezoideSmallBase);

    for (int iC=0; iC<nCells; iC++) {
        std::vector<std::shared_ptr<Vertex<2>>> localVertices; // tempoarary container vertices for (iX,iY) cell creation
        for (int iV=0; iV<4; iV++) {
            double r;
            double theta;
            if (iV < 2){
                r = ringRadius;
            }
            else {
                r = ringRadius - cellHeight;    
            }
            if ((iV == 0) or (iV == 3)) {
                theta = iC*2*pi/nCells;
            }
            else{
                theta = (1+iC)*2*pi/nCells;     
            }

            Array<double,2> pos; //define the position of each vertex ...
            pos[0] = r * cos(theta);
            pos[1] = r * sin(theta); 

            bool alreadyExisting = false; 
            for (auto &v : vertices) { // check if it hasn't already be created ...
                if (arrayOperations<double,2>::isEqual(v->getPosition(),pos,std::numeric_limits<double>::epsilon()*100.0*A0*nCells)) {
                    alreadyExisting = true;
                    localVertices.push_back(v); // if yes then add the existing vertex to the cell list of vertices ...
                    break;
                }
            }
            if (!alreadyExisting) { // if not add the vertex to the global list and to the cell list of vertices
                auto v = std::shared_ptr<Vertex<2>>(new Vertex<2>(pos,dyn));
                vertices.push_back(v);
                localVertices.push_back(v);
            }
        }
        
        // create the cell
        cells.push_back(cellHelpers<2>::createCell(localVertices, K, A0, gamma, lambda, sdGamma, sdLambda, iT)); 

    }

    // vector containing the connectivity between cells and the connecting edges
    connectivity.resize(cells.size());
    for (unsigned iC = 0; iC < cells.size(); ++iC) { // for each pair of cells
        for (unsigned iC2 = iC+1; iC2 < cells.size(); ++iC2) { 
            // if edge is duplicated (means that there is a common edge and therefore the cells are neighrbors)
            if (cellHelpers<2>::removeDuplicatedEdge(cells[iC],cells[iC2])) { 

                auto ce = cellHelpers<2>::getCommonEdge(cells[iC],cells[iC2]);
                // add cell and edge to connectivity
                connectivity[iC].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC2],ce)); 
                connectivity[iC2].push_back(std::tuple<std::shared_ptr<Cell<2>>,std::shared_ptr<Edge<2>>>(cells[iC],ce));
            }
        }
    }

    // sanity check (max neighbors is 6 for hexagons). on boundaries can be less
    for (auto c: connectivity) {
        EPC_ASSERT(c.size() <= 6 && "max number of neighbors is 6. error here there are more.");
    }
    // add edges in the edge container
    for (auto c: cells) {
        for (auto &e: c->getEdges()) {
            edgeHelpers<2>::appendIfNotInVector(e, edges);
        }
    }
}

}; // end cellCluster helper impl

}  // end namespace epc

#endif
