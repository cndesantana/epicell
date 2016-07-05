#ifndef CELL_CLUSTER_HH
#define CELL_CLUSTER_HH

#include "vertex.h"
#include "edge.h"
#include "cellClusterHelpers.h"
#include "cellCluster.h"


namespace epc {

// =============== Constuctors ========================= //
template<int d>
CellCluster<d>::CellCluster()
{
    cells.clear();
    edges.clear();
    vertices.clear();
    resetForce();
}

// ========== Interface functions ================ //

template<int d>
std::vector<std::shared_ptr<Cell<d>>> const& CellCluster<d>::getCells() const {
    return cells;
}

template<int d>
std::vector<std::shared_ptr<Cell<d>>> & CellCluster<d>::getCells() {
    return cells;
}

template<int d>
std::vector<std::shared_ptr<Cell<d>>> CellCluster<d>::getCellsOfType(std::shared_ptr<CellType> type_) {
    std::vector<std::shared_ptr<Cell<d>>> cellsOfGivenType;
    for (auto c: cells){
        if (c->getType().getId()==type_->getId()){
            cellsOfGivenType.push_back(c);
        }
    }
    return cellsOfGivenType;
}

template<int d>
std::vector<std::shared_ptr<Edge<d>>> const& CellCluster<d>::getEdges() const {
    return edges;
}

template<int d>
std::vector<std::shared_ptr<Edge<d>>> & CellCluster<d>::getEdges() {
    return edges;
}

template<int d>
std::vector<std::shared_ptr<Vertex<d>>> const& CellCluster<d>::getVertices() const {
    return vertices;
}

template<int d>
std::vector<std::shared_ptr<Vertex<d>>> & CellCluster<d>::getVertices() {
    return vertices;
}

template<int d>
std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> const& CellCluster<d>::getConnectivity() const {
    return connectivity;
}

template<int d>
std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> & CellCluster<d>::getConnectivity() {
    return connectivity;
}


template<int d>
CellTypesManager & CellCluster<d>::getCellTypesManager() {
    return cellTypesManager;
}

template<int d>
double const& CellCluster<d>::getBoundaryLambda() const {
    return boundaryLambda;
}

template<int d>
std::vector<std::shared_ptr<Vertex<d>>> const& CellCluster<d>::getBoundaryVertices() const {
    return boundaryVertices;
}

// ========== Line tensions ================ //
// Set boundary line tension.
template<int d>
void CellCluster<d>::setBoundaryLambda(    double boundaryLambda_){
    boundaryLambda = boundaryLambda_;
    // Set the lambda of boundary edges.
    for (auto boundaryEdge : epitheliumHelpers<d>::getAllBoundaryEdges(*this)){
        boundaryEdge->setLambda(boundaryLambda_);//standard deviation is '0'
    }

}

// Set the list of boundary vertices.
template<int d>
void CellCluster<d>::setBoundaryVertices(){
    boundaryVertices = epitheliumHelpers<d>::getAllBoundaryVertices(*this);
}


// Add a vertex to the list of boundary vertices.
template<int d>
void CellCluster<d>::addToBoundaryVertices(const std::shared_ptr<Vertex<d>> &v){
    vertexHelpers<d>::appendIfNotInVector(v, boundaryVertices);
}

// Remove a vertex from the list of boundary vertices.
template<int d>
void CellCluster<d>::removeFromBoundaryVertices(const std::shared_ptr<Vertex<d>> &v){
    for (unsigned bv = 0; bv < boundaryVertices.size(); bv++)
    {
        if (boundaryVertices[bv]->getId() == v->getId()){
            boundaryVertices.erase(boundaryVertices.begin()+bv);
            break;
        } 
    }
    EPC_ASSERT("Asked to remove a vertex from the list of boundary vertices. BUT it was not found in this list.");
}

// ========== geometric functions ================ //
template<int d>
Array<double,d> CellCluster<d>::computeCentroid() const {
    auto centroid = Array<double,d>(double());
    for ( auto v : vertices ) {
        centroid += v->getPosition();
    }
    centroid /= (double)vertices.size();
    return centroid;
}

template<int d>
double CellCluster<d>::computeArea() const {
    double area = double();
    for ( auto cell : cells ) {
        area += cell->computeArea();
    }

    return area;
}

// ========== dynamic functions ================ //

template<int d>
void CellCluster<d>::advance(double dt) {
    for (auto vertex : vertices) {
        vertex->advance(dt);
    }
}

template<int d>
void CellCluster<d>::resetForce() {
    Array<double,d> zero(0.0);
    for (auto vertex : vertices) {
        vertex->setForce(zero);
    }
}

template<int d>
void CellCluster<d>::updateForce() {
    for (auto cell : cells) {
        cell->computeAndUpdateForce();
    }
    for (auto edge : edges) {
        edge->computeAndUpdateForce();
    }
}

template<int d>
void CellCluster<d>::remCell(std::shared_ptr<Cell<d>> &mycell){
//    std::cout << "Removing cell c" << mycell->getId() << std::endl;
    for(int i=0;i<(int)cells.size();i++){
        if (cells[i]->getId() == mycell->getId()){
                cells.erase(cells.begin()+i);
                return;
        }
    }
    return;
}    

/// Remove a vertex from the CellCluster
template<int d>
void CellCluster<d>::removeVertex(std::shared_ptr<Vertex<d>> &myvertex){
    for(int i=0;i<(int)vertices.size();i++){
        if (vertices[i]->getId() == myvertex->getId()){
                vertices.erase(vertices.begin()+i);
                break;
        }
    }
}    


/// Remove an edge from the CellCluster
template<int d>
void CellCluster<d>::removeEdge(std::shared_ptr<Edge<d>> &myedge){
    for(int i=0;i<(int)edges.size();i++){
        if (edges[i]->getId() == myedge->getId()){
                edges.erase(edges.begin()+i);
                break;
        }
    }
}    

// Add the connection between two cells through the given edge.
template<int d>
void CellCluster<d>::addConnectivityBetweenCellsThroughEdge(    std::shared_ptr<Cell<d>> c1,
                                                                std::shared_ptr<Cell<d>> c2, 
                                                                std::shared_ptr<Edge<d>> e){
    //find the position of c1 in the list of cells of the cell-cluster, and add it a connection to c2 through e.   
    unsigned c1Pos = (std::find(cells.begin(), cells.end(),c1) - cells.begin() ); 
    connectivity[c1Pos].push_back(std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>(c2,e));
    
    //find the position of c2 in the list of cells of the cell-cluster, and add it a connection to c1 through e.
    unsigned c2Pos = (std::find(cells.begin(), cells.end(),c2) - cells.begin() ); 
    connectivity[c2Pos].push_back(std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>(c1,e)); 
    
}

// Remove the connection between two cells through the given edge.
template<int d>
bool CellCluster<d>::removeConnectivityBetweenCellsThroughEdge( std::shared_ptr<Cell<d>> c1,
                                                                std::shared_ptr<Cell<d>> c2, 
                                                                std::shared_ptr<Edge<d>> e){
    //find the position of c1 and c2 in the list of cells of the cell-cluster
    unsigned c1Pos = (std::find(cells.begin(), cells.end(),c1) - cells.begin() ); 
    unsigned c2Pos = (std::find(cells.begin(), cells.end(),c2) - cells.begin() ); 

    bool foundConnection = false;
    for (int iB=0; iB < connectivity[c1Pos].size(); iB++){
        if (std::get<0>(connectivity[c1Pos][iB])->getId() == c2->getId()){
            if (std::get<1>(connectivity[c1Pos][iB])->getId() == e->getId()){
                connectivity[c1Pos].erase(connectivity[c1Pos].begin()+iB);
                foundConnection = true;
                break;
            } else{
                std::cout << "ERROR: c1 -> c2 are connected but through a different edge than the one given!" << std::endl;
            }
        }
    }

    if (foundConnection == false){
        std::cout << "ERROR: c1 -> c2 are not connected!" << std::endl;
    }

    foundConnection = false;
    for (int iB=0; iB < connectivity[c2Pos].size(); iB++){
        if (std::get<0>(connectivity[c2Pos][iB])->getId() == c1->getId()){
            if (std::get<1>(connectivity[c2Pos][iB])->getId() == e->getId()){
                connectivity[c2Pos].erase(connectivity[c2Pos].begin()+iB);
                foundConnection = true;
                break;
            } else{
                std::cout << "ERROR: c2 -> c1 are connected but through a different edge than the one given!" << std::endl;
            }
        }
    }

    if (foundConnection == false){
        std::cout << "ERROR: c2 -> c1 are not connected!" << std::endl;
    }
    return foundConnection;
}


/// Remove all connectivities with mycell from the cellCluster
template<int d>
void CellCluster<d>::removeConnectivity(std::shared_ptr<Cell<d>> &mycell){
//    std::cout << "BEGIN : Remove connectivity with cell " << mycell->getId() << std::endl;

    // vector conatining the connectivity between cells and the connecting edges
    connectivity.resize(cells.size());
    unsigned mycellPos = (std::find(this->getCells().begin(), this->getCells().end(),mycell) - this->getCells().begin() ); //find the position of mycell in the list of cells of the cell-cluster.
    unsigned numberOfNeighbors = connectivity[mycellPos].size();

//    std::cout << "Initial number of neighbors = " << numberOfNeighbors << std::endl;
//    std::cout << "Initial number of edges = " << mycell->getEdges().size() << std::endl;

    while (numberOfNeighbors > 0) { // for each pair of cells      
        connectivity[mycellPos].erase(connectivity[mycellPos].begin());
        numberOfNeighbors = connectivity[mycellPos].size();
//        std::cout << "Intermediate Number of neighbors = " << numberOfNeighbors << std::endl;

    }
//    std::cout << "END : Remove connectivity with cell " << mycell->getId() << std::endl;

}    

/// Replace all connectivities with mycell by connectivity with cell1 or cell2 in the cellCluster.
template<int d>
void CellCluster<d>::replaceConnectivity(std::shared_ptr<Cell<d>> &mycell, std::shared_ptr<Cell<d>> &cell1, std::shared_ptr<Cell<d>> &cell2){

    
}    

/// Every cell which preferred area A0 increases is included in a list of "Cells on Mitosis"    
template<int d>
void CellCluster<d>::addCellOnMitosis(std::shared_ptr<Cell<d>> &cell){
	if(cellsOnMitosis.size() > 0){
		auto posCell = std::find(std::begin(cellsOnMitosis), std::end(cellsOnMitosis), cell);
		if (posCell == std::end(cellsOnMitosis)){//if cell is not in the list of cells on Mitosis yet
			cell->setOnMitosis(true);
            cell->setAreaBeforeMitosis(cell->computeArea());
			cellsOnMitosis.push_back(cell);
		}        
	}
	else{
        cell->setOnMitosis(true);
        cell->setAreaBeforeMitosis(cell->computeArea());
		cellsOnMitosis.push_back(cell);
//		initialAreaCOM.push_back(cell->computeArea());
	}
}

template<int d>
std::vector<std::shared_ptr<Cell<d>>> CellCluster<d>::getCellsToDivide(){
    std::vector<std::shared_ptr<Cell<d>>> cellsToDivide;

    for(int i=0;i<int(cellsOnMitosis.size());i++){
      if(cellsOnMitosis[i]->computeArea() >= (2.0*cellsOnMitosis[i]->getAreaBeforeMitosis())){
            cellsToDivide.push_back(cellsOnMitosis[i]);
            cellsOnMitosis.erase(cellsOnMitosis.begin()+i);//removing the cell-to-divide from the list of cells on mitosis
    	}
    }
    return(cellsToDivide);
}

// ======================================================================== //
// ========================HexagonalCellsRectangularCluster============================ //
// ======================================================================== //

template<int d>
HexagonalCellsRectangularCluster<d>::HexagonalCellsRectangularCluster(
    int nx, int ny, double circumRadius, double K, double gamma, 
    double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &&dyn, bool isLarge,unsigned long iT) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createRectangularPavementOfHexagonalCells(nx, ny, circumRadius, K, gamma, lambda, sdGamma, sdLambda,  
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn, isLarge,iT);
}

template<int d>
HexagonalCellsRectangularCluster<d>::HexagonalCellsRectangularCluster(
    int nx, int ny, double circumRadius, double K, double gamma, 
    double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &dyn, bool isLarge,unsigned long iT) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createRectangularPavementOfHexagonalCells(nx, ny, circumRadius, K, gamma, lambda, sdGamma, sdLambda,  
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn, isLarge,iT);
}

template<int d>
HexagonalCellsRectangularCluster<d>::HexagonalCellsRectangularCluster(
    int nx, int ny, double initialCircumRadius, double K, double A0, double gamma, 
    double lambda, std::shared_ptr<Dynamics<d>> &&dyn, bool isLarge) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createRectangularPavementOfHexagonalCells(nx, ny, initialCircumRadius, K, A0, gamma, lambda, 
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn, isLarge);
}

template<int d>
HexagonalCellsRectangularCluster<d>::HexagonalCellsRectangularCluster(
    int nx, int ny, double initCircumRadius, double K, double A0, double gamma, 
    double lambda, std::shared_ptr<Dynamics<d>> &dyn, bool isLarge) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createRectangularPavementOfHexagonalCells(nx, ny, initCircumRadius, K, A0, gamma, lambda, 
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn, isLarge);
}

// ======================================================================== //
// ========================SquareCellsRectangularCluster============================ //
// ======================================================================== //

template<int d>
SquareCellsRectangularCluster<d>::SquareCellsRectangularCluster(
    int nx, int ny, double length, double K, double gamma, 
    double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &&dyn,unsigned long iT) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createRectangularPavementOfSquareCells(nx, ny, length, K, gamma, lambda, sdGamma, sdLambda,  
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn,iT);
}

template<int d>
SquareCellsRectangularCluster<d>::SquareCellsRectangularCluster(
    int nx, int ny, double length, double K, double gamma, 
    double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &dyn,unsigned long iT) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createRectangularPavementOfSquareCells(nx, ny, length, K, gamma, lambda, sdGamma, sdLambda,  
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn,iT);
}

template<int d>
SquareCellsRectangularCluster<d>::SquareCellsRectangularCluster(
    int nx, int ny, double initLength, double K, double A0, double gamma, 
    double lambda, std::shared_ptr<Dynamics<d>> &&dyn) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createRectangularPavementOfSquareCells(nx, ny, initLength, K, A0, gamma, lambda, 
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn);
}

template<int d>
SquareCellsRectangularCluster<d>::SquareCellsRectangularCluster(
    int nx, int ny, double initLength, double K, double A0, double gamma, 
    double lambda, std::shared_ptr<Dynamics<d>> &dyn) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createRectangularPavementOfSquareCells(nx, ny, initLength, K, A0, gamma, lambda, 
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn);
}

// ======================================================================== //
// ========================HexagonalCellsCircularCluster============================ //
// ======================================================================== //


template<int d>
HexagonalCellsCircularCluster<d>::HexagonalCellsCircularCluster(
    int n, double circumRadius, double K, double gamma, 
    double lambda, std::shared_ptr<Dynamics<d>> &&dyn) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createCircularPavementOfHexagonalCells(n, circumRadius, K, gamma, lambda, 
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn);
}

template<int d>
HexagonalCellsCircularCluster<d>::HexagonalCellsCircularCluster(
    int n, double circumRadius, double K, double gamma, 
    double lambda, std::shared_ptr<Dynamics<d>> &dyn) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createCircularPavementOfHexagonalCells(n, circumRadius, K, gamma, lambda, 
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn);
}

// ======================================================================== //
// ========================SquareCellsCircularCluster============================ //
// ======================================================================== //

template<int d>
SquareCellsCircularCluster<d>::SquareCellsCircularCluster(
    int n, double length, double K, double gamma, 
    double lambda, std::shared_ptr<Dynamics<d>> &&dyn) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createCircularPavementOfSquareCells(n, length, K, gamma, lambda, 
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn);
}

template<int d>
SquareCellsCircularCluster<d>::SquareCellsCircularCluster(
    int n, double length, double K, double gamma, 
    double lambda, std::shared_ptr<Dynamics<d>> &dyn) : CellCluster<d>()
{ 
    cellClusterHelpers<d>::createCircularPavementOfSquareCells(n, length, K, gamma, lambda, 
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn);
}

// ======================================================================== //
// ========================TrapezoidalCellsRingCluster============================ //
// ======================================================================== //

template<int d>
TrapezoidalCellsRingCluster<d>::TrapezoidalCellsRingCluster(int nCells, double cellHeight, double ringRadius, double K, 
    double gamma, double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &&dyn,unsigned long iT)
{ 
    // Create Ring CellCluster
    cellClusterHelpers<d>::createRingOfTrapezoides(nCells,cellHeight, ringRadius, K, gamma, lambda, sdGamma, sdLambda,
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn, iT);
}

template<int d>
TrapezoidalCellsRingCluster<d>::TrapezoidalCellsRingCluster(int nCells, double cellHeight, double ringRadius, double K, 
    double gamma, double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &dyn,unsigned long iT)
{ 
    // Create Ring CellCluster
    cellClusterHelpers<d>::createRingOfTrapezoides(nCells,cellHeight, ringRadius, K, gamma, lambda, sdGamma, sdLambda,
        this->getVertices(), this->getEdges(), this->getCells(), this->getConnectivity(), dyn, iT);
}



}  // end namespace epc

#endif
