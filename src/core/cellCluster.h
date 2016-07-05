#ifndef CELL_CLUSTER_H
#define CELL_CLUSTER_H

#include "vertex.h"
#include "edge.h"
#include "cellClusterHelpers.h"
#include "cellTypes/cellTypesManager.h"


namespace epc {

template<int d>
class CellCluster {
public:
    /// Default constructor (everything is set to zero)
    CellCluster();
    CellCluster(int nx, int ny, double circumRadius, double K, double gamma, double lambda, 
        double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &&dyn,unsigned long iT);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");
public:
    // ========== Interface functions ================ //

    /// get a non-modifiable reference to cells
    std::vector<std::shared_ptr<Cell<d>>> const &getCells() const;  
    /// get a modifiable reference to cells
    std::vector<std::shared_ptr<Cell<d>>> &getCells();
    /// get a modifiable reference to cells of given type
    std::vector<std::shared_ptr<Cell<d>>> getCellsOfType(std::shared_ptr<CellType> type_);
    /// get a non-modifiable reference to edges
    std::vector<std::shared_ptr<Edge<d>>> const &getEdges() const;
    /// get a modifiable reference to edges
    std::vector<std::shared_ptr<Edge<d>>> &getEdges();
    /// get a non-modifiable reference to vertices
    std::vector<std::shared_ptr<Vertex<d>>> const &getVertices() const;
    /// get a modifiable reference to vertices
    std::vector<std::shared_ptr<Vertex<d>>> &getVertices();
    /// get a non-modifiable reference to connectivity
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> const &getConnectivity() const;
    /// get a modifiable reference to connectivity
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> &getConnectivity();

    /// get a modifiable reference to connectivity
    CellTypesManager &getCellTypesManager();
    /// Get a non-modifiable reference to boubdaryLambda.
    double const& getBoundaryLambda() const;
    // Get a non-modifiable reference to boubdaryLambda.
    std::vector<std::shared_ptr<Vertex<d>>> const &getBoundaryVertices() const;

    /// Every cell which preferred area A0 increases is included in a list of "Cells on Mitosis"    
    std::vector<std::shared_ptr<Cell<d>>> &getCellsOnMitosis(){return(cellsOnMitosis);};
    // Set the boundary vertices.
    void setBoundaryVertices();

    /// Every cell which preferred area A0 increases is included in a list of "Cells on Mitosis"    
    void addCellOnMitosis(std::shared_ptr<Cell<d>> &cell);
    // Add the given vertex to the list of boundary vertices.
    void addToBoundaryVertices(const std::shared_ptr<Vertex<d>> &v);
    void removeFromBoundaryVertices(const std::shared_ptr<Vertex<d>> &v);

    // ========== Line tensions ================ //
    // Set boundary line tension.
    void setBoundaryLambda(    double boundaryLambda_);

    // ========== Geometric functions ================ //

    /// Computed the area of the CellCluster
    Array<double,d> computeCentroid() const;
    /// Computed the area of the CellCluster
    double computeArea() const;

    // ========== dynamic functions ================ //

    /// advance all the vertices according to their 
    void advance(double dt);
    /// advance all the vertices according to their and reset force value
    void resetForce();
    /// Computes the force and update the value for each vertex
    void updateForce();

    /// Add a cell to the cellCluster
    void addCell(std::shared_ptr<Cell<d>> &cell){cells.push_back(cell);};
    /// Remove a cell from the cellCluster
    void remCell(std::shared_ptr<Cell<d>> &cell);
    /// Remove one edge from the CellCluster
    void removeEdge(std::shared_ptr<Edge<d>> &myedge);
    /// Remove one vertex from the CellCluster
    void removeVertex(std::shared_ptr<Vertex<d>> &myvertex);
    
    // Remove the connection between two cells through the given edge.
    bool removeConnectivityBetweenCellsThroughEdge( std::shared_ptr<Cell<d>> c1,
                                                                std::shared_ptr<Cell<d>> c2, 
                                                                std::shared_ptr<Edge<d>> e);
    // Add the connection between two cells through the given edge.
    void addConnectivityBetweenCellsThroughEdge(    std::shared_ptr<Cell<d>> c1,
                                                                std::shared_ptr<Cell<d>> c2, 
                                                                std::shared_ptr<Edge<d>> e);
    

    /// Remove all connectivities through with given cell from the cellCluster
    void removeConnectivity(std::shared_ptr<Cell<d>> &cell);
    /// Replace all connectivities of mycell with either a connectivity cell1 or cell2 from the cellCluster
    void replaceConnectivity(std::shared_ptr<Cell<d>> &mycell, std::shared_ptr<Cell<d>> &cell1, std::shared_ptr<Cell<d>> &cell2);
    


    
    /// Get the list of cells that have a preferred area A0 higher than its area A  
    std::vector<std::shared_ptr<Cell<d>>> getCellsToDivide();

private:
    std::vector<std::shared_ptr<Vertex<d>>> boundaryVertices;//the boundary vertices
    std::vector<std::shared_ptr<Cell<d>>> cellsOnMitosis;//the cells on Mitosis
    std::vector<double> initialAreaCOM;//the initial area of cells on Mitosis
    std::vector<std::shared_ptr<Cell<d>>> cells;      // complete list of cells in the tissue
    std::vector<std::shared_ptr<Edge<d>>> edges;     // complete list of edges in the tissue
    std::vector<std::shared_ptr<Vertex<d>>> vertices; // complete list of vertices in the tissue
    /// determines the neighborhood of each cell and which edge connects the different cells
    std::vector<std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>>> connectivity; 
    CellTypesManager cellTypesManager;
    double boundaryLambda;
};

// =========== HexagonalCells in a rectangular cluster ====================== //

template<int d>
class HexagonalCellsRectangularCluster final : public CellCluster<d> {
public:
    HexagonalCellsRectangularCluster(int nx, int ny, double circumRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &&dyn, bool isLarge, unsigned long iT);
    HexagonalCellsRectangularCluster(int nx, int ny, double circumRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &dyn, bool isLarge, unsigned long iT);
    HexagonalCellsRectangularCluster(int nx, int ny, double initCircumRadius, double K, double A0, double gamma, double lambda, std::shared_ptr<Dynamics<d>> &&dyn, bool isLarge);
    HexagonalCellsRectangularCluster(int nx, int ny, double initCircumRadius, double K, double A0, double gamma, double lambda, std::shared_ptr<Dynamics<d>> &dyn, bool isLarge);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");

};

// =========== Square Cells in a rectangular cluster ====================== //

template<int d>
class SquareCellsRectangularCluster final : public CellCluster<d> {
public:
    SquareCellsRectangularCluster(int nx, int ny, double circumRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &&dyn,unsigned long iT);
    SquareCellsRectangularCluster(int nx, int ny, double circumRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &dyn,unsigned long iT);
    SquareCellsRectangularCluster(int nx, int ny, double initCircumRadius, double K, double A0, double gamma, double lambda, std::shared_ptr<Dynamics<d>> &&dyn);
    SquareCellsRectangularCluster(int nx, int ny, double initCircumRadius, double K, double A0, double gamma, double lambda, std::shared_ptr<Dynamics<d>> &dyn);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");

};

// =========== HexagonalCells in a rectangular cluster ====================== //

template<int d>
class HexagonalCellsCircularCluster final : public CellCluster<d> {
public:
    HexagonalCellsCircularCluster(int n, double circumRadius, double K, double gamma, double lambda, std::shared_ptr<Dynamics<d>> &&dyn);
    HexagonalCellsCircularCluster(int n, double circumRadius, double K, double gamma, double lambda, std::shared_ptr<Dynamics<d>> &dyn);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");

};

// =========== Square Cells in a rectangular cluster ====================== //

template<int d>
class SquareCellsCircularCluster final : public CellCluster<d> {
public:
    SquareCellsCircularCluster(int n, double circumRadius, double K, double gamma, double lambda, std::shared_ptr<Dynamics<d>> &&dyn);
    SquareCellsCircularCluster(int n, double circumRadius, double K, double gamma, double lambda, std::shared_ptr<Dynamics<d>> &dyn);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");

};


// =========== Ring cluster of trapezoidal ====================== //

template <int d>
class TrapezoidalCellsRingCluster final : public CellCluster<d> {
public:
    TrapezoidalCellsRingCluster(int nCells, double cellHeight, double ringRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &&dyn, unsigned long iT);
    TrapezoidalCellsRingCluster(int nCells, double cellHeight, double ringRadius, double K, double gamma, double lambda, double sdGamma, double sdLambda, std::shared_ptr<Dynamics<d>> &dyn, unsigned long iT);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");
};

}  // end namespace epc

#endif
