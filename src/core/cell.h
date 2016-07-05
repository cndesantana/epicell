#ifndef CELL_H
#define CELL_H

#include "vertex.h"
#include "edge.h"
#include "cellTypes/cellType.h"
#include "cellHelpers.h"
#include "vertexHelpers.h"
#include "mathematics/functions.h"

namespace epc {

template<int d>
class Cell {
public:
    /// Default constructor (everything is set to zero)
    Cell();
    Cell(std::vector<std::shared_ptr<Edge<d>>> &edges_, double K_, double A0_, double gamma_, double sdGamma_, unsigned long iT,  bool onMitosis_ = false);
    Cell(std::vector<std::shared_ptr<Edge<d>>> &edges_, double K_, double A0_, double gamma_, double sdGamma_, CellType type_, unsigned long iT,  bool onMitosis_ = false);

private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");
    /// increments the id of a cell
    void createNewId();
    void addBirthDate(unsigned long iT);
public:
    // ========== Interface functions ================ //
    /// Get a reference to non-modifiable id
    unsigned long const& getId() const;
    unsigned long const& getBirthDate() const;
    double const& getAreaBeforeMitosis() const;

    /// get a non-modifiable reference to edges
    std::vector<std::shared_ptr<Edge<d>>> const &getEdges() const;
    /// get a modifiable reference to edges
    std::vector<std::shared_ptr<Edge<d>>> &getEdges();
    std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>> getNeighbors(){
        return neighbors;
    }
    /// Add an edge to the cellCluster
    void addEdge(std::shared_ptr<Edge<d>> &edge);
    /// Remove one edge from the cellCluster
    void removeEdge(std::shared_ptr<Edge<d>> &edge);
    /// Get a non-modifiable reference to K
    double const& getK() const;
    /// Get a non-modifiable reference to A0
    double const& getA0() const;
    /// Get a non-modifiable reference to the mean of gamma distribution
    double const& getGamma() const;
    /// Get a non-modifiable reference to the standard deviation of gamma distribution
    double const& getSdGamma() const;
    /// Get a non-modifiable reference to areaSteady
    double const& getAreaSteady() const;
    /// Get a non-modifiable reference to perimeterSteady
    double const& getPerimeterSteady() const;
    /// Get a non-modifiable reference to gamma
    double const& getGammaInit() const;
    /// Get a non-modifiable reference to KInit
    double const& getKInit() const;
    /// Check if the cell is on mitosis. The retrieved reference is non-modifiable.
    bool const& isOnMitosis() const;
    /// Get a non-modifiable reference to type
    CellType const& getType() const;
    /// Get a non-modifiable reference to signal
    double const& getSignal() const;
    /// Get a non-modifiable reference to oldSignal
    double const& getOldSignal() const;
    /// Get a non-modifiable reference to receivedSignal
    double const& getReceivedSignal() const;

    /// Set Cell Area Before Mitosis 
    void setAreaBeforeMitosis(double areaBeforeMitosis_);
    /// Set a value to K
    void setK(double K_);
    /// Set a value to A0
    void setA0(double A0_);
    /// Set a value to A0
    void setAreaSteady(double areaSteady_);
    /// Set a value to perimeterSteady
    void setPerimeterSteady(double perimeterSteady_);
    /// Set a value to gamma (and sd)
    void setGamma(double gamma_, double sdGamma_ = 0.0);
    /// Set a value to gammaINitit
    void setGammaInit(double KInit_);
    /// Set a value to gammaInit
    void setKInit(double KInit_);
    /// Set a value to onMitosis
    void setOnMitosis(bool onMitosis_);
    /// Set a value to type
    void setType(CellType type_);
    /// Set a value to signal
    void setSignal(double signal_);
    /// Set a value to oldSignal
    void setOldSignal(double oldSignal_);
    /// Set a value to receivedSignal
    void setReceivedSignal(double receivedSignal_);
    void addToReceivedSignal(double receivedSignal_);
    // ========== geometric functions ================ //
    /// Computed the centroid (geometric center) of the Cell
    Array<double,d> computeCentroid() const;
    /// Computed the area (geometric center) of the Cell
    double computeArea() const;
    /// Computed the degree of the cell
    double computeDegree();
    /// Compute vector from vertex 2 to 1
    double computePerimeter();

    // =============== dynamics functions ================== //
    void computeAndUpdateForce();
private:
    std::vector<std::shared_ptr<Edge<d>>> edges;
    std::vector<std::tuple<std::shared_ptr<Cell<d>>,std::shared_ptr<Edge<d>>>> neighbors; // a neighbhor is a cell and the border edge
    double K, A0, gamma, sdGamma; // K(strength of area force), A0(preferref area), gamma(strengh of perim force)
    double KInit, gammaInit, perimeterSteady, areaSteady;
    CellType type;
    unsigned long id; // unique identifier for each Cell TODO: rethink that maybe
    unsigned long birth; //birth date of the cell
    bool onMitosis;
    double areaBeforeMitosis;
    
    double cellSignal;
    double oldCellSignal;
    double receivedSignal;
};


}  // end namespace epc

#endif
