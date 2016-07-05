#ifndef EPITHELIUM_H
#define EPITHELIUM_H

#include "vertex.h"
#include "cellCluster.h"
#include "epitheliumHelpers.h"


namespace epc {

template<int d>
class Epithelium {
public:
    Epithelium(CellCluster<d> &&cellCluster_);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");
public:
    // ========== Interface functions ================ //
    /// get a non-modifiable reference to edges
    CellCluster<d> const &getCellCluster() const;
    /// get a modifiable reference to edges
    CellCluster<d> &getCellCluster();
    /// get a non-modifiable reference to the vector of moved vertices
    std::vector<std::vector<std::shared_ptr<Vertex<d>>>> const &getMovedVertices() const;
    /// get a non-modifiable reference to the vector of forced vertices
    std::vector<std::vector<std::shared_ptr<Vertex<d>>>> const &getForcedVertices() const;
    /// get a non-modifiable reference to a component of moved vertices
    std::vector<std::shared_ptr<Vertex<d>>> const &getMovedVertices(unsigned iA) const;
    /// get a non-modifiable reference to a component of forced vertices
    std::vector<std::shared_ptr<Vertex<d>>> const &getForcedVertices(unsigned iA) const;
    /// add vertices and function to special vertices and force function vectors
    template<class T, class U>
    void addSpecialVerticesAndForceFunction(T &&vertices, U &&foo);
    // add vertices and function to forced vertices and force function vectors
    template<class T, class U>
    void addSpecialVerticesAndMovementFunction(T &&vertices, U &&foo);
    // add vertices which must keep their Y position unchanged
    template<class T>
    void addFixedOnYVertices(T &&vertices);
    /// remove all moved vertices and movement function vectors
    void removeAllMovedVerticesAndMovementFunction();
    /// remove all forced vertices and force function vectors
    void removeAllForcedVerticesAndForceFunction();
     /// remove all vertices which must keep their Y position unchanged
    void removeAllFixedOnYVertices();
    // ========== dynamic functions ================ //
    /// update the force, advance the vertices 
    void performTimeStep(double dt);
    /// compute and update the force of all vertices
    void updateForce();
    /// reset the force
    void resetTheForce();
    /// update the cell mechanics of all edges and all cells 
    void updateCellMechanics(double gamma, double lambda, double sdGamma, double sdLambda);
    
private:
    // ========== dynamic functions ================ //
    /// add force which is a function to the special vertices to be submitted to an external force
    void addForceOnSpecialVertices();
    /// set a null force on the moved vertices
    void setNullForceOnMovedVertices();
    /// set a null force along direction Y on the vertices that should keep their position Y.
    void setNullForceYOnFixedYVertices();
    /// add position movement which is a function to the special vertices to be moved
    void moveSpecialVertices(double dt);
private:
    /// cellCluster composing the epithelium
    CellCluster<d> cellCluster;
    /// complete list of forced vertices in the epithelium ()
    std::vector<std::vector<std::shared_ptr<Vertex<d>>>> forcedVertices; 
    /// complete list of moved vertices in the epithelium ()
    std::vector<std::vector<std::shared_ptr<Vertex<d>>>> movedVertices; 
    /// complete list of vertices whch must stay "fixed along the Y" direction in the epithelium ()
    std::vector<std::shared_ptr<Vertex<d>>> fixedOnYVertices; 
    /// complete list of force functions that act on the special vertices
    std::vector<std::function<Array<double,2> (const Vertex<2> &)>> forceFunctions;
    /// complete list of movement functions that act on the special vertices
    std::vector<std::function<Array<double,2> (const Vertex<2> &)>> movementFunctions;
};


}  // end namespace epc

#endif
