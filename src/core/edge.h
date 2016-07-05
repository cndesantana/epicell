#ifndef EDGE_H
#define EDGE_H

#include "vertex.h"

namespace epc {

template<int d>
class Edge {
public:
    /// Default constructor (everything is set to zero)
    Edge();
    Edge(std::shared_ptr<Vertex<d> > &vertexOne_, std::shared_ptr<Vertex<d> > &vertexTwo_, double lambda_, double sdLambda_);
private:
    static_assert(d <= 3 && d >= 2, "d must be 2 or 3!");
    /// increments the id of a cell
    void createNewId();
public:
    // ========== Interface functions ================ //
    /// Get a reference to non-modifiable id
    unsigned long const& getId() const;

    /// Get a reference to non-modifiable pointer to vertex 1
    std::shared_ptr<Vertex<d>> const& getVertexOne() const;
    /// Get a reference to modifiable pointer to vertex 1
    std::shared_ptr<Vertex<d>> & getVertexOne();

    /// Get a reference to modifiable pointer to vertex 2
    std::shared_ptr<Vertex<d>> & getVertexTwo();
    /// Get a reference to non-modifiable pointer to vertex 2
    std::shared_ptr<Vertex<d>> const& getVertexTwo() const;

    /// Get a non-modifiable reference to the mean of lambda distribution
    double const& getLambda() const;
    /// Get a non-modifiable reference to standard deviation of lambda distribution
    double const& getSdLambda() const;

    /// Set a value to lambda
    void setLambda(double lambda_,double sdLambda_ = 0);
    /// Set an Id
    void setId(unsigned long lambda_);
   

    // ========== geometric functions ================ //
    /// Computed the length of the edge
    double computeLength();
    /// Compute vector from vertex 1 to 2
    Array<double,d> computeVectorFromOneToTwo();
    /// Compute vector from vertex 2 to 1
    Array<double,d> computeVectorFromTwoToOne();
    // ========== dynamic functions ================ //
    void computeAndUpdateForce();
private:
    std::tuple<std::shared_ptr<Vertex<d> > , std::shared_ptr<Vertex<d> > > vertices;
    double lambda; // constant linked with the force related to the edge
    double sdLambda; // constant linked with the variability of lambda
    unsigned long id; // unique identifier for each Edge TODO: rethink that maybe
};



}  // end namespace epc

#endif
