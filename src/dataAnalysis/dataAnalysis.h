#ifndef DATA_ANALYSIS_H
#define DATA_ANALYSIS_H

#include "mathematics/constants.h"
#include "dataAnalysis/statistics.h"
#include "core/cellCluster.h"
#include <cmath>

namespace epc {

namespace dataAnalysis {

// ============== vertex list analysis ================== //
template<int d>
Array<double,d> computeCentroid(const std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    Array<double,d> centroid = 0.0;
    for (const auto &v : vertices) {
        centroid += v->getPosition();
    }
    return centroid / (double)vertices.size();
}

template<int d>
double computeTotalForceNorm(const std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    double force = 0.0;
    for (const auto &v : vertices) {
        force += v->getForce().norm();
    }
    return force;
}

template<int d>
double computeAverageForceNorm(const std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    return computeTotalForceNorm(vertices) / vertices.size();
}

template<int d>
Array<double,d> computeNetForce(const std::vector<std::shared_ptr<Vertex<d>>> &vertices) {
    Array<double,d> force(double(0)); 
    for (const auto &v : vertices) {
        force += v->getForce();
    }
    return force;
}

template<int d>
double computeNetEnergy(const CellCluster<d> &cellCluster) {
   double energy = double(0); 
    for (const auto &c : cellCluster.getCells()) {
        energy += 0.5*c->getK()*(c->computeArea()-c->getA0())*(c->computeArea()-c->getA0());
        energy += 0.5*c->getGamma()*c->computePerimeter()*c->computePerimeter();
    }
    for (const auto &e : cellCluster.getEdges()) {
        energy += e->getLambda()*e->computeLength();
    }
    return energy;
}

// ============== cell list analysis ================== //
template<int d>
double computeArea(const std::vector<std::shared_ptr<Cell<d>>> &cells) {
    double area = 0.0; 
    for (const auto &c : cells) {
        area += c->computeArea();
    }
    return area;
}

template<int d>
double computeAverageCellArea(const std::vector<std::shared_ptr<Cell<d>>> &cells) {
    return computeArea(cells)/cells.size();
}

template<int d>
double computeAverageCellPerimeter(const std::vector<std::shared_ptr<Cell<d>>> &cells) {
    double avgPerimeter = 0.0; 
    for (const auto &c : cells) {
        avgPerimeter += c->computePerimeter();
    }
    return avgPerimeter/cells.size();
}

// ============== cell cluster analysis ================== //
template<int d>
double computeTotalForceNorm(const CellCluster<d> &cellCluster) {
    return computeTotalForceNorm(cellCluster.getVertices());
}

template<int d>
double computeAverageForceNorm(const CellCluster<d> &cellCluster) {
    return computeAverageForceNorm(cellCluster.getVertices());
}

template<int d>
Array<double,d> computeNetForce(const CellCluster<d> &cellCluster) {
    return computeNetForce(cellCluster.getVertices());
}

template<int d>
double computeArea(const CellCluster<d> &cellCluster) {
    return computeArea(cellCluster.getCells());
}

template<int d>
double computeAverageCellArea(const CellCluster<d> &cellCluster) {   
    return computeArea(cellCluster.getCells())/cellCluster.getCells().size();
}

template<int d>
double computeAverageCellPerimeter(const CellCluster<d> &cellCluster) {   
    return computeAverageCellPerimeter(cellCluster.getCells());
}

template<int d>
double computeRegularness(const CellCluster<d> &cellCluster) {   
    auto nEdges = cellCluster.getEdges().size();
    statistics::Convergence conv(double(mathConstants::ForceEpsilon), nEdges);
    for(auto iE: cellCluster.getEdges()){
        conv.takeValue(iE->computeLength());//get the length of the edges
    }
    return(conv.getNormalizedStd()); //get the std/mean of the length of the edges
}

template<int d>
std::vector<double> computeDistributionOfSides(const CellCluster<d> &cellCluster) {   
    int sizeOfVec = 7;//7 possible sides: 3 and lower, 4, 5, 6, 7, 8, 9 and higher
    std::vector<int> vecOfSides(sizeOfVec, 0);
    std::vector<double> vecOfProportion(sizeOfVec, 0.0);
    auto ncells = cellCluster.getCells().size();
    for(auto iC: cellCluster.getCells()){
        auto nsides = iC->getEdges().size();//a value higher or equal to 3
        if(nsides <= 9){//nsides 3 refers to position 0, nsides 4 refers to position 1, ..., nsides 9 refers to position 6
           vecOfSides[nsides-3]++;
        }
        else{//nsides higher than 9 refers to position 6
           vecOfSides[6]++;
        }
    }
    for(auto iV=0; iV < sizeOfVec; iV++){
        vecOfProportion[iV] = double(vecOfSides[iV])/ncells; 
    }

    return(vecOfProportion); //get the std/mean of the length of the edges
}

} // end namespace dataAnalysis


}  // end namespace epc

#endif
