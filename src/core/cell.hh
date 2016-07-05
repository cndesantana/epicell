#ifndef CELL_HH
#define CELL_HH

#include "cell.h"
#include "vertex.h"
#include "edge.h"
#include "cellHelpers.h"
#include "vertexHelpers.h"
#include "mathematics/functions.h"

namespace epc {
// =============== Constuctors ========================= //
template<int d>
Cell<d>::Cell() : K(double()), A0(double()), gamma(double()), sdGamma(double())
{
    edges.resize(0);
    createNewId();
}

template<int d>
Cell<d>::Cell(std::vector<std::shared_ptr<Edge<d>>> &edges_, double K_, double A0_, double gamma_, double sdGamma_, unsigned long iT, bool onMitosis_) : 
        edges(edges_), K(K_), A0(A0_), gamma(gamma_), sdGamma(sdGamma_), onMitosis(onMitosis_),cellSignal(double(0.0)),oldCellSignal(double(0.0)),receivedSignal(double(0.0)){ 
    createNewId();
    setGamma(gamma_, sdGamma_);
    addBirthDate(iT);
}

template<int d>
Cell<d>::Cell(std::vector<std::shared_ptr<Edge<d>>> &edges_, double K_, double A0_, double gamma_, double sdGamma_, CellType type_, unsigned long iT, bool onMitosis_) : 
        edges(edges_), K(K_), A0(A0_), gamma(gamma_), sdGamma(sdGamma_), type(type_), onMitosis(onMitosis_),cellSignal(double(0.0)),oldCellSignal(double(0.0)),receivedSignal(double(0.0)){ 
    createNewId();
    setGamma(gamma_, sdGamma_);
    addBirthDate(iT);
}

template<int d>
void Cell<d>::createNewId() {
    static unsigned long allIds = 0; 
    id = allIds++;
}

template<int d>
void Cell<d>::addBirthDate(unsigned long iT) {
    birth = iT;
}

// ========== Interface functions ================ //
template<int d>
unsigned long const& Cell<d>::getId() const {
    return id;
}

template<int d>
unsigned long const& Cell<d>::getBirthDate() const {
    return birth;
}

template<int d>
std::vector<std::shared_ptr<Edge<d>>> const& Cell<d>::getEdges() const {
    return edges;
}

template<int d>
std::vector<std::shared_ptr<Edge<d>>> & Cell<d>::getEdges() {
    return edges;
}

template<int d>
double const& Cell<d>::getAreaBeforeMitosis() const {
    return areaBeforeMitosis;
}

template<int d>
double const& Cell<d>::getK() const {
    return K;
}
template<int d>
double const& Cell<d>::getKInit() const {
    return KInit;
}
template<int d>
double const& Cell<d>::getA0() const {
    return A0;
}
template<int d>
double const& Cell<d>::getAreaSteady() const {
    return areaSteady;
}
template<int d>
double const& Cell<d>::getPerimeterSteady() const {
    return perimeterSteady;
}

template<int d>
double const& Cell<d>::getGamma() const {
    return gamma;
}

template<int d>
double const& Cell<d>::getSdGamma() const {
    return sdGamma;
}

template<int d>
bool const& Cell<d>::isOnMitosis() const {
    return onMitosis;
}

template<int d>
CellType const& Cell<d>::getType() const {
    return type;
}

template<int d>
double const& Cell<d>::getSignal() const {
    return cellSignal;
}

template<int d>
double const& Cell<d>::getOldSignal() const {
    return oldCellSignal;
}

template<int d>
double const& Cell<d>::getReceivedSignal() const {
    return receivedSignal;
}

template<int d>
void Cell<d>::setAreaBeforeMitosis(double areaBeforeMitosis_){
    areaBeforeMitosis = areaBeforeMitosis_;
}

template<int d>
double const& Cell<d>::getGammaInit() const {
    return gammaInit;
}

template<int d>
void Cell<d>::setK(double K_) {
    K = K_;
}

template<int d>
void Cell<d>::setKInit(double KInit_) {
    KInit = KInit_;
}

template<int d>
void Cell<d>::setA0(double A0_) {
    A0 = A0_;
}
template<int d>
void Cell<d>::setAreaSteady(double areaSteady_) {
    areaSteady = areaSteady_;
}
template<int d>
void Cell<d>::setPerimeterSteady(double perimeterSteady_) {
    perimeterSteady = perimeterSteady_;
}

template<int d>
void Cell<d>::setGamma(double gamma_, double sdGamma_) {
    sdGamma = sdGamma_;
    gamma = NRNG(gamma,sdGamma_).ngenerate();
    // std::cout << "Gamma = " << gamma << std::endl;
}

template<int d>
void Cell<d>::setGammaInit(double gammaInit_) {
    gammaInit = gammaInit_;
}

template<int d>
void Cell<d>::setOnMitosis(bool onMitosis_) {
    onMitosis = onMitosis_;
}

template<int d>
void Cell<d>::setType(CellType type_) {
    type = type_;
}

template<int d>
void Cell<d>::setSignal(double cellSignal_) {
    cellSignal = cellSignal_;
}

template<int d>
void Cell<d>::setOldSignal(double oldCellSignal_) {
    oldCellSignal = oldCellSignal_;
}

template<int d>
void Cell<d>::setReceivedSignal(double receivedSignal_) {
    receivedSignal = receivedSignal_;
}

template<int d>
void Cell<d>::addToReceivedSignal(double receivedSignal_) {
    receivedSignal += receivedSignal_;
}

// ========== geometric functions ================ //
template<int d>
Array<double,d> Cell<d>::computeCentroid() const {
    Array<double,d> centroid(double(0.0));
    for ( auto edge : edges ) {
        centroid += edge->getVertexOne()->getPosition();
        centroid += edge->getVertexTwo()->getPosition();
    }
    centroid /= (double)edges.size();
    return (double)0.5*centroid; // all vertices are added twice (must divide by 2)
}

template<int d>
double Cell<d>::computeDegree() {
    double degree = double (edges.size());

    return degree;
}


template<int d>
double Cell<d>::computeArea() const {
    auto centroid = computeCentroid();
    double area = double();
    for ( auto edge : edges ) {
        Array<double,d> v1 = edge->getVertexOne()->getPosition()-centroid;
        Array<double,d> v2 = edge->getVertexTwo()->getPosition()-centroid;

        area += arrayOperations<double,d>::crossProductNorm(v1,v2);
    }

    return area*(double)0.5;
}

template<int d>
void Cell<d>::computeAndUpdateForce() {
    auto centroid = computeCentroid();
    auto area = double();
    auto perimeter = double();

    std::vector<int> sign(edges.size());
    std::vector<Array<double,d>> vecAlongEdges(edges.size());
    std::vector<Array<double,d>> vecNormEdges(edges.size());
    std::vector<double> lengths(edges.size());
    for ( unsigned iA = 0; iA < edges.size(); ++iA ) { // for each edge in cell
        Array<double,d> v1 = edges[iA]->getVertexOne()->getPosition()-centroid;
        Array<double,d> v2 = edges[iA]->getVertexTwo()->getPosition()-centroid;

        vecAlongEdges[iA] = v1-v2;
        lengths[iA] = vecAlongEdges[iA].norm(); // used for the computation of the unit vector for perimeter force

        auto tmpA = arrayOperations<double,d>::crossProduct(v1,v2);
        sign[iA] = fun::sgn(tmpA); // used for the computation of the direction of the area force

        perimeter += lengths[iA]; // compute perimeter
        area += (double)sign[iA]*tmpA; // compute area
    }
    area *= 0.5;

    double preFactorArea = -0.5*K*(area-A0);
    for ( unsigned iA = 0; iA < edges.size(); ++iA ) {
        Array<double,d> areaForce = preFactorArea * arrayOperations<double,d>::rotateCounterClockwise(vecAlongEdges[iA]) * (double)sign[iA]; // TODO make sur the signs are correct
        edges[iA]->getVertexOne()->addToForce(areaForce); // area force terms
        edges[iA]->getVertexTwo()->addToForce(areaForce);

        Array<double,d> perimForce = -vecAlongEdges[iA]/lengths[iA]*gamma*perimeter; // TODO make sur the signs are correct
        edges[iA]->getVertexOne()->addToForce(perimForce); // perimeter force terms
        edges[iA]->getVertexTwo()->addToForce(-perimForce);
    }
}

template<int d>
double Cell<d>::computePerimeter() {
    double perimeter = double(0.0);
    for ( auto edge : edges ) {
        perimeter += edge->computeLength();
    }     
    return perimeter;
}

/// Add an edge to the Cell
template<int d>
void Cell<d>::addEdge(std::shared_ptr<Edge<d>> &myedge){
	edges.push_back(myedge);
}     

/// Remove an edge from the Cell
template<int d>
void Cell<d>::removeEdge(std::shared_ptr<Edge<d>> &myedge){
    auto nEdges = edges.size();
    for(int i=0;i<nEdges;i++){
        if (edges[i]->getId() == myedge->getId()){
                edges.erase(edges.begin()+i);
                break;
        }
    }
}    

}  // end namespace epc

#endif
