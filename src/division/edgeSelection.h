#ifndef EDGE_SELECTION_H
#define EDGE_SELECTION_H

#include "core/epcDebug.h"
#include "core/edge.h"
#include "mathematics/rng.h"

namespace epc {

template<int d>
class 	EdgeSelection {
public:
//	virtual int select(std::vector<std::shared_ptr<Edge<d>>> &cells) = 0;
	virtual ~EdgeSelection(){};
	virtual int select(const std::shared_ptr<Cell<d>> &cell) = 0;
};

template<int d>
class RandomEdgeSelection final : public EdgeSelection<d> {
public:
	virtual int select(const std::shared_ptr<Cell<d>> &cell) override;
//	virtual int select(std::vector<std::shared_ptr<Edge<d>>> &cells) override;
};

template<int d>
int RandomEdgeSelection<d>::select(const std::shared_ptr<Cell<d>> &cell) {
//int RandomEdgeSelection<d>::select(std::vector<std::shared_ptr<Edge<d>>> &cells) {
	// randomly choose a cell and tell it to start the division...
	int nedges = cell->getEdges().size();

    	RNG myirn(0,(nedges-1));
	int edgeSelected = myirn.igenerate();//generate a random integer number

//	std::cout << "Random edge selected: e" << cell->getEdges()[edgeSelected]->getId() << std::endl;
//	std::cout << "IN: Edge e" << cell->getEdges()[edgeSelected]->getId() << " = (v" << cell->getEdges()[edgeSelected]->getVertexOne()->getId() << ",v" << cell->getEdges()[edgeSelected]->getVertexTwo()->getId() << ")" << std::endl;
	// now to increase the size of "cellSelected".	
	return(edgeSelected);
}

template<int d>
class EldestEdgeSelection final : public EdgeSelection<d> {
public:
	virtual int select(const std::shared_ptr<Cell<d>> &cell) override;
//	virtual int select(std::vector<std::shared_ptr<Edge<d>>> &cells) override;
};

template<int d>
int EldestEdgeSelection<d>::select(const std::shared_ptr<Cell<d>> &cell) {
//int EldestEdgeSelection<d>::select(std::vector<std::shared_ptr<Edge<d>>> &cells) {
	// chose a cell according to its age (eldest cell first) and tell it to start the division...
	return(0);
}


}  // end namespace epc

#endif
