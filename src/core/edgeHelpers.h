#ifndef EDGE_HELPERS_H
#define EDGE_HELPERS_H

#include "vertex.h"
#include "edge.h"

namespace epc {

/// equality of edges is the equality of the id of the vertices composing it
template<int d>
inline bool operator==(Edge<d> const& e1, Edge<d> const& e2) {
	auto id1_e1 = e1.getVertexOne()->getId();
	auto id2_e1 = e1.getVertexTwo()->getId();
	auto id1_e2 = e2.getVertexOne()->getId();
	auto id2_e2 = e2.getVertexTwo()->getId();

	if (((id1_e1 == id1_e2) && (id2_e1 == id2_e2)) || 
		((id1_e1 == id2_e2) && (id2_e1 == id1_e2)) ) 
	{
		return true;
	} else {
		return false;
	}
}

template<int d>
struct edgeHelpers {

static bool isInVector(Edge<d> const &e, const std::vector<std::shared_ptr<Edge<d>>> &edges) {
    for (auto eTmp : edges) {
        if (e == *eTmp) return true;
    }
    return false;
}

static void appendIfNotInVector(std::shared_ptr<Edge<d>> &e, std::vector<std::shared_ptr<Edge<d>>> &edges) {
    if (!isInVector(*e, edges)) edges.push_back(e);
}

};

}  // end namespace epc

#endif
