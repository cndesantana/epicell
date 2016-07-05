#ifndef SIGNALING_HELPERS_H
#define SIGNALING_HELPERS_H


#include "core/cellCluster.h"

namespace epc {
  
template<int d>
struct signalingHelpers {

static bool isSignalDistributionStable(const std::vector<std::shared_ptr<Cell<d>>> & cells, double signal_eps ) {
    for (auto c: cells){
        if (fabs(c->getOldSignal() - c->getSignal()) > signal_eps){
            return false;
        }
    }
    return true;
}


};


}  // end namespace epc

#endif
