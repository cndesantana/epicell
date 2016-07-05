#ifndef _CELLSELECTIONHELPERS_
#define _CELLSELECTIONHELPERS_

#include "mathematics/array.h"
#include "core/cell.h"
#include "core/vertex.h"
#include "core/edge.h"
#include "core/cellHelpers.h"
#include "core/vertexHelpers.h"
#include <vector> 
#include <iostream>
#include <algorithm>

namespace epc {

template<int d>
class CellSelectionHelpers {
public:
    std::vector<float> CumulativeSum(std::vector<unsigned long> sequence);
    std::vector<float> CumulativeSum(std::vector<double> sequence);
    std::vector<unsigned long> getAgeAllCells(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT);
    std::vector<double> getSignalAllCells(const std::vector<std::shared_ptr<Cell<d>>> &cells);
    std::vector<double> getSqrSignalAllCells(const std::vector<std::shared_ptr<Cell<d>>> &cells);
    std::vector<double> getCubicSignalAllCells(const std::vector<std::shared_ptr<Cell<d>>> &cells);
    int Multinomial(std::vector<float> Vc,float myrand);
};

template<int d>
std::vector<float> CellSelectionHelpers<d>::CumulativeSum(std::vector<unsigned long> sequence){
    std::vector<float> res (sequence.size());
    
    res[0] = (float) sequence[0];
    float previous;

    for(auto pos=1;pos<(long int)sequence.size();pos++){
        res[pos] += previous + (float)sequence[pos];
        previous=res[pos];
    }        
    auto maximum = *std::max_element(res.begin(),res.end());
    for(auto i=0;i<(int)res.size();i++){
    	res[i] = res[i]/maximum;
    }		

    return res;
}

template<int d>
std::vector<float> CellSelectionHelpers<d>::CumulativeSum(std::vector<double> sequence){
    std::vector<float> res(sequence.size());
    
    res[0] = (float)sequence[0];
    float previous = res[0];
    for(auto pos=1;pos<(long int)sequence.size();pos++){
       res[pos] += previous + (float)sequence[pos];
       previous=res[pos];
    }        
    float maximum = previous;

    if (maximum == 0){
        return res;
    }

    for(auto i=0;i<(int)res.size();i++){
        res[i] = res[i]/maximum;
    }       

    return res;
}

/// CHARLES
/*template<int d>
int CellSelectionHelpers<d>::Multinomial(std::vector<float> Vc,float myrand){

	int pos = 0;
	for(int i=0;i<Vc.size();i++){
		if(myrand >= Vc[i]){
			pos = i;
		}
	}
	return(pos);
}*/

/// AZIZA
template<int d>
int CellSelectionHelpers<d>::Multinomial(std::vector<float> Vc,float myrand){

    int i = 0;
    while(myrand > Vc[i]){
        i += 1;
    }
    return(i);
}

template<int d>
std::vector<unsigned long> CellSelectionHelpers<d>::getAgeAllCells(const std::vector<std::shared_ptr<Cell<d>>> &cells,unsigned long iT){
    std::vector<unsigned long> res;

    for(auto c: cells){
	   res.push_back(c->getBirthDate()); 
    }        
    for(auto i=0;i<(long int)res.size();i++){//so the oldest cells will have the the highest age. And the minimum age is 1
	   res[i] = iT-res[i] + 1; 
    }        

    return res;
}

template<int d>
std::vector<double> CellSelectionHelpers<d>::getSignalAllCells(const std::vector<std::shared_ptr<Cell<d>>> &cells){
    std::vector<double> res;
    for(auto c: cells){
       res.push_back(c->getSignal()); 
    }              
    return res;
}

template<int d>
std::vector<double> CellSelectionHelpers<d>::getSqrSignalAllCells(const std::vector<std::shared_ptr<Cell<d>>> &cells){
    std::vector<double> res;
    for(auto c: cells){
       res.push_back(pow(c->getSignal(),2)); 
    }              
    return res;
}

template<int d>
std::vector<double> CellSelectionHelpers<d>::getCubicSignalAllCells(const std::vector<std::shared_ptr<Cell<d>>> &cells){
    std::vector<double> res;
    for(auto c: cells){
       res.push_back(pow(c->getSignal(),3)); 
    }              
    return res;
}

}//end namespace epd

#endif
