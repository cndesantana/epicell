#ifndef VECTOR_HELPERS_H
#define VECTOR_HELPERS_H

#include <iostream>
#include <cmath>
#include <vector>
#include "core/epcDebug.h"

namespace epc {

template <typename T>
std::vector<T> operator+(std::vector<T> const& a, std::vector<T> const& b) {
    EPC_ASSERT(a.size() == b.size() && "size of a must be equal to b for addition");
    std::vector<T> result;
    for (unsigned i=0; i<a.size; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator+(std::vector<T> const& a, T alpha) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = a[i] + alpha;
    }
    return result;
}

template<typename T>
std::vector<T> operator+(T alpha, std::vector<T> const& a) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = alpha + a[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator-(std::vector<T> const& a, std::vector<T> const& b) {
    EPC_ASSERT(a.size() == b.size() && "size of a must be equal to b for subtraction");
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator-(std::vector<T> const& a) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = -a[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator-(std::vector<T> const& a, T alpha) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = a[i] - alpha;
    }
    return result;
}

template<typename T>
std::vector<T> operator-(T alpha, std::vector<T> const& a) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = alpha - a[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator*(std::vector<T> const& a, std::vector<T> const& b) {
    EPC_ASSERT(a.size() == b.size() && "size of a must be equal to b for elementwise multiplication");
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator*(std::vector<T> const& a, T alpha) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = a[i] * alpha;
    }
    return result;
}

template<typename T>
std::vector<T> operator*(T alpha, std::vector<T> const& a) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = alpha * a[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator/(std::vector<T> const& a, std::vector<T> const& b) {
    EPC_ASSERT(a.size() == b.size() && "size of a must be equal to b for elementwise division");
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator/(std::vector<T> const& a, T alpha) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = a[i] / alpha;
    }
    return result;
}

template<typename T>
inline bool operator==(std::vector<T> const& A, std::vector<T> const& B) {
    EPC_ASSERT(A.size() == B.size() && "size of A must be equal to B for equality test");
    for (unsigned iA = 0; iA < A.size(); ++iA) {
        if (A[iA] != B[iA]) return false;
    }
    return true;
}

template<typename T>
inline bool operator!=(std::vector<T> const& A, std::vector<T> const& B) {
    EPC_ASSERT(A.size() == B.size() && "size of A must be equal to B for inequality test");
    for (unsigned iA = 0; iA < A.size(); ++iA) {
        if (A[iA] != B[iA]) return true;
    }
    return false;
}

template<typename T>
inline bool operator<(std::vector<T> const& A, std::vector<T> const& B) {
    EPC_ASSERT(A.size() == A.size() && "size of A must be equal to B for smaller than test");
    for (unsigned iA=0; iA<A.size(); ++iA) {
        if (A[iA] < B[iA]) return true;
        if (A[iA] > B[iA]) return false;
    }
    return false;
}

template<typename T>
std::vector<T> operator/(T alpha, std::vector<T> const& a) {
    std::vector<T> result;
    for (unsigned i=0; i<a.size(); ++i) {
        result[i] = alpha / a[i];
    }
    return result;
}




} // end namespace epc

#endif

