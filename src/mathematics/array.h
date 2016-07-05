#ifndef ARRAY_H
#define ARRAY_H

#include <iostream>
#include <cmath>
#include <vector>
#include <limits>
#include "core/epcDebug.h"


namespace epc {

/** A simple array class, which is slightly more convenient
 *  to use than the pure C-array. It can be assigned from an
 *  Array of different type, contains bound verifications for
 *  debugging.
 */
template<typename T, unsigned size>
class Array {
public:
    Array() { }

    Array(T v1, T v2) { 
        EPC_ASSERT(size == 2);
        data[0] = v1;
        data[1] = v2;
    }

    Array(T const* cArray) { 
        from_cArray(cArray);
    }

    Array(const std::vector<T> &vec) { 
        from_cVector(vec);
    }

    Array(T val) { 
        resetTo(val);
    }

    /// Copy-construction is only allowed with same-size Arrays,
    /// to prevent bugs due to cutting off one of the arrays.
    template<typename U>
    Array(Array<U,size> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+size, data);
    }
    /// Assignment is only allowed with same-size Arrays,
    /// to prevent bugs due to cutting off one of the arrays.
    template<typename U>
    Array<T,size>& operator=(Array<U,size> const& rhs) {
        std::copy(&rhs[0], &rhs[0]+size, data);
        return *this;
    }
    T& operator[](unsigned index) {
        EPC_ASSERT((index >= 0  && index<size) && "the index must be smaller than the dimension of the array and bigger than 0");
        return data[index];
    }
    T const& operator[](unsigned index) const {
        EPC_ASSERT((index >= 0  && index<size) && "the index must be smaller than the dimension of the array and bigger than 0");
        return data[index];
    }

    void from_cArray(T const* cArray) {
        std::copy(cArray, cArray+size, data);
    }
    void from_cVector(const std::vector<T> &vec) {
        data = &vec[0];
    }
    void add_from_cArray(T const* cArray) {
        for (unsigned i=0; i<size; ++i) {
            data[i] += cArray[i];
        }
    }
    void add_from_cArray(T const* cArray, T factor) {
        for (unsigned i=0; i<size; ++i) {
            data[i] += factor*cArray[i];
        }
    }
    void to_cArray(T* cArray) const {
        std::copy(data, data+size, cArray);
    }
    void add_to_cArray(T* cArray) const {
        for (unsigned i=0; i<size; ++i) {
            cArray[i] += data[i];
        }
    }
    void add_to_cArray(T* cArray, T factor) const {
        for (unsigned i=0; i<size; ++i) {
            cArray[i] += factor*data[i];
        }
    }
    void resetTo(T val = T()) {
        std::fill_n(data, size, val);
    }
    Array<T,size>& operator += (Array<T,size> const& b) {
        for (unsigned i=0; i<size; ++i) {
            data[i] += b[i];
        }
        return *this;
    }
    Array<T,size>& operator += (T alpha) {
        for (unsigned i=0; i<size; ++i) {
            data[i] += alpha;
        }
        return *this;
    }
    Array<T,size>& operator -= (Array<T,size> const& b) {
        for (unsigned i=0; i<size; ++i) {
            data[i] -= b[i];
        }
        return *this;
    }
    Array<T,size>& operator -= (T alpha) {
        for (unsigned i=0; i<size; ++i) {
            data[i] -= alpha;
        }
        return *this;
    }
    Array<T,size>& operator *= (Array<T,size> const& b) {
        for (unsigned i=0; i<size; ++i) {
            data[i] *= b[i];
        }
        return *this;
    }
    Array<T,size>& operator *= (T alpha) {
        for (unsigned i=0; i<size; ++i) {
            data[i] *= alpha;
        }
        return *this;
    }
    Array<T,size>& operator /= (Array<T,size> const& b) {
        for (unsigned i=0; i<size; ++i) {
            data[i] /= b[i];
        }
        return *this;
    }
    Array<T,size>& operator /= (T alpha) {
        for (unsigned i=0; i<size; ++i) {
            data[i] /= alpha;
        }
        return *this;
    }

    T normSqr() const {
        T nSqr = data[0]*data[0];
        for (unsigned iD = 1; iD < size; ++iD) {
            nSqr += data[iD]*data[iD];
        }
        return nSqr;
    }
    T norm() const {
        return sqrt(normSqr());
    }
private:
    T data[size];
};

// =========== array-array operations ==================== //

template<typename T, unsigned size> struct arrayOperationsImpl;

template <typename T, unsigned size>
struct arrayOperations {
    static T crossProductNorm(Array<T,size> const& a, Array<T,size> const& b) {
        return arrayOperationsImpl<T,size>::crossProductNorm(a, b);
    }

    static T crossProduct(Array<T,size> const& a, Array<T,size> const& b) {
        return arrayOperationsImpl<T,size>::crossProduct(a, b);
    }

    static void crossProduct(Array<T,size> const& a, Array<T,size> const& b, Array<T,size> &res) {
        arrayOperationsImpl<T,size>::crossProduct(a, b, res);
    }

    static Array<T,size> rotateCounterClockwise(Array<T,size> const& a) {
        return arrayOperationsImpl<T,size>::rotateCounterClockwise(a);
    }

    static bool isEqual(Array<T,size> const& a, Array<T,size> const& b, T eps = 100.0*std::numeric_limits<T>::epsilon()) {
        for (unsigned i=0; i<size; ++i) {
            if (std::fabs(a[i]-b[i]) > eps) {
                return false;
            }
        }
        return true;
    }
};

// =========== array-array operations ==================== //

template <typename T, unsigned size>
Array<T,size> operator+(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}


template<typename T, unsigned size>
Array<T,size> operator-(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<typename T, unsigned size>
Array<T,size> operator-(Array<T,size> const& a) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = -a[i];
    }
    return result;
}

template<typename T, unsigned size>
Array<T,size> operator*(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template<typename T, unsigned size>
Array<T,size> operator/(Array<T,size> const& a, Array<T,size> const& b) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}

template<typename T, unsigned size>
inline bool operator==(Array<T,size> const& A, Array<T,size> const& B) {
    for (unsigned iA = 0; iA < size; ++iA) {
        if (A[iA] != B[iA]) return false;
    }
    return true;
}

template<typename T, unsigned size>
inline bool operator!=(Array<T,size> const& A, Array<T,size> const& B) {
    for (unsigned iA = 0; iA < size; ++iA) {
        if (A[iA] != B[iA]) return true;
    }
    return false;
}

template<typename T, unsigned size>
inline bool operator<(Array<T,size> const& A, Array<T,size> const& B) {
    for (unsigned iA=0; iA<size; ++iA) {
        if (A[iA] < B[iA]) return true;
        if (A[iA] > B[iA]) return false;
    }
    return false;
}

template<typename T, unsigned size>
inline bool operator>(Array<T,size> const& A, Array<T,size> const& B) {
    for (unsigned iA=0; iA<size; ++iA) {
        if (A[iA] > B[iA]) return true;
        if (A[iA] < B[iA]) return false;
    }
    return false;
}

// ======================================================== //
// =========== array-vector operations ==================== //
// ======================================================== //

// =================== addition =========================== //

template <typename T, unsigned size>
Array<T,size> operator+(Array<T,size> const& a, std::vector<T> const& b) {
    ELC_ASSERT(b.size() == size && "size of the array must be equal to the size of the vector");
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template <typename T, unsigned size>
std::vector<T> operator+(Array<T,size> const& a, std::vector<T> const& b) {
    ELC_ASSERT(b.size() == size && "size of the array must be equal to the size of the vector");
    std::vector<T> result(size);
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template <typename T, unsigned size>
Array<T,size> operator+(std::vector<T> const& a, Array<T,size> const& b) {
    ELC_ASSERT(a.size() == size && "size of the array must be equal to the size of the vector");
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

template <typename T, unsigned size>
std::vector<T> operator+(std::vector<T> const& a, Array<T,size> const& b) {
    ELC_ASSERT(a.size() == size && "size of the array must be equal to the size of the vector");
    std::vector<T> result(size);
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// =================== subtraction =========================== //

template <typename T, unsigned size>
Array<T,size> operator-(Array<T,size> const& a, std::vector<T> const& b) {
    ELC_ASSERT(b.size() == size && "size of the array must be equal to the size of the vector");
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template <typename T, unsigned size>
std::vector<T> operator-(Array<T,size> const& a, std::vector<T> const& b) {
    ELC_ASSERT(b.size() == size && "size of the array must be equal to the size of the vector");
    std::vector<T> result(size);
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template <typename T, unsigned size>
Array<T,size> operator-(std::vector<T> const& a, Array<T,size> const& b) {
    ELC_ASSERT(a.size() == size && "size of the array must be equal to the size of the vector");
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template <typename T, unsigned size>
std::vector<T> operator-(std::vector<T> const& a, Array<T,size> const& b) {
    ELC_ASSERT(a.size() == size && "size of the array must be equal to the size of the vector");
    std::vector<T> result(size);
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

// =================== multiplication =========================== //

template <typename T, unsigned size>
Array<T,size> operator*(Array<T,size> const& a, std::vector<T> const& b) {
    ELC_ASSERT(b.size() == size && "size of the array must be equal to the size of the vector");
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template <typename T, unsigned size>
std::vector<T> operator*(Array<T,size> const& a, std::vector<T> const& b) {
    ELC_ASSERT(b.size() == size && "size of the array must be equal to the size of the vector");
    std::vector<T> result(size);
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template <typename T, unsigned size>
Array<T,size> operator*(std::vector<T> const& a, Array<T,size> const& b) {
    ELC_ASSERT(a.size() == size && "size of the array must be equal to the size of the vector");
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

template <typename T, unsigned size>
std::vector<T> operator*(std::vector<T> const& a, Array<T,size> const& b) {
    ELC_ASSERT(a.size() == size && "size of the array must be equal to the size of the vector");
    std::vector<T> result(size);
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] * b[i];
    }
    return result;
}

// =================== division =========================== //

template <typename T, unsigned size>
Array<T,size> operator/(Array<T,size> const& a, std::vector<T> const& b) {
    ELC_ASSERT(b.size() == size && "size of the array must be equal to the size of the vector");
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}

template <typename T, unsigned size>
std::vector<T> operator/(Array<T,size> const& a, std::vector<T> const& b) {
    ELC_ASSERT(b.size() == size && "size of the array must be equal to the size of the vector");
    std::vector<T> result(size);
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}

template <typename T, unsigned size>
Array<T,size> operator/(std::vector<T> const& a, Array<T,size> const& b) {
    ELC_ASSERT(a.size() == size && "size of the array must be equal to the size of the vector");
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}

template <typename T, unsigned size>
std::vector<T> operator/(std::vector<T> const& a, Array<T,size> const& b) {
    ELC_ASSERT(a.size() == size && "size of the array must be equal to the size of the vector");
    std::vector<T> result(size);
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] / b[i];
    }
    return result;
}


// =========== scalar-array operations ==================== //

template<typename T, unsigned size>
Array<T,size> operator+(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] + alpha;
    }
    return result;
}

template<typename T, unsigned size>
Array<T,size> operator+(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = alpha + a[i];
    }
    return result;
}


template<typename T, unsigned size>
Array<T,size> operator-(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] - alpha;
    }
    return result;
}

template<typename T, unsigned size>
Array<T,size> operator-(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = alpha - a[i];
    }
    return result;
}

template<typename T, unsigned size>
Array<T,size> operator*(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] * alpha;
    }
    return result;
}

template<typename T, unsigned size>
Array<T,size> operator*(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = alpha * a[i];
    }
    return result;
}



template<typename T, unsigned size>
Array<T,size> operator/(Array<T,size> const& a, T alpha) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = a[i] / alpha;
    }
    return result;
}

template<typename T, unsigned size>
Array<T,size> operator/(T alpha, Array<T,size> const& a) {
    Array<T,size> result;
    for (unsigned i=0; i<size; ++i) {
        result[i] = alpha / a[i];
    }
    return result;
}

// ==== template specific array implementation === //

template <typename T, unsigned size>
struct arrayOperationsImpl {
    static T crossProductNorm(Array<T,size> const& a, Array<T,size> const& b) {
        EPC_ASSERT(false && "No generic implementation for cross product exists.");
        return NAN;
    }

    static T crossProduct(Array<T,size> const& a, Array<T,size> const& b) {
        EPC_ASSERT(false && "No generic implementation for cross product exists.");
        return NAN;
    }

    static void crossProduct(Array<T,size> const& a, Array<T,size> const& b, Array<T,size> &res) {
        EPC_ASSERT(false && "No generic implementation for cross product exists.");
        return Array<T,size>(NAN);
    }

    static Array<T,size> rotateCounterClockwise(Array<T,size> const& a) {
        EPC_ASSERT(false && "no generic implementation existing for rotateCounterClockwise for the moment");
        return Array<T,size>(NAN);
    }

};

template<typename T>
struct arrayOperationsImpl<T,2> {
    static T crossProduct(Array<T,2> const& a, Array<T,2> const& b) {
        return a[0]*b[1]-a[1]*b[0];
    }

    static void crossProduct(Array<T,2> const& a, Array<T,2> const& b, Array<T,2> &res) {
        EPC_ASSERT(false && "no implementation existing for corss product that fills a vector in 2d");
    }

    static T crossProductNorm(Array<T,2> const& a, Array<T,2> const& b) {
        return fabs(crossProduct(a, b));
    }

    static Array<T,2> rotateCounterClockwise(Array<T,2> const& a) {
        Array<T,2> res;
        res[0] = -a[1];
        res[1] = a[0];
        return res;
    }

};

template<typename T>
struct arrayOperationsImpl<T,3> {
    static T crossProduct(Array<T,3> const& a, Array<T,3> const& b) {
        EPC_ASSERT(false && "no implementation existing for corss product that return scalar in 3d");
    }

    static void crossProduct(Array<T,3> const& a, Array<T,3> const& b, Array<T,3> &res) {
        res[0] = a[1]*b[2]-a[2]*b[1];
        res[1] = -(a[0]*b[2]-a[2]*b[0]);
        res[2] = a[0]*b[1]-a[1]*b[0];
    }

    static T crossProductNorm(Array<T,3> const& a, Array<T,3> const& b) {
        Array<T,3> c;
        crossProduct(a, b, c);
        return c.norm();
    }

    static Array<T,3> rotateCounterClockwise(Array<T,3> const& a) {
        EPC_ASSERT(false && "no implementation existing for rotateCounterClockwise in 3d for the moment");
        return Array<T,3>(NAN);
    }
};




} // end namespace epc

#endif

