#ifndef __quat_H__
#define __quat_H__

#include <iostream>
#include <cmath>

template<class T> class quat {

    public:
        T w, x, y, z;

        // constructors
        quat();                     // default constructor
        quat(T, T, T, T);           // custom constructor
        
        // convenience stuff
        void print();                     // print out quat
        
        // interaction with scalars and quats
        quat<T> smul(T);            // multiplication with scalar
        quat<T> sdiv(T);            // division by scalar
        quat<T> mul(quat<T>); // exterior clifford product
        quat<T> add(quat<T>); // componentwise addition
        quat<T> sub(quat<T>); // componentwise substraction
        T dot (quat<T>);            // inner product
        
        // self interaction
        quat<T> conjugate();        // conjugation
        quat<T> inverse();          // inverse quat
        quat<T> normalize();        // project onto unit sphere
        quat<T> project();          // project onto upper plane
        quat<T> exp(T);             // exponential map
        quat<T> log();              // logarithm map

};


template<class T> class dualQuat {

    public:
        quat<T> q, Q;

        // constructors
        dualQuat();                         // default constructor
        dualQuat(T, T, T, T, T, T, T, T);   // custom constructor
        dualQuat(quat<T>, quat<T>);
        
        // convenience stuff
        void print();                             // print out quat
        
        // interaction with scalars and quats
        dualQuat<T> smul(T);                // multiplication with scalar
        dualQuat<T> sdiv(T);                // division by scalar
        dualQuat<T> mul(dualQuat<T>);       // exterior clifford product
        dualQuat<T> add(dualQuat<T>);       // componentwise addition
        dualQuat<T> sub(dualQuat<T>);       // componentwise substraction
        T dot (dualQuat<T>);                // inner product
        
        // self interaction
        dualQuat<T> fullConjugate();        // dual and complex conjugation
        dualQuat<T> compConjugate();        // complex conjugation
        dualQuat<T> dualConjugate();        // dual conjugation
        dualQuat<T> inverse();              // inverse quat
        dualQuat<T> normalize();            // project onto unit sphere
        dualQuat<T> project();              // project onto upper plane
        dualQuat<T> exp(T);                 // exponential map
        dualQuat<T> log(T);                 // logarithm map

};
#endif
