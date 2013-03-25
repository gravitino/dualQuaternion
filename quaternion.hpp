#ifndef __QUATERNION_H__
#define __QUATERNION_H__

#include <iostream>
#include <cmath>

template<class T> class quaternion {

    public:
        T w, x, y, z;

        // constructors
        quaternion();                     // default constructor
        quaternion(T, T, T, T);           // custom constructor
        
        // convenience stuff
        void print();                     // print out quaternion to shell
        
        // interaction with scalars and quaternions
        quaternion<T> smul(T);            // multiplication with scalar
        quaternion<T> sdiv(T);            // division by scalar
        quaternion<T> mul(quaternion<T>); // exterior clifford product
        quaternion<T> add(quaternion<T>); // componentwise addition
        quaternion<T> sub(quaternion<T>); // componentwise substraction
        T dot (quaternion<T>);            // inner product
        
        // self interaction
        quaternion<T> conjugate();        // conjugation
        quaternion<T> inverse();          // inverse quaternion
        quaternion<T> normalize();        // project onto unit sphere S^3
        quaternion<T> project();          // project onto upper half plane in R^4

};

#endif
