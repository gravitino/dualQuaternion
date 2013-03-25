#include "quaternion.hpp"

template<class T> quaternion<T>::quaternion(){
    this->w = 1; this->x = 0; this->y = 0; this->z = 0;
}

template<class T> quaternion<T>::quaternion(T w, T x, T y, T z){
    this->w = w; this->x = x; this->y = y; this->z = z;
}

template<class T> void quaternion<T> ::print(){
    
    std::cout << "Quaternion [" 
              << this->w << ", " << this->x << ", "
              << this->y << ", " << this->z << "]" << std::endl;
}

template<class T> quaternion<T> quaternion<T> ::smul(T alpha){
    return quaternion<T>(this->w * alpha, this->x * alpha, 
                         this->y * alpha, this->z * alpha);
}

template<class T> quaternion<T> quaternion<T> ::sdiv(T alpha){
    return quaternion<T>(this->w / alpha, this->x / alpha, 
                         this->y / alpha, this->z / alpha);
}

template<class T> quaternion<T> quaternion<T> ::mul(quaternion<T> q){
    return quaternion<T>(
        this->w*q.w - this->x*q.x - this->y*q.y - this->z*q.z,
        this->w*q.x + this->x*q.w + this->y*q.z - this->z*q.y,
        this->w*q.y - this->x*q.z + this->y*q.w + this->z*q.x,
        this->w*q.z + this->x*q.y - this->y*q.x + this->z*q.w);
}

template<class T> quaternion<T> quaternion<T> ::add(quaternion<T> q){
    return quaternion<T>(this->w + q.w, this->x + q.x, 
                         this->y + q.y, this->z + q.z);
}

template<class T> quaternion<T> quaternion<T> ::sub(quaternion<T> q){
    return quaternion<T>(this->w - q.w, this->x - q.x, 
                         this->y - q.y, this->z - q.z);
}

template<class T> T quaternion<T> ::dot(quaternion<T> q){
    return this->w*q.w + this->x*q.x + this->y*q.y + this->z*q.z;
}

template<class T> quaternion<T> quaternion<T> ::inverse(){
    return this->conjugate().sdiv(this->dot(*this));
}

template<class T> quaternion<T> quaternion<T> ::conjugate(){
    return quaternion<T>(this->w, -(this->x), -(this->y), -(this->z));
}

template<class T> quaternion<T> quaternion<T> ::normalize(){
    return this->sdiv(sqrt(this->dot(*this)));
}

template<class T> quaternion<T> quaternion<T> ::project(){

    if (this->w < 0)
        return this->smul(-1);

    return *this;
}


int main() {

    quaternion<double> Q = quaternion<double> (-1,1, 0, 0);
    Q.mul(Q.inverse()).print();
}

