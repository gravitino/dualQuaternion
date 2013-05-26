#include "quaternion.hpp"

template<class T> quat<T>::quat(){
    this->w = 1; this->x = 0; this->y = 0; this->z = 0;
}

template<class T> quat<T>::quat(T w, T x, T y, T z){
    this->w = w; this->x = x; this->y = y; this->z = z;
}

template<class T> void quat<T> ::print(){
    
    std::cout << "quat [" 
              << w << ", " << x << ", "
              << y << ", " << z << "]" << std::endl;
}

template<class T> quat<T> quat<T> ::smul(T alpha){
    return quat<T>(w*alpha, x*alpha, y*alpha, z*alpha);
}

template<class T> quat<T> quat<T> ::sdiv(T alpha){
    return quat<T>(w/alpha, x/alpha, y/alpha, z/alpha);
}

template<class T> quat<T> quat<T> ::mul(quat<T> q){
    return quat<T>(
        w*q.w - x*q.x - y*q.y - z*q.z,
        w*q.x + x*q.w + y*q.z - z*q.y,
        w*q.y - x*q.z + y*q.w + z*q.x,
        w*q.z + x*q.y - y*q.x + z*q.w);
}

template<class T> quat<T> quat<T> ::add(quat<T> q){
    return quat<T>(w+q.w, x+q.x, y+q.y, z+q.z);
}

template<class T> quat<T> quat<T> ::sub(quat<T> q){
    return quat<T>(w-q.w, x-q.x, y-q.y, z-q.z);
}

template<class T> T quat<T> ::dot(quat<T> q){
    return w*q.w + x*q.x + y*q.y + z*q.z;
}

template<class T> quat<T> quat<T> ::inverse(){
    return conjugate().sdiv(dot(*this));
}

template<class T> quat<T> quat<T> ::conjugate(){
    return quat<T>(w, -x, -y, -z);
}

template<class T> quat<T> quat<T> ::normalize(){
    return sdiv(sqrt(dot(*this)));
}

template<class T> quat<T> quat<T> ::project(){

    if (w < 0)
        return smul(-1);

    return *this;
}

template<class T> quat<T> quat<T> ::log(){

    T qnorm = sqrt(dot(*this));
    T vnorm = sqrt(x*x + y*y + z*z);
    T angle = acos(w/qnorm);
      angle = 0 == angle ? 0 : angle/vnorm;

    return quat<T> (::log(qnorm), x*angle, y*angle, z*angle);
}

template<class T> quat<T> quat<T> ::exp(T prec=0){

    quat<T> sum = quat<T>();
    quat<T> pow = *this;
    unsigned int i = 1;
    
    while (pow.dot(pow) > prec) {
        
        sum = sum.add(pow);
        pow = pow.mul(*this).sdiv(++i);
    }
    
    return sum;
}

template<class T> dualQuat<T>::dualQuat(){
    
    q = quat<T>(1, 0, 0, 0);
    Q = quat<T>(0, 0, 0, 0);
}

template<class T> dualQuat<T>::dualQuat(T w, T x, T y, T z, T W, T X, T Y ,T Z){
    q = quat<T>(w, x, y, z);
    Q = quat<T>(W, X, Y, Z);
}

template<class T> dualQuat<T>::dualQuat(quat<T> q, quat<T> Q){
    this->q = q;
    this->Q = Q;
}

template<class T> void dualQuat<T> ::print(){
    
    std::cout << "dualQuat [" 
              << q.w << ", " << q.x << ", "
              << q.y << ", " << q.z << ", " 
              << Q.w << ", " << Q.x << ", "
              << Q.y << ", " << Q.z << "]" << std::endl;
}

template<class T> dualQuat<T> dualQuat<T> ::smul(T alpha){
    return dualQuat<T>(q.smul(alpha), Q.smul(alpha));
}

template<class T> dualQuat<T> dualQuat<T> ::sdiv(T alpha){
    return dualQuat<T>(q.sdiv(alpha), Q.sdiv(alpha));
}

template<class T> dualQuat<T> dualQuat<T> ::mul(dualQuat<T> d){
    return dualQuat<T>(q.mul(d.q), q.mul(d.Q).add(Q.mul(d.q)));
}

template<class T> dualQuat<T> dualQuat<T> ::add(dualQuat<T> d){
    return dualQuat<T>(q.add(d.q), Q.add(d.Q));
}

template<class T> dualQuat<T> dualQuat<T> ::sub(dualQuat<T> d){
    return dualQuat<T>(q.sub(d.q), Q.sub(d.Q));
}

template<class T> T dualQuat<T> ::dot(dualQuat<T> d){
    return q.dot(d.q) + Q.dot(d.Q);
}

template<class T> dualQuat<T> dualQuat<T> ::inverse(){
    
    T qq = q.dot(q);
    T qQ = q.dot(Q);
    
    quat<T> q0 = q.conjugate().sdiv(qq);
    quat<T> qe = Q.conjugate().sdiv(qq).sub(q0.smul(2*qQ).sdiv(qq));
       
    return dualQuat<T>(q0, qe);
}

template<class T> dualQuat<T> dualQuat<T> ::fullConjugate(){
    return dualQuat<T>(q.conjugate(), Q.conjugate().smul(-1));
}

template<class T> dualQuat<T> dualQuat<T> ::compConjugate(){
    return dualQuat<T>(q.conjugate(), Q.conjugate());
}

template<class T> dualQuat<T> dualQuat<T> ::dualConjugate(){
    return dualQuat<T>(q, Q.smul(-1));
}

template<class T> dualQuat<T> dualQuat<T> ::normalize(){
    
    T qq = q.dot(q);
    T sq = sqrt(qq);
    T qQ = q.dot(Q);
    
    quat<T> q0 = q.sdiv(sq);
    quat<T> qe = Q.sdiv(qq).sub(q0.smul(qQ).sdiv(qq*sq));
       
    return dualQuat<T>(q0, qe);
}

template<class T> dualQuat<T> dualQuat<T> ::project(){

    if (q.w < 0)
        return *this.smul(-1);

    return *this;
}

template<class T> dualQuat<T> dualQuat<T> ::log(T prec=0){

    dualQuat<T> sum = dualQuat<T>(0, 0, 0, 0, 0, 0, 0, 0);
    dualQuat<T> pow = (*this).sub(dualQuat<T>());
    dualQuat<T> fac = (*this).sub(dualQuat<T>()).smul(-1);
    unsigned int i = 0;
    
    while (pow.dot(pow) > prec) {
        
        sum = sum.add(pow.sdiv(++i));
        pow = pow.mul(fac);
    }
    
    return sum;
}

template<class T> dualQuat<T> dualQuat<T> ::exp(T prec=0){

    dualQuat<T> sum = dualQuat<T>();
    dualQuat<T> pow = *this;
    unsigned int i = 1;
    
    while (pow.dot(pow) > prec) {
        
        sum = sum.add(pow);
        pow = pow.mul(*this).sdiv(++i);
    }
    
    return sum;
}

int main() {

    // vector x, rotation R, translation T
    dualQuat<double> x = dualQuat<double>(1, 0, 0, 0, 0, 1, 1, 0);
    dualQuat<double> R = dualQuat<double>(1.0/sqrt(2), 0, 0, 1.0/sqrt(2), 0, 0, 0, 0);
    dualQuat<double> T = dualQuat<double>(1, 0, 0, 0, 0, 0.5*100, 0.5*100, 0.5*100);
    dualQuat<double> Q = T.mul(R);
    
    // R_z(90°) * x + (100, 100, 100)
    std::cout << "vector x: "; x.print();
    std::cout << "R*x+T:    "; Q.mul(x).mul(Q.fullConjugate()).print();
    std::cout << std::endl;
    
    // inverse
    std::cout << "Q*Q^(-1):        "; Q.mul(Q.inverse()).print();
    std::cout << "Qnorm*Qnorm^(*): "; Q.normalize().mul(Q.normalize().compConjugate()).print();
    std::cout << std::endl;
    
    // logarithm is broken for dualQuat but works for quat
    std::cout << "ordinary quat q: "; Q.q.print();
    std::cout << "exp(log(q))=q:   "; Q.q.log().exp().print();
    std::cout << std::endl;
    
    std::cout << "dual quat Q:     "; Q.print();
    std::cout << "exp(log(Q)) = Q: "; Q.log().exp().print();
    std::cout << std::endl;
}

