
#ifndef DUALQUAT_HPP
#define DUALQUAT_HPP

#include <iostream>
#include <limits>
#include <cmath>
#include <assert.h>

#define EPS(value_t) (std::numeric_limits<value_t>::epsilon())

template <class value_t>
struct quat {

    value_t w, x, y, z;

    quat() : w(1), x(0), y(0), z(0) {}

    quat(value_t w_, 
         value_t x_, 
         value_t y_, 
         value_t z_) : w(w_), x(x_) , y(y_), z(z_) {}

    void print () const {
        std::cout << "(" << w << "," << x 
                  << "," << y << "," << z 
                  << ")" << std::endl;
    }

    quat<value_t> operator+(const quat<value_t>& other) const {
        return quat<value_t>(w+other.w, x+other.x,
                             y+other.y, z+other.z);
    }

    quat<value_t>& operator+=(const quat<value_t>& other) {
        w += other.w; x += other.x; y += other.y; z += other.z;
        return *this;
    }

    quat<value_t> operator-(const quat<value_t>& other) const {
        return quat<value_t>(w-other.w, x-other.x,
                             y-other.y, z-other.z);
    }

    quat<value_t>& operator-=(const quat<value_t>& other) {
        w -= other.w; x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    quat<value_t> operator^(const quat<value_t>& other) const {
        return quat<value_t>(w*other.w-x*other.x-y*other.y-z*other.z,
                             w*other.x+other.w*x+y*other.z-z*other.y,
                             w*other.y+other.w*y+z*other.x-x*other.z,
                             w*other.z+other.w*z+x*other.y-y*other.x);
    }

    quat<value_t>& operator^=(const quat<value_t>& other) {
        *this = (*this)^other;
        return *this;
    }

    quat<value_t> operator*(const value_t alpha) const {
        return quat<value_t>(w*alpha, x*alpha,
                             y*alpha, z*alpha);
    }

    quat<value_t> operator/(const value_t alpha) const {
        return quat<value_t>(w/alpha, x/alpha,
                             y/alpha, z/alpha);
    }

    quat<value_t> N() const {
        const value_t rho = sqrt(w*w+x*x+y*y+z*z);
        return quat<value_t>(w/rho, x/rho, y/rho, z/rho);
    }

    quat<value_t> C() const {
        return quat<value_t>(w, -x, -y, -z);
    }

    quat<value_t> I() const {
        const value_t rho2 = w*w+x*x+y*y+z*z;
        return quat<value_t>(w/rho2, -x/rho2, -y/rho2, -z/rho2);
    }

    quat<value_t> exp () const {

        assert(w*w < EPS(value_t));

        value_t dst = std::sqrt(x*x+y*y+z*z);
        value_t fac = std::sin(dst)/dst;        

        return quat<value_t>(std::cos(dst), x*fac, y*fac, z*fac);
    }

    quat<value_t> numexp (value_t eps = EPS(value_t)) const {

        quat<value_t> pow;
        quat<value_t> sum;
        size_t i = 1;

        while (pow.dot(pow) > eps) {
            pow ^= (*this)/i++;
            sum += pow;
        }

        return sum;
    }

    quat<value_t> log () const {

        assert(isunit());

        const value_t inv = 1.0/sqrt(x*x+y*y+z*z+EPS(value_t));
        const value_t fac = std::acos(w > 1 ? 1 : w < -1 ? -1 : w)*inv;

        return quat<value_t> (0, x*fac, y*fac, z*fac);
    }

    value_t dot (const quat<value_t> other) const {
        return w*other.w+x*other.x+y*other.y+z*other.z;
    }

    value_t eucnorm () const {
        return dot(*this);
    }

    value_t lognorm () const {
        const auto& generator = log();
        return generator.dot(generator);
    }

    value_t eucdist (const quat<value_t>& other) const {
        
        if (dot(other) < 0.0) {
            const auto& difference = (*this)+other;
            return difference.dot(difference);  
        }

        const auto& difference = (*this)-other;
        return difference.dot(difference);
    }

    value_t logdist (const quat<value_t>& other) const {
        const auto& difference = ((*this)^(other.C())).log();
        return difference.dot(difference);
    }

    value_t baddist (const quat<value_t>& other) const {

        const bool I_know_what_I_do = false;
        assert(I_know_what_I_do);

        if (dot(other) < 0.0) {
            const auto& difference = log()-(other*(-1.0)).log();
            return difference.dot(difference);  
        }

        const auto& difference = log()-other.log();
        return difference.dot(difference);
    }

    bool isunit () const {
        const value_t residue = w*w+x*x+y*y+z*z-1.0;
        return residue*residue < EPS(value_t);
    }
};

template <class value_t>
struct dualquat {
    
    value_t w, x, y, z, W, X, Y, Z;

    dualquat() : w(1) , x(0), y(0), z(0), W(0), X(0), Y(0), Z(0) {};

    dualquat(value_t w_, 
             value_t x_, 
             value_t y_, 
             value_t z_,
             value_t W_, 
             value_t X_, 
             value_t Y_, 
             value_t Z_) : w(w_), x(x_), y(y_), z(z_), 
                           W(W_), X(X_), Y(Y_), Z(Z_) {}

    dualquat(const quat<value_t>& real, const quat<value_t>& dual) {
        w = real.w; x = real.x; y = real.y; z = real.z;
        W = dual.w; X = dual.x; Y = dual.y; Z = dual.z; 
    }

    quat<value_t> real () const {
        return quat<value_t>(w, x, y, z);
    }

    quat<value_t> dual () const {
        return quat<value_t>(W, X, Y, Z);
    }

    void print () const {
        std::cout << "(" << w << "," << x 
                  << "," << y << "," << z 
                  << "," << W << "," << X
                  << "," << Y << "," << Z
                  << ")" << std::endl;
    }

    dualquat<value_t> operator+(const dualquat<value_t>& other) const {
        return dualquat<value_t>(w+other.w, x+other.x,
                                 y+other.y, z+other.z,
                                 W+other.W, X+other.X,
                                 Y+other.Y, Z+other.Z);
    }

    dualquat<value_t>& operator+=(const dualquat<value_t>& other) {
        w += other.w; x += other.x; y += other.y; z += other.z;
        W += other.W; X += other.X; Y += other.Y; Z += other.Z;
        return *this;
    }

    dualquat<value_t> operator-(const dualquat<value_t>& other) const {
        return dualquat<value_t>(w-other.w, x-other.x,
                                 y-other.y, z-other.z,
                                 W-other.W, X-other.X,
                                 Y-other.Y, Z-other.Z);
    }

    dualquat<value_t>& operator-=(const dualquat<value_t>& other) {
        w -= other.w; x -= other.x; y -= other.y; z -= other.z;
        W -= other.W; X -= other.X; Y -= other.Y; Z -= other.Z;
        return *this;
    }

    // TODO: expand this
    dualquat<value_t> operator^(const dualquat<value_t>& other) const {

        const auto& a = real();
        const auto& A = dual();
        const auto& b = other.real();
        const auto& B = other.dual();

        return dualquat<value_t>(a^b, (a^B)+(A^b));
    }

    dualquat<value_t>& operator^=(const dualquat<value_t>& other) {
        *this = (*this)^other;
        return *this;
    }

    dualquat<value_t> operator*(const value_t alpha) const {
        return dualquat<value_t>(w*alpha, x*alpha,
                                 y*alpha, z*alpha,
                                 W*alpha, X*alpha,
                                 Y*alpha, Z*alpha);
    }

    dualquat<value_t> operator/(const value_t alpha) const {
        return dualquat<value_t>(w/alpha, x/alpha,
                                 y/alpha, z/alpha,
                                 W/alpha, X/alpha,
                                 Y/alpha, Z/alpha);
    }

    dualquat<value_t> N() const {
        
        const value_t qq = w*w+x*x+y*y+z*z+EPS(value_t);
        const value_t qQ = w*W+x*X+y*Y+z*Z;
        const value_t invqq = 1.0/qq;
        const value_t invsq = 1.0/std::sqrt(qq);
        const value_t alpha = qQ*invqq*invqq;

        return dualquat<value_t>(w*invsq, x*invsq, 
                                 y*invsq, z*invsq,
                                 W*invqq-w*alpha, X*invqq-x*alpha,
                                 Y*invqq-y*alpha, Z*invqq-z*alpha);
    }

    dualquat<value_t> C() const {
        return dualquat<value_t>(w, -x, -y, -z, W, -X, -Y, -Z);
    }

    dualquat<value_t> D() const {
        return dualquat<value_t>(w, x, y, z, -W, -X, -Y, -Z);
    }

    dualquat<value_t> I() const {

        const value_t qq = w*w+x*x+y*y+z*z;
        const value_t qQ = w*W+x*X+y*Y+z*Z;
        const value_t invqq = 1.0/qq;
        const value_t alpha = 2.0*qQ*invqq*invqq;

        return dualquat<value_t>(w*invqq, -x*invqq, 
                                -y*invqq, -z*invqq,
                                 W*invqq-w*alpha, -X*invqq-x*alpha,
                                -Y*invqq-y*alpha, -Z*invqq-z*alpha);
    }

    dualquat<value_t> exp() const {

        assert(w*w < EPS(value_t) && W*W < EPS(value_t));

        const value_t theta = 2.0*std::sqrt(x*x+y*y+z*z+EPS(value_t));
        const value_t invvv = 2.0/theta;
        const value_t lx = x*invvv;
        const value_t ly = y*invvv;
        const value_t lz = z*invvv;

        const value_t pitch = 2.0*(lx*X+ly*Y+lz*Z);
        const value_t mx = (2.0*X-pitch*lx)/theta;
        const value_t my = (2.0*Y-pitch*ly)/theta;
        const value_t mz = (2.0*Z-pitch*lz)/theta;

        assert((lx*mx+ly*my+lz*mz)*(lx*mx+ly*my+lz*mz) < EPS(value_t));

        const value_t cost2 = std::cos(0.5*theta);
        const value_t sint2 = std::sin(0.5*theta);
        const value_t alpha = 0.5*pitch*cost2;


        return dualquat<value_t> (cost2, 
                                  sint2*lx,
                                  sint2*ly, 
                                  sint2*lz, 
                                 -0.5*pitch*sint2,
                                  sint2*mx+alpha*lx, 
                                  sint2*my+alpha*ly, 
                                  sint2*mz+alpha*lz);
    }

    dualquat<value_t> numexp (value_t eps = EPS(value_t)) const {

        dualquat<value_t> pow;
        dualquat<value_t> sum;
        size_t i = 1;

        while (pow.dot(pow) > eps) {
            pow ^= (*this)/i++;
            sum += pow;
        }

        return sum;
    }

    dualquat<value_t> log () const {

        assert(isunit());
        
        const value_t theta =  2.0*std::acos(w > 1 ? 1 : w < -1 ? -1 : w);
        const value_t invvv =  0.5/std::sqrt(x*x+y*y+z*z+EPS(value_t));
        const value_t pitch = -4.0*W*invvv;
        const value_t alpha =  pitch*w;

        const value_t lx = x*invvv;
        const value_t ly = y*invvv;
        const value_t lz = z*invvv;
        const value_t mx = (X-lx*alpha)*invvv;
        const value_t my = (Y-ly*alpha)*invvv;
        const value_t mz = (Z-lz*alpha)*invvv;

        assert((lx*mx+ly*my+lz*mz)*(lx*mx+ly*my+lz*mz) < EPS(value_t));

        return dualquat<value_t> (0, theta*lx, 
                                     theta*ly, 
                                     theta*lz, 
                                  0, (pitch*lx+theta*mx), 
                                     (pitch*ly+theta*my), 
                                     (pitch*lz+theta*mz));
    
    }

    value_t dot(const dualquat<value_t> other) const {
        return w*other.w+x*other.x+y*other.y+z*other.z+
               W*other.W+X*other.X+Y*other.Y+Z*other.Z;
    }

    value_t eucnorm () const {
        return dot(*this);
    }

    value_t lognorm () const {
        const auto& generator = log();
        return generator.dot(generator);
    }

    value_t eucdist (const dualquat<value_t>& other) const {
        
        // TODO: could still be errornous
        if (dot(other) < 0.0) {
            const auto& difference = (*this)+other;
            return difference.dot(difference);  
        }

        const auto& difference = (*this)-other;
        return difference.dot(difference);
    }

    value_t logdist (const dualquat<value_t>& other) const {
        const auto& difference = ((*this)^(other.C())).log();
        return difference.dot(difference);
    }

    value_t baddist (const dualquat<value_t>& other) const {

        const bool I_know_what_I_do = false;
        assert(I_know_what_I_do);

        // TODO: could still be errornous
        if (dot(other) < 0.0) {
            const auto& difference = log()-(other*(-1.0)).log();
            return difference.dot(difference);  
        }

        const auto& difference = log()-other.log();
        return difference.dot(difference);
    }

    bool isunit() const {
    
        const value_t residue0 = w*w+x*x+y*y+z*z-1.0;
        const value_t residue1 = w*W+x*X+y*Y+z*Z;
    
        return residue0*residue0 < EPS(value_t) &&
               residue1*residue1 < EPS(value_t);
    }

};

#endif