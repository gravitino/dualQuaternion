#include "dualquat.hpp"
#include <random>

#define NUM_TESTS (1<<22)
#define THRESHOLD (1E-6)

int main () {

    std::mt19937 engine;
    engine.seed(0);
    std::uniform_real_distribution<double> dist(-1, 1);

    for (size_t i = 0; i < NUM_TESTS; i++) {
        quat<double> Q(dist(engine), dist(engine), dist(engine), dist(engine)), R;
        Q = Q.N();

        // check if unit
        if(!Q.isunit()) {
            std::cout << "not unit" << std::endl;
        }

        // check +,-,*,/
        R = (Q-Q*4+Q)/2+Q;
        if (R.dot(R) > THRESHOLD) {
            std::cout << "+,-,*,/ " << R.dot(R) << std::endl;
            R.print();
        }

        // check left inverse
        R = (Q.I()^Q)-quat<double>();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "left inverse" << R.dot(R) << std::endl;
            R.print();
        }

        // check right inverse
        R = (Q^Q.I())-quat<double>();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "right inverse" << R.dot(R) << std::endl;
            R.print();
        }

        // check left conjugate inverse
        R = (Q.C()^Q)-quat<double>();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "left conjugate inverse" << R.dot(R) << std::endl;
            R.print();
        }

        // check right conjugate inverse
        R = (Q^Q.C())-quat<double>();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "right conjugate inverse" << R.dot(R) << std::endl;
            R.print();
        }

        // check log exp
        R = Q.log().exp()-Q;
        if (R.dot(R) > THRESHOLD) {
            std::cout << "log exp " << R.dot(R) << std::endl;
            R.print();
        }

        // check numexp exp
        R = Q.log().numexp()-Q.log().exp();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "numexp exp " << R.dot(R) << std::endl;
        }

        const double eucdist = Q.eucdist(Q*(-1.0)),
                     logdist = Q.logdist(Q*(-1.0));
    
        if (eucdist > THRESHOLD || logdist > THRESHOLD)
            std::cout << "distances: " <<  eucdist << ", " << logdist << std::endl;

    }

    for (size_t i = 0; i < NUM_TESTS; i++) {
        dualquat<double> Q(dist(engine), dist(engine), dist(engine), dist(engine),
                           dist(engine), dist(engine), dist(engine), dist(engine)), R;
        Q = Q.N();

        // check if unit
        if(!Q.isunit()) {
            std::cout << "not unit" << std::endl;
        }

        // check +,-,*,/
        R = (Q-Q*4+Q)/2+Q;
        if (R.dot(R) > THRESHOLD) {
            std::cout << "+,-,*,/ " << R.dot(R) << std::endl;
            R.print();
        }

        // check left inverse
        R = (Q.I()^Q)-dualquat<double>();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "left inverse" << R.dot(R) << std::endl;
            R.print();
        }

        // check right inverse
        R = (Q^Q.I())-dualquat<double>();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "right inverse" << R.dot(R) << std::endl;
            R.print();
        }

        // check left conjugate inverse
        R = (Q.C()^Q)-dualquat<double>();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "left conjugate inverse" << R.dot(R) << std::endl;
            R.print();
        }

        // check right conjugate inverse
        R = (Q^Q.C())-dualquat<double>();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "right conjugate inverse" << R.dot(R) << std::endl;
            R.print();
        }

        // check log exp
        R = Q.log().exp()-Q;
        if (R.dot(R) > THRESHOLD) {
            std::cout << "log exp " << R.dot(R) << std::endl;
            R.print();
        }

        // check numexp exp
        R = Q.log().numexp()-Q.log().exp();
        if (R.dot(R) > THRESHOLD) {
            std::cout << "numexp exp " << R.dot(R) << std::endl;
        }
    
        const double eucdist = Q.eucdist(Q*(-1.0)),
                     logdist = Q.logdist(Q*(-1.0)),
                     baddist = Q.logdist(Q*(-1.0));
    
        if (eucdist > THRESHOLD || logdist > THRESHOLD || baddist > THRESHOLD)
            std::cout << "distances: " <<  eucdist << ", " << logdist << std::endl;

    }

    
}