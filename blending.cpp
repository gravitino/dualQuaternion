#include "dualquat.hpp"
#include <random>

int main () {

    using namespace average;

    std::mt19937 engine;
    engine.seed(0);
    std::uniform_real_distribution<double> dist(-1, 1);
    
    size_t N = 10000;
    
    auto q0 = quat<double>(0, 1, 1, 1).exp(); q0.print();

    std::vector<double> weights(N, 1.0/N);
    std::vector<quat<double>> quats;

    for (size_t i = 0; i < N; i++)
        quats.push_back( q0*quat<double>(1+dist(engine), dist(engine), 
                                           dist(engine), dist(engine)).N());

    QLA(quats).print();
    std::cout << QLA(quats).logdist(q0) << std::endl;
    QLB(quats, weights).print();
    std::cout << QLB(quats, weights).logdist(q0) << std::endl;
    QIA(quats).print();
    std::cout << QIA(quats).logdist(q0) << std::endl;
    QIB(quats, weights).print();
    std::cout << QIB(quats, weights).logdist(q0) << std::endl;

    std::cout << std::endl;

    auto Q0 = dualquat<double>(0, 1, 1, 1, 0, 1, 2, 3).exp(); Q0.print();
    std::vector<double> Weights(N, 1.0/N);
    std::vector<dualquat<double>> Quats;

    for (size_t i = 0; i < N; i++)
        Quats.push_back(Q0*dualquat<double>(1+dist(engine), dist(engine), 
                                              dist(engine), dist(engine),
                                              dist(engine), dist(engine),
                                              dist(engine), dist(engine)).N());

    DLA(Quats).print();
    std::cout << DLA(Quats).logdist(Q0) << std::endl;
    DLB(Quats, Weights).print();
    std::cout << DLB(Quats, Weights).logdist(Q0) << std::endl;
    DIA(Quats).print();
    std::cout << DIA(Quats).logdist(Q0) << std::endl;
    DIB(Quats, Weights).print();
    std::cout << DIB(Quats, Weights).logdist(Q0) << std::endl;

}