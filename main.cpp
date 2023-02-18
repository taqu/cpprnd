#include <iostream>
#include <chrono>
#include "Random.h"

int main(void)
{
    using namespace rnd;
    static const uint32_t Samples = 1000;
    {//SFMT
    uint32_t* samples = new uint32_t[Samples];
    SFMT19937 engine;
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    for(uint32_t i=0; i<Samples; ++i) {
        samples[i] = engine.rand();
    }
    std::chrono::high_resolution_clock::duration duration = std::chrono::high_resolution_clock::now() - start;
    uint64_t total = 0;
    for(uint32_t i=0; i<Samples; ++i){
        //std::cout << samples[i] << '\n';
        total += samples[i];
    }
    std::cout.flush();
    std::cout << total << std::endl;
    std::cout << "total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() << "ns" << std::endl;
    std::chrono::high_resolution_clock::duration per1 = duration/Samples;
    std::cout << "/1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(per1).count() << "ns" << std::endl;
    delete[] samples;
    }
    {
        uint64_t* samples = new uint64_t[Samples];
        MELG19937 engine;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < Samples; ++i) {
            samples[i] = engine.rand();
        }
        std::chrono::high_resolution_clock::duration duration = std::chrono::high_resolution_clock::now() - start;
        uint64_t total = 0;
        for(uint32_t i = 0; i < Samples; ++i) {
            std::cout << samples[i] << '\n';
            total += samples[i];
        }
        std::cout.flush();
        std::cout << total << std::endl;
        std::cout << "total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() << "ns" << std::endl;
        std::chrono::high_resolution_clock::duration per1 = duration / Samples;
        std::cout << "/1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(per1).count() << "ns" << std::endl;
        delete[] samples;
    }
    return 0;
}
