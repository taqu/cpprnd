#include "cpprnd.h"
#include <chrono>
#include <iostream>

int main(void)
{
    using namespace cpprnd;
    static const uint32_t Samples = 1000;
    {//PCG32
        uint32_t* samples = new uint32_t[Samples];
        PCG32 engine;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < Samples; ++i) {
            samples[i] = engine.rand();
        }
        std::chrono::high_resolution_clock::duration duration = std::chrono::high_resolution_clock::now() - start;
        uint64_t total = 0;
        for(uint32_t i = 0; i < Samples; ++i) {
            // std::cout << samples[i] << '\n';
            total += samples[i];
        }
        std::cout.flush();
        std::cout << "PCG32: " << total << std::endl;
        std::cout << "total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() << "ns" << std::endl;
        std::chrono::high_resolution_clock::duration per1 = duration / Samples;
        std::cout << "/1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(per1).count() << "ns" << std::endl;
        delete[] samples;
    }
    {//PCG64
        uint64_t* samples = new uint64_t[Samples];
        PCG64 engine;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < Samples; ++i) {
            samples[i] = engine.rand();
        }
        std::chrono::high_resolution_clock::duration duration = std::chrono::high_resolution_clock::now() - start;
        uint64_t total = 0;
        for(uint32_t i = 0; i < Samples; ++i) {
            //std::cout << samples[i] << '\n';
            total += samples[i];
        }
        std::cout.flush();
        std::cout << "PCG64: " << total << std::endl;
        std::cout << "total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() << "ns" << std::endl;
        std::chrono::high_resolution_clock::duration per1 = duration / Samples;
        std::cout << "/1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(per1).count() << "ns" << std::endl;
        delete[] samples;
    }
    {//SFMT
        uint32_t* samples = new uint32_t[Samples];
        SFMT19937 engine(1234UL);
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < Samples; ++i) {
            samples[i] = engine.rand();
        }
        std::chrono::high_resolution_clock::duration duration = std::chrono::high_resolution_clock::now() - start;
        uint64_t total = 0;
        for(uint32_t i = 0; i < Samples; ++i) {
            // std::cout << samples[i] << '\n';
            total += samples[i];
        }
        std::cout.flush();
        std::cout << "SFMT: " << total << std::endl;
        std::cout << "total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() << "ns" << std::endl;
        std::chrono::high_resolution_clock::duration per1 = duration / Samples;
        std::cout << "/1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(per1).count() << "ns" << std::endl;
        delete[] samples;
    }
    {//MELG
        uint64_t init[4]={0x12345ULL, 0x23456ULL, 0x34567ULL, 0x45678ULL};
        uint64_t* samples = new uint64_t[Samples];
        MELG19937 engine(4, init);
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        for(uint32_t i = 0; i < Samples; ++i) {
            samples[i] = engine.rand();
        }
        std::chrono::high_resolution_clock::duration duration = std::chrono::high_resolution_clock::now() - start;
        uint64_t total = 0;
        for(uint32_t i = 0; i < Samples; ++i) {
            //std::cout << samples[i] << '\n';
            total += samples[i];
        }
        std::cout.flush();
        std::cout << "MELG: " << total << std::endl;
        std::cout << "total: " << std::chrono::duration_cast<std::chrono::nanoseconds>(duration).count() << "ns" << std::endl;
        std::chrono::high_resolution_clock::duration per1 = duration / Samples;
        std::cout << "/1: " << std::chrono::duration_cast<std::chrono::nanoseconds>(per1).count() << "ns" << std::endl;
        delete[] samples;
    }
    return 0;
}
