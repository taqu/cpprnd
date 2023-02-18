#ifndef INC_CPPRND_H_
#define INC_CPPRND_H_
#include <cstdint>

namespace cpprnd
{
#if defined(_MSC_VER_)
#    define CPPRND_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
#    define CPPRND_RESTRICT __restrict
#else
#    define CPPRND_RESTRICT
#endif

//--- PCG32
//---------------------------------------------------------
class PCG32
{
public:
    PCG32();
    explicit PCG32(uint64_t seed);
    ~PCG32();

    void srand(uint64_t seed);

    uint32_t rand();

    float frand();

private:
    inline static constexpr uint64_t Increment = 1442695040888963407ULL;
    inline static constexpr uint64_t Multiplier = 6364136223846793005ULL;
    uint64_t state_;
};

class PCGS32
{
public:
    inline static constexpr uint64_t DefaultStream = 1442695040888963407ULL;
    PCGS32();
    explicit PCGS32(uint64_t seed, uint64_t stream = DefaultStream);
    ~PCGS32();
    void srand(uint64_t seed, uint64_t stream = DefaultStream);
    uint32_t rand();
    float frand();

private:
    PCGS32(const PCGS32&) = delete;
    PCGS32& operator=(const PCGS32&) = delete;
    inline static constexpr uint64_t Multiplier = 6364136223846793005ULL;
    uint64_t increment_;
    uint64_t state_;
};

//--- PCG64
//---------------------------------------------------------
class PCG64
{
public:
    PCG64();
    PCG64(uint64_t seed0, uint64_t seed1);
    ~PCG64();

    void srand(uint64_t seed0, uint64_t seed1);

    uint64_t rand();

    double frand();

private:
    struct UInt128
    {
        uint64_t high_;
        uint64_t low_;
    };
    inline static UInt128 add(const UInt128& x0, const UInt128& x1);
    inline static UInt128 mul(uint64_t x0, uint64_t x1);
    inline static UInt128 mul(const UInt128& x0, const UInt128& x1);
    inline void next();

    inline static constexpr UInt128 Multiplier = {2549297995355413924ULL, 4865540595714422341ULL};
    UInt128 increment_;
    UInt128 state_;
};

//--- SplitMix
//---------------------------------------------------------
class SplitMix
{
public:
    static uint64_t next(uint64_t& state);
};

//--- RandWELL
//---------------------------------------------------------
class RandWELL
{
public:
    RandWELL();
    explicit RandWELL(uint32_t seed);
    ~RandWELL();

    void srand(uint32_t seed);

    uint32_t rand();

    float frand();

private:
    static const uint32_t N = 16;

    uint32_t state_[N];
    uint32_t index_;
};

//--- SFMT
//---------------------------------------------------------
struct w128_t
{
    uint32_t u_[4];
};

struct sfmt_param19937
{
    inline static constexpr uint32_t SFMT_MEXP = 19937;
    inline static constexpr uint32_t SFMT_N = SFMT_MEXP / 128 + 1;
    inline static constexpr uint32_t SFMT_N32 = SFMT_N * 4;

    inline static constexpr uint32_t SFMT_POS1 = 122;
    inline static constexpr uint32_t SFMT_SL1 = 18;
    inline static constexpr uint32_t SFMT_SL2 = 1;
    inline static constexpr uint32_t SFMT_SR1 = 11;
    inline static constexpr uint32_t SFMT_SR2 = 1;
    inline static constexpr uint32_t SFMT_MSK1 = 0xDFFFFFEFU;
    inline static constexpr uint32_t SFMT_MSK2 = 0xDDFECB7FU;
    inline static constexpr uint32_t SFMT_MSK3 = 0xBFFAFFFFU;
    inline static constexpr uint32_t SFMT_MSK4 = 0xBFFFFFF6U;
    inline static constexpr uint32_t SFMT_PARITY1 = 0x00000001UL;
    inline static constexpr uint32_t SFMT_PARITY2 = 0x00000000UL;
    inline static constexpr uint32_t SFMT_PARITY3 = 0x00000000UL;
    inline static constexpr uint32_t SFMT_PARITY4 = 0x13C9E684UL;

    inline static constexpr uint32_t SFMT_ALTI_SL1[] = {SFMT_SL1, SFMT_SL1, SFMT_SL1, SFMT_SL1};
    inline static constexpr uint32_t SFMT_ALTI_SR1[] = {SFMT_SR1, SFMT_SR1, SFMT_SR1, SFMT_SR1};
    inline static constexpr uint32_t SFMT_ALTI_MSK[] = {SFMT_MSK1, SFMT_MSK2, SFMT_MSK3, SFMT_MSK4};
    inline static constexpr uint32_t SFMT_ALTI_MSK64[] = {SFMT_MSK2, SFMT_MSK1, SFMT_MSK4, SFMT_MSK3};
    inline static constexpr uint32_t SFMT_ALTI_SL2_PERM[] = {1, 2, 3, 23, 5, 6, 7, 0, 9, 10, 11, 4, 13, 14, 15, 8};
    inline static constexpr uint32_t SFMT_ALTI_SL2_PERM64[] = {1, 2, 3, 4, 5, 6, 7, 31, 9, 10, 11, 12, 13, 14, 15, 0};
    inline static constexpr uint32_t SFMT_ALTI_SR2_PERM[] = {7, 0, 1, 2, 11, 4, 5, 6, 15, 8, 9, 10, 17, 12, 13, 14};
    inline static constexpr uint32_t SFMT_ALTI_SR2_PERM64[] = {15, 0, 1, 2, 3, 4, 5, 6, 17, 8, 9, 10, 11, 12, 13, 14};

    struct state_t
    {
        w128_t state_[SFMT_N];
        uint32_t index_;
    };
};

class SFMT19937
{
public:
    SFMT19937();
    explicit SFMT19937(uint32_t seed);
    ~SFMT19937();

    void srand(uint32_t seed);
    uint32_t rand();
    float frand();

private:
    sfmt_param19937::state_t state_;
};

//--- MELG
//---------------------------------------------------------
struct melg_param19937
{
    inline static constexpr uint64_t NN = 311;
    inline static constexpr uint64_t MM = 81;
    inline static constexpr uint64_t MATRIX_A = 0x5c32e06df730fc42ULL;
    inline static constexpr uint64_t P = 33;
    inline static constexpr uint64_t W = 64;
    inline static constexpr uint64_t MASKU = (0xffffffffffffffffULL << (W - P));
    inline static constexpr uint64_t MASKL = (~MASKU);
    inline static constexpr uint64_t LAG1 = 19;
    inline static constexpr uint64_t SHIFT1 = 16;
    inline static constexpr uint64_t MASK1 = 0x6AEDE6FD97B338ECULL;
    inline static constexpr uint64_t LAG1over = 292;
    inline static constexpr uint64_t mag01[2] = {0ULL, MATRIX_A};

    struct state_t
    {
        uint64_t lung_;
        uint64_t state_[NN];
        uint32_t index_;
        uint64_t (*function_)(state_t*);
    };
};

class MELG19937
{
public:
    MELG19937();
    explicit MELG19937(uint64_t seed);
    MELG19937(uint64_t size, const uint64_t* seeds);
    ~MELG19937();

    void srand(uint64_t seed);
    void srand(uint64_t size, const uint64_t* seeds);
    uint64_t rand();
    double frand();

private:
    melg_param19937::state_t state_;
};
} // namespace cpprnd
#endif // INC_CPPRND_H_