#ifndef INC_CPPRND_H_
#define INC_CPPRND_H_
/**
 * @mainpage
 * # cpprnd
 *
 */
#include <cassert>
#include <cstdint>
#include <utility>

namespace cpprnd
{
#if defined(_MSC_VER_)
#    define CPPRND_RESTRICT __restrict
#elif defined(__GNUC__) || defined(__clang__)
#    define CPPRND_RESTRICT __restrict
#else
#    define CPPRND_RESTRICT
#endif

inline static constexpr uint32_t CPPRNG_DEFAULT_SEED32 = 12345UL;
inline static constexpr uint64_t CPPRNG_DEFAULT_SEED64 = 12345ULL;

//--- PCG32
//---------------------------------------------------------
/**
 * @brief A fast 32 bit PRNG
 *
 * | Feature |      |
 * | :------ | :--- |
 * | Bits    | 32   |
 * | Period  | 2^64 |
 * | Streams | 1    |
 */
class PCG32
{
public:
    /**
     * @brief Initialize with CPPRNG_DEFAULT_SEED64
     */
    PCG32();

    /**
     * @brief Initialize with a seed
     * @param [in] seed ... initialize states with
     */
    explicit PCG32(uint64_t seed);
    ~PCG32();

    /**
     * @brief Initialize states with a seed
     * @param [in] seed
     */
    void srand(uint64_t seed);

    /**
     * @brief Generate a 32bit unsigned value
     * @return
     */
    uint32_t rand();

    /**
     * @brief Generate a 32bit real number
     * @return [0 1)
     */
    float frand();

private:
    inline static constexpr uint64_t Increment = 1442695040888963407ULL;
    inline static constexpr uint64_t Multiplier = 6364136223846793005ULL;
    uint64_t state_;
};

/**
 * @brief A fast 32 bit PRNG
 *
 * | Feature |      |
 * | :------ | :--- |
 * | Bits    | 32   |
 * | Period  | 2^64 |
 * | Streams | 2^63 |
 */
class PCGS32
{
public:
    inline static constexpr uint64_t DefaultStream = 1442695040888963407ULL;
    /**
     * @brief Initialize with CPPRNG_DEFAULT_SEED64 and DefaultStream
     */
    PCGS32();

    /**
     * @brief Initialize with a seed and a stream
     * @param [in] seed ... initialize states with
     * @param [in] stream
     */
    explicit PCGS32(uint64_t seed, uint64_t stream = DefaultStream);
    ~PCGS32();

    /**
     * @brief Initialize with a seed and a stream
     * @param [in] seed ... initialize states with
     * @param [in] stream
     */
    void srand(uint64_t seed, uint64_t stream = DefaultStream);

    /**
     * @brief Generate a 32bit unsigned value
     * @return
     */
    uint32_t rand();

    /**
     * @brief Generate a 32bit real number
     * @return [0 1)
     */
    float frand();

private:
    inline static constexpr uint64_t Multiplier = 6364136223846793005ULL;
    uint64_t increment_;
    uint64_t state_;
};

//--- PCG64
//---------------------------------------------------------
/**
 * @brief A 64 bit PRNG
 *
 * | Feature |      |
 * | :------ | :--- |
 * | Bits    | 64   |
 * | Period  | 2^128|
 * | Streams | 1    |
 */
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
/**
 * @brief A fast 64 bit PRNG
 *
 * | Feature |      |
 * | :------ | :--- |
 * | Bits    | 64   |
 * | Period  | 2^64 |
 * | Streams | 1    |
 */
class SplitMix
{
public:
    static uint64_t next(uint64_t& state);
};

//--- RandWELL512
//---------------------------------------------------------
/**
 * @brief A high-dimensional uniform distributed 32 bit PRNG
 *
 * | Feature |       |
 * | :------ | :---- |
 * | Bits    | 32    |
 * | Period  | 2^512 |
 * | Streams | 1     |
 */
class RandWELL512
{
public:
    /**
     * @brief Initialize with CPPRNG_DEFAULT_SEED32
     */
    RandWELL512();
    /**
     * @brief Initialize with a seed
     * @param [in] seed ... initialize states with
     */
    explicit RandWELL512(uint32_t seed);
    ~RandWELL512();

    /**
     * @brief Initialize with a seed
     * @param [in] seed ... initialize states with
     */
    void srand(uint32_t seed);

    /**
     * @brief Generate a 32bit unsigned value
     * @return
     */
    uint32_t rand();

    /**
     * @brief Generate a 32bit real number
     * @return [0 1)
     */
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

/**
 * @brief A long term 32 bit PRNG
 *
 * | Feature |         |
 * | :------ | :------ |
 * | Bits    | 32      |
 * | Period  | 2^19937 |
 * | Streams | 1       |
 */
class SFMT19937
{
public:
    /**
     * @brief Initialize with CPPRNG_DEFAULT_SEED32
     */
    SFMT19937();
    /**
     * @brief Initialize with a seed
     * @param [in] seed ... initialize states with
     */
    explicit SFMT19937(uint32_t seed);
    ~SFMT19937();

    /**
     * @brief Initialize with a seed
     * @param [in] seed ... initialize states with
     */
    void srand(uint32_t seed);

    /**
     * @brief Generate a 32bit unsigned value
     * @return
     */
    uint32_t rand();

    /**
     * @brief Generate a 32bit real number
     * @return [0 1)
     */
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

/**
 * @brief A high-dimensional uniform distributed 64 bit PRNG
 *
 * | Feature |         |
 * | :------ | :------ |
 * | Bits    | 64      |
 * | Period  | 2^19937 |
 * | Streams | 1       |
 */
class MELG19937
{
public:
    /**
     * @brief Initialize with CPPRNG_DEFAULT_SEED64
     */
    MELG19937();
    /**
     * @brief Initialize with a seed
     * @param [in] seed ... initialize states with
     */
    explicit MELG19937(uint64_t seed);
    /**
     * @brief Initialize with seeds
     * @param [in] size ...  the number of seeds
     * @param [in] seeds ... seeds for initialization
     */
    MELG19937(uint64_t size, const uint64_t* seeds);
    ~MELG19937();

    /**
     * @brief Initialize with a seed
     * @param [in] seed ... initialize states with
     */
    void srand(uint64_t seed);

    /**
     * @brief Initialize with seeds
     * @param [in] size ...  the number of seeds
     * @param [in] seeds ... seeds for initialization
     */
    void srand(uint64_t size, const uint64_t* seeds);

    /**
     * @brief Generate a 64bit unsigned value
     * @return
     */
    uint64_t rand();

    /**
     * @brief Generate a 64bit real number
     * @return [0 1)
     */
    double frand();

private:
    melg_param19937::state_t state_;
};

#ifdef CPPRND_CRPNG
//--- CPRNG
//---------------------------------------------------------
uint32_t crypt_rand32();
uint64_t crypt_rand64();
bool crypt_rand(uint32_t size, void* buffer);
#endif

//--- Common algorithms
//---------------------------------------------------------
/**
 * @brief Return [vmin, vmax)
 */
template<class T, class U>
U range_ropen(T& random, U vmin, U vmax)
{
    assert(vmin <= vmax);
    return static_cast<U>(random.frand() * (vmax - vmin)) + vmin;
}

/**
 * @brief [0, v)
 */
template<class T, class U>
U range_ropen(T& random, U v)
{
    assert(0 <= v);
    return static_cast<U>(random.frand() * v);
}

/**
 * @brief return [0, s)
 * [Fast Random Integer Generation in an Interval](https://arxiv.org/abs/1805.10941)
 */
template<class T>
uint32_t range_ropen(T& random, uint32_t s)
{
    uint32_t x = random.rand();
    uint64_t m = static_cast<uint64_t>(x) * static_cast<uint64_t>(s);
    uint32_t l = static_cast<uint32_t>(m);
    if(l < s) {
        uint32_t t = static_cast<uint32_t>(-static_cast<int32_t>(s) % s);
        while(l < t) {
            x = random.rand();
            m = static_cast<uint64_t>(x) * static_cast<uint64_t>(s);
            l = static_cast<uint32_t>(m);
        }
    }
    return static_cast<uint32_t>(m >> 32);
}

/**
 * @brief return [0, s)
 */
template<class T>
int32_t range_ropen(T& random, int32_t s)
{
    assert(0 <= s);
    return static_cast<int32_t>(range_ropen(random, static_cast<uint32_t>(s)));
}

/**
 * @brief Return [vmin, vmax)
 */
template<class T>
uint32_t range_ropen(T& random, uint32_t vmin, uint32_t vmax)
{
    assert(vmin <= vmax);
    return range_ropen(random, vmax - vmin) + vmin;
}

/**
 * @brief return [vim, vmax)
 */
template<class T>
int32_t range_ropen(T& random, int32_t vmin, int32_t vmax)
{
    assert(vmin <= vmax);
    return range_ropen(random, vmax - vmin) + vmin;
}

/**
 * @brief return [0, v)
 * [Fast Random Integer Generation in an Interval](https://arxiv.org/abs/1805.10941)
 */
template<class T>
uint64_t range_ropen(T& random, uint64_t s)
{
    assert(0 <= s);
    uint64_t t = (~s + 1) % s;
    uint64_t x;
    do {
        x = random.rand();
    } while(x < t);
    return x % s;
}

template<class T>
int64_t range_ropen(T& random, int64_t s)
{
    assert(0 <= s);
    return static_cast<int64_t>(range_ropen(random, static_cast<uint64_t>(s)));
}

/**
 * @brief return [vmin, vmax)
 */
template<class T>
uint64_t range_ropen(T& random, uint64_t vmin, uint64_t vmax)
{
    assert(vmin <= vmax);
    return range_ropen(random, vmax - vmin) + vmin;
}

template<class T>
int64_t range_ropen(T& random, int64_t vmin, int64_t vmax)
{
    assert(vmin <= vmax);
    return range_ropen(random, vmax - vmin) + vmin;
}

template<class T, class U>
void shuffle(T& random, U* start, U* end)
{
    std::ptrdiff_t size = static_cast<std::ptrdiff_t>(end - start);
    for(; 1 < size; --size) {
        std::size_t offset = range_ropen(random, size);
        std::swap(*(start + size - 1), *(start + offset));
    }
}

template<class T, class U>
void shuffle(T& random, std::size_t num, U* v)
{
    for(std::size_t i = num; 1 < i; --i) {
        std::size_t offset = range_ropen(random, i);
        std::swap(v[i - 1], v[offset]);
    }
}

//--- RandomBinarySelect
//---------------------------------------------------------
class RandomBinarySelect
{
public:
    RandomBinarySelect();
    ~RandomBinarySelect();
    uint32_t size() const;
    void build(uint32_t size, const float* weights);
    template<class T>
    uint32_t select(T& random) const;

private:
    RandomBinarySelect(const RandomBinarySelect&) = delete;
    RandomBinarySelect& operator=(const RandomBinarySelect&) = delete;
    uint32_t capacity_;
    uint32_t size_;
    float* weights_;
};

template<class T>
uint32_t RandomBinarySelect::select(T& random) const
{
    const float* first = weights_;
    float x = random.frand();
    int32_t count = size_;
    while(0 < count) {
        int32_t d = count >> 1;
        const float* m = first + d;
        if(*m < x) {
            first = m + 1;
            count -= d + 1;
        } else {
            count = d;
        }
    }
    uint32_t result = static_cast<uint32_t>(first - weights_);
    return size_ <= result ? size_ - 1 : result;
}

//--- RandomAliasSelect
//---------------------------------------------------------
class RandomAliasSelect
{
public:
    RandomAliasSelect();
    ~RandomAliasSelect();
    uint32_t size() const;
    void build(uint32_t size, float* weights);
    template<class T>
    uint32_t select(T& random) const;

private:
    RandomAliasSelect(const RandomAliasSelect&) = delete;
    RandomAliasSelect& operator=(const RandomAliasSelect&) = delete;
    uint32_t capacity_;
    uint32_t size_;
    float* weights_;
    uint32_t* aliases_;
};

template<class T>
uint32_t RandomAliasSelect::select(T& random) const
{
    uint32_t index = range_ropen(random, size_);
    float w = random.frand();
    return w < weights_[index] ? index : aliases_[index];
}
} // namespace cpprnd
#endif // INC_CPPRND_H_
