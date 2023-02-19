#include "cpprnd.h"
#if defined(__linux__) || defined(__unix__)
#    ifdef CPPRND_CRPNG
#        include <fcntl.h>
#        include <sys/types.h>
#        include <unistd.h>
#    endif
#endif
#if defined(_WIN32)
#    include <Windows.h>
#endif

#if defined(SFMT_SSE)
#    include <immintrin.h>
#endif

namespace cpprnd
{
namespace
{
    /**
     * @brief 32 bit right rotation
     * @param [in] x ... input
     * @param [in] r ... count of rotation
     * @return rotated
     */
    inline uint32_t rotr32(uint32_t x, uint32_t r)
    {
        return (x >> r) | (x << ((~r + 1) & 31U));
    }

    /**
     * @brief 64 bit right rotation
     * @param [in] x ... input
     * @param [in] r ... count of rotation
     * @return rotated
     */
    inline uint64_t rotr64(uint64_t value, uint32_t rot)
    {
        return (value >> rot) | (value << ((~rot + 1) & 63U));
    }

    /**
     * @brief Convert to a [0 1) real number
     * @param [in] x
     * @return a [0 1) real number
     */
    inline float to_real32(uint32_t x)
    {
        return static_cast<float>((x >> 8) * (1.0 / 16777216.0));
    }

    /**
     * @brief Convert to a [0 1) real number
     * @param [in] x
     * @return a [0 1) real number
     */
    inline double to_real64(uint64_t x)
    {
        return (x >> 11) * (1.0 / 9007199254740992.0);
    }

    /**
     * @brief Scramble an input
     * @param [in] x
     * @return scrambled
     */
    uint32_t scramble(uint32_t x)
    {
        x += 0x7f4A7C15UL;
        uint32_t t = x;
        t = (t ^ (t >> 15)) * 0x1CE4E5B9UL;
        t = (t ^ (t >> 13)) * 0x133111EBUL;
        return t ^ (t >> 15);
    }

    /**
     * @brief Scramble an input
     * @param [in] x
     * @return scrambled
     */
    uint64_t scramble(uint64_t x)
    {
        x += 0x9E3779B97f4A7C15ULL;
        uint64_t t = x;
        t = (t ^ (t >> 30)) * 0xBF58476D1CE4E5B9ULL;
        t = (t ^ (t >> 27)) * 0x94D049BB133111EBULL;
        return t ^ (t >> 31);
    }
} // namespace

//--- PCG32
//------------------------------------------------------------
PCG32::PCG32()
    : state_{CPPRNG_DEFAULT_SEED64}
{
}

PCG32::PCG32(uint64_t seed)
{
    srand(seed);
}

PCG32::~PCG32()
{
}

void PCG32::srand(uint64_t seed)
{
    state_ = SplitMix::next(seed);
    while(0 == state_) {
        state_ = SplitMix::next(state_);
    }
}

uint32_t PCG32::rand()
{
    uint64_t x = state_;
    uint32_t c = static_cast<uint32_t>(x >> 59);
    state_ = x * Multiplier + Increment;
    x ^= x >> 18;
    return rotr32(static_cast<uint32_t>(x >> 27), c);
}

float PCG32::frand()
{
    return to_real32(rand());
}

PCGS32::PCGS32()
    : state_(CPPRNG_DEFAULT_SEED64)
    , increment_(DefaultStream)
{
}

PCGS32::PCGS32(uint64_t seed, uint64_t stream)
{
    srand(seed);
    increment_ = stream | 1ULL;
}

PCGS32::~PCGS32()
{
}

void PCGS32::srand(uint64_t seed, uint64_t stream)
{
    state_ = SplitMix::next(seed);
    while(0 == state_) {
        state_ = SplitMix::next(state_);
    }
    increment_ = stream | 1ULL;
}

uint32_t PCGS32::rand()
{
    uint64_t x = state_;
    uint32_t c = static_cast<uint32_t>(x >> 59);
    state_ = x * Multiplier + increment_;
    x ^= x >> 18;
    return rotr32(static_cast<uint32_t>(x >> 27), c);
}

float PCGS32::frand()
{
    return to_real32(rand());
}

//--- PCG64
//------------------------------------------------------------
inline PCG64::UInt128 PCG64::add(const UInt128& x0, const UInt128& x1)
{
    uint64_t low = x0.low_ + x1.low_;
    uint64_t high = x0.high_ + x1.high_ + (low < x1.low_);
    return {high, low};
}

inline PCG64::UInt128 PCG64::mul(uint64_t x0, uint64_t x1)
{
    uint64_t low = x0 * x1;

    uint64_t a0 = x0 & 0xFFFFFFFFULL;
    uint64_t a1 = x0 >> 32U;
    uint64_t b0 = x1 & 0xFFFFFFFFULL;
    uint64_t b1 = x1 >> 32U;
    uint64_t c0 = a0 * b0;
    uint64_t t = a1 * b0 + (c0 >> 32U);
    uint64_t c1 = t & 0xFFFFFFFFULL;
    uint64_t c2 = t >> 32;
    c1 += a0 * b1;
    uint64_t high = a1 * b1 + c2 + (c1 >> 32U);
    return {high, low};
}

inline PCG64::UInt128 PCG64::mul(const UInt128& x0, const UInt128& x1)
{
    uint64_t high = x0.high_ * x1.low_ + x0.low_ * x1.high_;
    UInt128 result = mul(x0.low_, x1.low_);
    result.high_ += high;
    return result;
}

inline void PCG64::next()
{
    state_ = add(mul(state_, Multiplier), increment_);
}

PCG64::PCG64()
    : increment_{6364136223846793005ULL, 1442695040888963407ULL}
    , state_{0x4D595DF4D0F33173ULL, 0xDA3E39CB94B95BDBULL}
{
}

PCG64::PCG64(uint64_t seed0, uint64_t seed1)
    : increment_{6364136223846793005ULL, 1442695040888963407ULL}
{
    srand(seed0, seed1);
}

PCG64::~PCG64()
{
}

void PCG64::srand(uint64_t seed0, uint64_t seed1)
{
    state_.low_ = SplitMix::next(seed0);
    state_.high_ = SplitMix::next(seed1);
    while(0 == state_.low_ && 0 == state_.high_) {
        state_.low_ = SplitMix::next(state_.low_);
        state_.high_ = SplitMix::next(state_.high_);
    }
}

uint64_t PCG64::rand()
{
    next();
    return rotr64(state_.high_ ^ state_.low_, state_.high_ >> 58U);
}

double PCG64::frand()
{
    return to_real64(rand());
}

//--- SplitMix
//--------------------------------------------
uint64_t SplitMix::next(uint64_t& state)
{
    state += 0x9E3779B97f4A7C15ULL;
    uint64_t t = state;
    t = (t ^ (t >> 30)) * 0xBF58476D1CE4E5B9ULL;
    t = (t ^ (t >> 27)) * 0x94D049BB133111EBULL;
    return t ^ (t >> 31);
}

//--- RandWELL512
//---------------------------------------------
RandWELL512::RandWELL512()
    : index_(0)
{
    srand(CPPRNG_DEFAULT_SEED32);
}

RandWELL512::RandWELL512(uint32_t seed)
    : index_(0)
{
    srand(seed);
}

RandWELL512::~RandWELL512()
{
}

void RandWELL512::srand(uint32_t seed)
{
    state_[0] = scramble(seed);
    for(uint32_t i = 1; i < N; ++i) {
        state_[i] = scramble(seed);
    }
}

uint32_t RandWELL512::rand()
{
    uint32_t a, b, c, d;

    a = state_[index_];
    c = state_[(index_ + 13) & 15];
    b = a ^ c ^ (a << 16) ^ (c << 15);
    c = state_[(index_ + 9) & 15];
    c ^= c >> 11;
    a = state_[index_] = b ^ c;
    d = a ^ ((a << 5) & 0xDA442D24UL);
    index_ = (index_ + 1) & 15;
    a = state_[index_];
    state_[index_] = a ^ b ^ d ^ (a << 2) ^ (b << 18) ^ (c << 28);
    return state_[index_];
}

float RandWELL512::frand()
{
    return to_real32(rand());
}

//--- SFMT
//---------------------------------------------------------
namespace
{
    inline void rshift128(w128_t* CPPRND_RESTRICT dst, const w128_t* CPPRND_RESTRICT src, uint32_t shift)
    {
        uint64_t th, tl, oh, ol;
        th = ((uint64_t)src->u_[3] << 32) | ((uint64_t)src->u_[2]);
        tl = ((uint64_t)src->u_[1] << 32) | ((uint64_t)src->u_[0]);

        oh = th >> (shift * 8);
        ol = tl >> (shift * 8);
        ol |= th << (64 - shift * 8);
        dst->u_[1] = (uint32_t)(ol >> 32);
        dst->u_[0] = (uint32_t)ol;
        dst->u_[3] = (uint32_t)(oh >> 32);
        dst->u_[2] = (uint32_t)oh;
    }

    inline void lshift128(w128_t* CPPRND_RESTRICT dst, const w128_t* CPPRND_RESTRICT src, uint32_t shift)
    {
        uint64_t th, tl, oh, ol;
        th = ((uint64_t)src->u_[3] << 32) | ((uint64_t)src->u_[2]);
        tl = ((uint64_t)src->u_[1] << 32) | ((uint64_t)src->u_[0]);

        oh = th << (shift * 8);
        ol = tl << (shift * 8);
        oh |= tl >> (64 - shift * 8);
        dst->u_[1] = (uint32_t)(ol >> 32);
        dst->u_[0] = (uint32_t)ol;
        dst->u_[3] = (uint32_t)(oh >> 32);
        dst->u_[2] = (uint32_t)oh;
    }

    template<class T>
    inline void do_recursion(w128_t* CPPRND_RESTRICT r, const w128_t* CPPRND_RESTRICT a, const w128_t* CPPRND_RESTRICT b, const w128_t* CPPRND_RESTRICT c, const w128_t* CPPRND_RESTRICT d)
    {
        w128_t x;
        w128_t y;

        lshift128(&x, a, T::SFMT_SL2);
        rshift128(&y, c, T::SFMT_SR2);
        r->u_[0] = a->u_[0] ^ x.u_[0] ^ ((b->u_[0] >> T::SFMT_SR1) & T::SFMT_MSK1)
                   ^ y.u_[0] ^ (d->u_[0] << T::SFMT_SL1);
        r->u_[1] = a->u_[1] ^ x.u_[1] ^ ((b->u_[1] >> T::SFMT_SR1) & T::SFMT_MSK2)
                   ^ y.u_[1] ^ (d->u_[1] << T::SFMT_SL1);
        r->u_[2] = a->u_[2] ^ x.u_[2] ^ ((b->u_[2] >> T::SFMT_SR1) & T::SFMT_MSK3)
                   ^ y.u_[2] ^ (d->u_[2] << T::SFMT_SL1);
        r->u_[3] = a->u_[3] ^ x.u_[3] ^ ((b->u_[3] >> T::SFMT_SR1) & T::SFMT_MSK4)
                   ^ y.u_[3] ^ (d->u_[3] << T::SFMT_SL1);
    }

#ifdef SFMT_SSE
    template<class T>
    inline __m128i do_recursion_sse(
        const __m128i& a,
        const __m128i& b,
        const __m128i& c,
        const __m128i& d)
    {
        __m128i v, x, y, z;
        y = _mm_srli_epi32(b, T::SFMT_SR1);
        z = _mm_srli_si128(c, T::SFMT_SR2);
        v = _mm_slli_epi32(d, T::SFMT_SL1);
        z = _mm_xor_si128(z, a);
        z = _mm_xor_si128(z, v);
        x = _mm_slli_si128(a, T::SFMT_SL2);
        y = _mm_and_si128(y, sse2_param_mask.si);
        z = _mm_xor_si128(z, x);
        z = _mm_xor_si128(z, y);
        return z;
    }
#endif

    template<class T>
    void genrand_all(typename T::state_t* state)
    {
        w128_t* r1 = &state->state_[T::SFMT_N - 2];
        w128_t* r2 = &state->state_[T::SFMT_N - 1];
        uint32_t i;
        for(i = 0; i < T::SFMT_N - T::SFMT_POS1; ++i) {
            do_recursion<T>(&state->state_[i], &state->state_[i],
                            &state->state_[i + T::SFMT_POS1], r1, r2);
            r1 = r2;
            r2 = &state->state_[i];
        }
        for(; i < T::SFMT_N; ++i) {
            do_recursion<T>(&state->state_[i], &state->state_[i],
                            &state->state_[i + T::SFMT_POS1 - T::SFMT_N], r1, r2);
            r1 = r2;
            r2 = &state->state_[i];
        }
    }

#ifdef SFMT_SSE
    template<class T>
    void genrand_all_sse(typename T::state_t* state)
    {
        w128_t* pstate = state->state_;
        __m128i r1 = _mm_loadu_si128(&pstate[T::SFMT_N - 2]);
        __m128i r2 = _mm_loadu_si128(&pstate[T::SFMT_N - 1]);
        uint32_t i;
        for(i = 0; i < T::SFMT_N - T::SFMT_POS1; ++i) {
            __m128i x0 = _mm_loadu_si128(&pstate[i]);
            __m128i x1 = _mm_loadu_si128(&pstate[i + T::SFMT_POS1]);
            __m128i r = do_recursion_sse<T>(x0, x1, r1, r2);
            r1 = r2;
            r2 = r;
        }
        for(; i < T::SFMT_N; ++i) {
            __m128i x0 = _mm_loadu_si128(&pstate[i]);
            __m128i x1 = _mm_loadu_si128(&pstate[i + T::SFMT_POS1 - T::SFMT_N]);
            __m128i r = do_recursion_sse<T>(x0, x1, r1, r2);
            r1 = r2;
            r2 = r;
        }
    }
#endif

    inline uint32_t indexof(uint32_t x)
    {
        return x;
    }

    template<class T>
    void period_certification(typename T::state_t* state)
    {
        uint32_t inner = 0;
        uint32_t i, j;
        uint32_t work;
        uint32_t* psfmt32 = &state->state_[0].u_[0];
        const uint32_t parity[4] = {T::SFMT_PARITY1, T::SFMT_PARITY2,
                                    T::SFMT_PARITY3, T::SFMT_PARITY4};

        for(i = 0; i < 4; ++i) {
            inner ^= psfmt32[indexof(i)] & parity[i];
        }
        for(i = 16; i > 0; i >>= 1) {
            inner ^= inner >> i;
        }
        inner &= 1;
        /* check OK */
        if(inner == 1) {
            return;
        }
        /* check NG, and modification */
        for(i = 0; i < 4; ++i) {
            work = 1;
            for(j = 0; j < 32; ++j) {
                if((work & parity[i]) != 0) {
                    psfmt32[indexof(i)] ^= work;
                    return;
                }
                work = work << 1;
            }
        }
    }

    template<class T>
    void init_rand(typename T::state_t* state, uint32_t seed)
    {
        uint32_t* psfmt32 = &state->state_[0].u_[0];
        psfmt32[indexof(0)] = seed;
        for(uint32_t i = 1; i < T::SFMT_N32; ++i) {
            psfmt32[indexof(i)] = 1812433253UL * (psfmt32[indexof(i - 1)] ^ (psfmt32[indexof(i - 1)] >> 30))
                                  + i;
        }
        state->index_ = T::SFMT_N32;
        period_certification<T>(state);
    }

    template<class T>
    inline uint32_t genrand_uint32(typename T::state_t* state)
    {
        if(T::SFMT_N32 <= state->index_) {
            state->index_ = 0;
#ifdef SFMT_SSE
            genrand_all_sse<T>(state);
#else
            genrand_all<T>(state);
#endif
        }
        uint32_t* psfmt32 = &state->state_[0].u_[0];
        uint32_t r = psfmt32[state->index_];
        ++state->index_;
        return r;
    }

} // namespace

SFMT19937::SFMT19937()
{
    init_rand<sfmt_param19937>(&state_, CPPRNG_DEFAULT_SEED32);
}

SFMT19937::SFMT19937(uint32_t seed)
{
    init_rand<sfmt_param19937>(&state_, seed);
}

SFMT19937::~SFMT19937()
{
}

void SFMT19937::srand(uint32_t seed)
{
    init_rand<sfmt_param19937>(&state_, seed);
}

uint32_t SFMT19937::rand()
{
    return genrand_uint32<sfmt_param19937>(&state_);
}

float SFMT19937::frand()
{
    return to_real32(rand());
}

//--- MELG
//---------------------------------------------------------
namespace
{
    inline uint64_t MAT3NEG(uint64_t t, uint64_t v)
    {
        return (v ^ (v << ((t))));
    }

    inline uint64_t MAT3POS(uint64_t t, uint64_t v)
    {
        return (v ^ (v >> ((t))));
    }

    template<class T>
    uint64_t melg_case_1(typename T::state_t* state);
    template<class T>
    uint64_t melg_case_2(typename T::state_t* state);
    template<class T>
    uint64_t melg_case_3(typename T::state_t* state);
    template<class T>
    uint64_t melg_case_4(typename T::state_t* state);

    template<class T>
    uint64_t melg_case_1(typename T::state_t* state)
    {
        uint64_t* melg = state->state_;
        uint64_t melgi = state->index_;
        uint64_t x = (melg[melgi] & T::MASKU) | (melg[melgi + 1] & T::MASKL);
        state->lung_ = (x >> 1) ^ T::mag01[(int32_t)(x & 1ULL)] ^ melg[melgi + T::MM] ^ MAT3NEG(23, state->lung_);
        melg[melgi] = x ^ MAT3POS(33, state->lung_);
        x = melg[melgi] ^ (melg[melgi] << T::SHIFT1);
        x = x ^ (melg[melgi + T::LAG1] & T::MASK1);
        ++state->index_;
        if((T::NN - T::MM) == state->index_) {
            state->function_ = melg_case_2<T>;
        }
        return x;
    }

    template<class T>
    uint64_t melg_case_2(typename T::state_t* state)
    {
        uint64_t* melg = state->state_;
        uint64_t melgi = state->index_;
        uint64_t x = (melg[melgi] & T::MASKU) | (melg[melgi + 1] & T::MASKL);
        state->lung_ = (x >> 1) ^ T::mag01[(int32_t)(x & 1ULL)] ^ melg[melgi + (T::MM - T::NN)] ^ MAT3NEG(23, state->lung_);
        melg[melgi] = x ^ MAT3POS(33, state->lung_);
        x = melg[melgi] ^ (melg[melgi] << T::SHIFT1);
        x = x ^ (melg[melgi + T::LAG1] & T::MASK1);
        ++state->index_;
        if(T::LAG1over == state->index_) {
            state->function_ = melg_case_3<T>;
        }
        return x;
    }

    template<class T>
    uint64_t melg_case_3(typename T::state_t* state)
    {
        uint64_t* melg = state->state_;
        uint64_t melgi = state->index_;
        uint64_t x = (melg[melgi] & T::MASKU) | (melg[melgi + 1] & T::MASKL);
        state->lung_ = (x >> 1) ^ T::mag01[(x & 1ULL)] ^ melg[melgi + (T::MM - T::NN)] ^ MAT3NEG(23, state->lung_);
        melg[melgi] = x ^ MAT3POS(33, state->lung_);
        x = melg[melgi] ^ (melg[melgi] << T::SHIFT1);
        x = x ^ (melg[melgi - T::LAG1over] & T::MASK1);
        ++state->index_;
        if((T::NN - 1) == state->index_) {
            state->function_ = melg_case_4<T>;
        }
        return x;
    }

    template<class T>
    uint64_t melg_case_4(typename T::state_t* state)
    {
        uint64_t* melg = state->state_;
        uint64_t melgi = state->index_;
        uint64_t x = (melg[T::NN - 1] & T::MASKU) | (melg[0] & T::MASKL);
        state->lung_ = (x >> 1) ^ T::mag01[(x & 1ULL)] ^ melg[T::MM - 1] ^ MAT3NEG(23, state->lung_);
        melg[T::NN - 1] = x ^ MAT3POS(33, state->lung_);
        x = melg[melgi] ^ (melg[melgi] << T::SHIFT1);
        x = x ^ (melg[melgi - T::LAG1over] & T::MASK1);
        state->index_ = 0;
        state->function_ = melg_case_1<T>;
        return x;
    }

    template<class T>
    void melg_init_rand(typename T::state_t* state, uint64_t seed)
    {
        uint64_t* melg = state->state_;
        uint64_t melgi;
        melg[0] = seed;
        for(melgi = 1; melgi < T::NN; melgi++) {
            melg[melgi] = (6364136223846793005ULL * (melg[melgi - 1] ^ (melg[melgi - 1] >> 62)) + melgi);
        }
        state->lung_ = (6364136223846793005ULL * (melg[melgi - 1] ^ (melg[melgi - 1] >> 62)) + melgi);
        state->index_ = 0;
        state->function_ = melg_case_1<T>;
    }
} // namespace

MELG19937::MELG19937()
{
    melg_init_rand<melg_param19937>(&state_, CPPRNG_DEFAULT_SEED64);
}

MELG19937::MELG19937(uint64_t seed)
{
    melg_init_rand<melg_param19937>(&state_, seed);
}
MELG19937::MELG19937(uint64_t size, const uint64_t* seeds)
{
    srand(size, seeds);
}

MELG19937::~MELG19937()
{
}

void MELG19937::srand(uint64_t seed)
{
    melg_init_rand<melg_param19937>(&state_, seed);
}

void MELG19937::srand(uint64_t size, const uint64_t* seeds)
{
    assert(nullptr != seeds);
    melg_init_rand<melg_param19937>(&state_, 19650218ULL);
    uint64_t i = 1;
    uint64_t j = 0;
    uint64_t k = size < melg_param19937::NN ? melg_param19937::NN : size;
    uint64_t* melg = state_.state_;

    for(; 0 < k; k--) {
        melg[i] = (melg[i] ^ ((melg[i - 1] ^ (melg[i - 1] >> 62)) * 3935559000370003845ULL))
                  + seeds[j] + j; /* non linear */
        i++;
        j++;
        if(melg_param19937::NN <= i) {
            melg[0] = melg[melg_param19937::NN - 1];
            i = 1;
        }
        if(size <= j) {
            j = 0;
        }
    }
    for(k = melg_param19937::NN - 1; 0 < k; k--) {
        melg[i] = (melg[i] ^ ((melg[i - 1] ^ (melg[i - 1] >> 62)) * 2862933555777941757ULL))
                  - i; /* non linear */
        i++;
        if(melg_param19937::NN <= i) {
            melg[0] = melg[melg_param19937::NN - 1];
            i = 1;
        }
    }
    state_.lung_ = (state_.lung_ ^ ((melg[melg_param19937::NN - 1] ^ (melg[melg_param19937::NN - 1] >> 62)) * 2862933555777941757ULL))
                   - melg_param19937::NN; /* non linear */
    melg[0] = (melg[0] | (1ULL << 63));   /* MSB is 1; assuring non-zero initial array. Corrected.  */
    state_.index_ = 0;
}

uint64_t MELG19937::rand()
{
    return state_.function_(&state_);
}

double MELG19937::frand()
{
    return to_real64(rand());
}

#ifdef CPPRND_CRPNG
//--- CPRNG
//---------------------------------------------------------
bool crypt_rand(uint32_t size, void* buffer)
{
    assert(nullptr != buffer);
#if _WIN32
    HMODULE handle = LoadLibraryA("Advapi32.dll");
    if(nullptr == handle) {
        return false;
    }
    FARPROC procAddress = GetProcAddress(handle, "SystemFunction036");
    BOOLEAN result = (*(BOOLEAN(*)(PVOID, ULONG))procAddress)(buffer, size);
    FreeLibrary(handle);
    return TRUE == result;
#else
    int32_t fd = open("/dev/random", O_RDONLY);
    if(fd < 0) {
        fd = open("/dev/urandom", O_RDONLY);
        if(fd < 0) {
            return false;
        }
    }
    read(fd, buffer, size);
    close(fd);
    return true;
#endif
}
#endif

//--- RandomBinarySelect
//---------------------------------------------------------
RandomBinarySelect::RandomBinarySelect()
    : capacity_(0)
    , size_(0)
    , weights_(nullptr)
{
}

RandomBinarySelect::~RandomBinarySelect()
{
    delete[] weights_;
    weights_ = nullptr;
}

uint32_t RandomBinarySelect::size() const
{
    return size_;
}

void RandomBinarySelect::build(uint32_t size, const float* weights)
{
    if(capacity_ < size) {
        while(capacity_ < size) {
            capacity_ += 16UL;
        }
        delete[] weights_;
        weights_ = new float[capacity_];
    }
    size_ = size;
    // Kahan's summation
    float sum = 0.0f;
    float c = 0.0f;
    for(uint32_t i = 0; i < size_; ++i) {
        float x = weights[i] - c;
        float t = sum + x;
        c = (t - sum) - x;
        sum = t;
        weights_[i] = sum;
    }
    if(sum < 1.0e-7f) {
        return;
    }
    float scale = 1.0f / sum;
    for(uint32_t i = 0; i < size_; ++i) {
        weights_[i] *= scale;
    }
}

RandomAliasSelect::RandomAliasSelect()
    : capacity_(0)
    , size_(0)
    , weights_(nullptr)
    , aliases_(nullptr)
{
}

RandomAliasSelect::~RandomAliasSelect()
{
    aliases_ = nullptr;
    delete[] weights_;
}

uint32_t RandomAliasSelect::size() const
{
    return size_;
}

void RandomAliasSelect::build(uint32_t size, float* weights)
{
    assert(0 <= size);
    if(capacity_ < size) {
        do {
            capacity_ += 16UL;
        } while(capacity_ < size);
        delete[] weights_;
        weights_ = new float[capacity_ * 3];
        aliases_ = reinterpret_cast<uint32_t*>(weights_ + capacity_);
    }
    size_ = size;

    // Kahan's summation
    float total = 0.0f;
    {
        float c = 0.0f;
        for(uint32_t i = 0; i < size; ++i) {
            float x = weights[i] - c;
            float t = total + x;
            c = (t - total) - x;
            total = t;
        }
    }

    float average = total / size_;
    float scale = (1.0e-7f < total) ? static_cast<float>(size_) / total : 0.0f;
    uint32_t* indices = reinterpret_cast<uint32_t*>(aliases_ + capacity_);

    int32_t underfull = -1;
    int32_t overfull = static_cast<int32_t>(size_);
    for(uint32_t i = 0; i < size_; ++i) {
        if(average <= weights[i]) {
            --overfull;
            indices[overfull] = i;
        } else {
            ++underfull;
            indices[underfull] = i;
        }
    }
    while(0 <= underfull && overfull < static_cast<int32_t>(size_)) {
        uint32_t under = indices[underfull];
        --underfull;
        uint32_t over = indices[overfull];
        ++overfull;
        aliases_[under] = over;
        weights_[under] = weights[under] * scale;
        weights[over] += weights[under] - average;
        if(weights[over] < average) {
            ++underfull;
            indices[underfull] = over;
        } else {
            --overfull;
            indices[overfull] = over;
        }
    }
    while(0 <= underfull) {
        weights_[indices[underfull]] = 1.0f;
        --underfull;
    }
    while(overfull < static_cast<int32_t>(size_)) {
        weights_[indices[overfull]] = 1.0f;
        ++overfull;
    }
}
} // namespace cpprnd
