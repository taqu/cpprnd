# cpprnd
Useful PRNGs with independent contexts, and cleaned codes for modern C++. 

| Name      | Bits   | Period  | Streams |
|:----------|:-------|:--------|:--------|
| PCG32     | 32bits | 2^64    | 1       |
| PCGS32    | 32bits | 2^64    | 2^31    |
| PCG64     | 64bits | 2^128   | 1       |
| SFMT19937 | 32bits | 2^19937 | 1       |
| MELG19937 | 64bits | 2^19937 | 1       |

# References

1. [SFMT: http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/SFMT/](http://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/SFMT/)
2. [MELG: https://github.com/sharase/melg-64](https://github.com/sharase/melg-64)
3. [PCG: https://www.pcg-random.org/index.html](https://www.pcg-random.org/index.html)
4. [WELL: http://www.iro.umontreal.ca/~panneton/WELLRNG.html](http://www.iro.umontreal.ca/~panneton/WELLRNG.html)
5. [SplitMix: https://dl.acm.org/doi/10.1145/2660193.2660195](https://dl.acm.org/doi/10.1145/2660193.2660195)

# License
A part of SFMT and MELG is destributed under other licences, see under the `doc`.
This software is distributed under two licenses 'The MIT License' or 'Public Domain', choose whichever you like.
