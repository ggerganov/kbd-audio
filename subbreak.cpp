/*! \file subbreak.cpp
 *  \brief Substitution cipher breaker
 *  \author Georgi Gerganov
 */

#include "constants.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <tuple>

using TCode = int32_t;
using TProb = double;
using TGramLen = int;
using TFreqMap = std::tuple<TGramLen, std::vector<TProb>>;
using TAlphabet = std::vector<char>;

TCode calcCode(const char * data, int n) {
    TCode res = 0;
    do { res <<= 5; res += kCharToInt[*data++]; } while (--n > 0);
    return res;
}

bool loadFreqMap(const char * fname, TFreqMap & res) {
    auto & [len, fmap] = res;
    len = 0;
    fmap.clear();

    printf("[+] Loading n-gram file '%s'\n", fname);
    std::ifstream fin(fname);

    if (fin.good() == false) {
        printf("    Failed to open file\n");
        return false;
    }

    std::string gram;
    int32_t nfreq = 0;
    int64_t nTotal = 0;

    while (true) {
        fin >> gram >> nfreq;
        if (fin.eof()) break;

        nTotal += nfreq;
    }

    fin.clear();
    fin.seekg(0);

    while (true) {
        fin >> gram >> nfreq;
        if (fin.eof()) break;

        if (len == 0) {
            len = gram.size();
            fmap.resize(1 << (5*len), 0.0);
        } else if (len != gram.size()) {
            printf("Error: loaded n-grams with vaying lengths\n");
            return false;
        }

        TCode idx = ::calcCode(gram.data(), len);
        if (fmap[idx] != 0.0) {
            printf("Error: duplicate n-gram '%s'\n", gram.c_str());
            return false;
        }
        fmap[idx] = std::log10(((double)(nfreq))/nTotal);
    }
    printf("    Total n-grams loaded = %g\n", (double) nTotal);

    double pMin = std::log(0.5/nTotal);
    printf("    P-min = %g\n", pMin);
    for (auto & p : fmap) {
        if (p >= 0.0) p = pMin;
    }

    return true;
}

TAlphabet getAlphabetRandom(int seed = 0) {
    TAlphabet res;

    for (int i = 0; i <= 26; ++i) res.push_back(i);

    std::random_device rd;
    std::mt19937 g(seed > 0 ? seed : rd());
    std::shuffle(res.begin() + 1, res.end(), g);

    return res;
}

bool encipher(const std::string & plain, TAlphabet & alphabet, std::string & enc) {
    alphabet = getAlphabetRandom(1234);

    enc = plain;
    int n = plain.length();
    for (int i = 0; i < n; ++i) {
        auto c = ::kCharToInt[plain[i]];
        enc[i] = enc[i] - c + alphabet[c];
    }

    return true;
}

bool translate(const TAlphabet & alphabet, const std::string & src, std::string & dst) {
    if (dst.size() != src.size()) return false;

    for (int i = 0; i < (int) src.size(); ++i) {
        if (src[i] == 0) continue;
        auto c = ::kCharToInt[src[i]];
        dst[i] = c > 0 ? src[i] - c + alphabet[c] : '.';
    }

    return true;
}

TProb calcScore0(const TFreqMap & freqMap, const std::string & txt) {
    TProb res = 0.0;

    auto len = txt.size();
    const auto & [n, fmap] = freqMap;

    int i1 = 0;

    int k = n;
    TCode curc = 0;
    while (k > 0) {
        if (i1 >= len) return 1e-100;
        auto c = ::kCharToInt[txt[i1++]];
        if (c > 0) {
            curc <<= 5;
            curc += c;
            --k;
        }
    }

    res += fmap[curc];
    while (true) {
        curc = curc & ((1 << 15) - 1);

        while (true) {
            if (i1 >= len) return res;
            auto c = ::kCharToInt[txt[i1++]];
            if (c > 0) {
                curc <<= 5;
                curc += c;
                break;
            }
        }

        res += fmap[curc];
    }

    return res;
}

TProb calcScore1(const TFreqMap & freqMap, const std::string & txt) {
    TProb res = 0.0;

    auto len = txt.size();
    const auto & [n, fmap] = freqMap;

    int i1 = 0;

    int k = n;
    TCode curc = 0;
    while (k > 0) {
        if (i1 >= len) return 1e-100;
        auto c = ::kCharToInt[txt[i1++]];
        curc <<= 5;
        curc += c;
        --k;
    }

    res += fmap[curc];
    while (true) {
        curc = curc & ((1 << 15) - 1);

        if (i1 >= len) return res;
        auto c = ::kCharToInt[txt[i1++]];
        curc <<= 5;
        curc += c;

        res += fmap[curc];
    }

    return res;
}

auto calcScore = calcScore0;

bool decrypt(const TFreqMap & freqMap, const std::string & enc, std::string & res) {
    TAlphabet besta = getAlphabetRandom();

    auto cure = enc;
    translate(besta, enc, cure);
    auto bestp = calcScore(freqMap, cure);
    while (true) {
        auto itera = getAlphabetRandom();
        translate(itera, enc, cure);
        auto iterp = calcScore(freqMap, cure);
        for (int i = 0; i < 1000; ++i) {
            auto cura = itera;
            int a0 = rand()%26 + 1;
            int a1 = rand()%26 + 1;
            while (a0 == a1) {
                a0 = rand()%26 + 1;
                a1 = rand()%26 + 1;
            }

            std::swap(cura[a0], cura[a1]);

            translate(cura, enc, cure);
            auto curp = calcScore(freqMap, cure);
            if (curp > iterp) {
                iterp = curp;
                itera = cura;
                i = 0;
            }
        }
        if (iterp > bestp) {
            besta = itera;
            bestp = iterp;
            translate(besta, enc, cure);

            printf("[+] Best score = %g\n", bestp);
            printf("    Alphabet:  '");
            for (int i = 'a'; i <= 'z'; ++i) {
                auto c = ::kCharToInt[i];
                printf("%c", i - c + besta[c]);
            }
            printf("'\n");
            printf("    Decrypted: '%s'\n", cure.c_str());
        }
    }

    return true;
}

int main(int argc, char ** argv) {
    printf("Usage: %s n-gram.txt\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    TFreqMap freqMap;
    if (loadFreqMap(argv[1], freqMap) == false) {
        return -1;
    }

    printf("[+] Score of 'XZZS'        = %g\n", calcScore(freqMap, "XZZS"));
    printf("[+] Score of 'HELLO WORLD' = %g\n", calcScore(freqMap, "HELLO WORLD"));

    TAlphabet alphabet;
    //std::string plain = "Hello world! By spite about do of do allow blush. Additions in conveying or collected objection in. Suffer few desire wonder her object hardly nearer.";
    std::string plain = R"(
    Is he staying arrival address earnest. To preference considered it themselves inquietude collecting estimating.
    View park for why gay knew face. Next than near to four so hand. Times so do he downs me would.
    Witty abode party her found quiet law. They door four bed fail now have.)";

    std::string enc;
    encipher(plain, alphabet, enc);

    printf("[+] Alphabet: 'abcdefghijklmnopqrstuvwxyz'\n");
    printf("[+] Alphabet: '");
    for (int i = 'a'; i <= 'z'; ++i) {
        auto c = ::kCharToInt[i];
        printf("%c", i - c + alphabet[c]);
    }
    printf("'\n");

    printf("[+] Plain:    '%s'\n", plain.c_str());
    printf("[+] Enc:      '%s'\n", enc.c_str());

    std::string decrypted;
    decrypt(freqMap, enc, decrypted);

    return 0;
}
