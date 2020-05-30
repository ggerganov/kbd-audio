/*! \file non-exact-subbreak2.cpp
 *  \brief New attempt at breaking non-exact substitution cipher
 *  \author Georgi Gerganov
 */

#include "subbreak2.h"

int main(int argc, char ** argv) {
    printf("Usage: %s n-gram.txt\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    srand(time(0));

    Cipher::TParameters params;

    Cipher::TFreqMap freqMap;
    if (Cipher::loadFreqMap(argv[1], freqMap) == false) {
        return -1;
    }

    std::string plain;

    plain = R"(
As far as services go, only two steady contributors of revenue
streams keep swelling without Apple having to charge
subscription fees. One is the money paid by Google parent Alphabet
for searches made through Apple products such as the Safari browser and Siri.
    )";

    TClusters enc;
    Cipher::encryptExact(params, plain, enc);

    params.maxClusters = 27;
    for (auto & c : enc) {
        params.maxClusters = std::max(params.maxClusters, c + 1);
    }

    TClusterToLetterMap clMap;
    Cipher::subbreak(params, freqMap, enc, clMap);

    Cipher::printText(enc, clMap);

    return 0;
}
