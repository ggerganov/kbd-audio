/*! \file non-exact-subbreak.cpp
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

    plain = R"(as far as services go only two stead8y contributo23rs of revenue
streams ke1ep swelling witho7ut Apple having to charge
subscription fees one is the4 money paid by google parent Alph-abet
for sear3ches made through apple prod6u4cts such as the Safa=ri browser a0nd Siri)";

    plain = R"(and shipping at Formosa and Luzon, the Third Fleet entered the South China Sea and aircraft from its carriers attacked Japanese shipping off)";

    plain = R"(resolutions and the protection of people's rights" and was rejected by 85% of voters.[1] The third was on a popular initiative on the)";

    plain = R"(The length of the forewings is about 8.7 mm. The ground colour of the forewings is whitish, strigulated (finely streaked) and suffused with pale orange)";

    //plain = R"(
    //As far as services go, only two steady contributors of revenue
    //streams keep swelling without Apple having to charge
    //subscription fees. One is the money paid by Google parent Alphabet
    //for searches made through Apple products such as the Safari browser and Siri
    //)";

    //std::string x;
    //for (int i = 0; i < (int) plain.size(); ++i) {
    //    while (rand()%1000 < 250) {
    //        x += rand()%10 + '1';
    //    }
    //    x += plain[i];
    //}

    //printf("x = '%s'\n", x.c_str());
    //plain = x;

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
