/*! \file subbreak.cpp
 *  \brief Substitution cipher breaker
 *  \author Georgi Gerganov
 */

#include "subbreak.h"

int main(int argc, char ** argv) {
    printf("Usage: %s n-gram.txt\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    srand(time(0));

    TFreqMap freqMap;
    if (loadFreqMap(argv[1], freqMap) == false) {
        return -1;
    }

    printf("[+] Score of 'XZZS'        = %g\n", calcScore(freqMap, "XZZS"));
    printf("[+] Score of 'HELLO WORLD' = %g\n", calcScore(freqMap, "HELLO WORLD"));

    TAlphabet alphabet;
    std::string plain = R"(
As far as services go, only two steady contributors of revenue
streams keep swelling without Apple having to charge
subscription fees. One is the money paid by Google parent Alphabet
for searches made through Apple products such as the Safari browser and Siri.
    )";

    std::string enc;
    encrypt(plain, alphabet, enc);

    printf("[+] Alphabet: 'abcdefghijklmnopqrstuvwxyz'\n");
    printf("[+] Alphabet: '");
    for (int i = 'a'; i <= 'z'; ++i) {
        auto c = ::myCharToInt[i];
        printf("%c", i - c + alphabet[c]);
    }
    printf("'\n");

    printf("[+] Plain:    '%s'\n", plain.c_str());
    {
        TAlphabet sola = alphabet;
        for (int i = 0; i <= kN; ++i) {
            sola[alphabet[i]] = i;
        }
        std::string solt = enc;
        translate(sola, enc, solt);
        printf("[+] Sol: ");
        printText(solt);
        printf("[+] Score = %g\n", calcScore(freqMap, solt));
    }
    printf("[+] Enc:      ");
    printText(enc);

    std::string decrypted;
    decrypt(freqMap, enc, decrypted);

    return 0;
}
