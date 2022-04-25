#include "subbreak3.h"

int main(int argc, char ** argv) {
    printf("Usage: %s n-gram.dat n-gram-compressed.dat\n", argv[0]);
    if (argc < 3) {
        return -1;
    }

    Cipher::TFreqMap freqMap;

    printf("[+] Reading n-grams from '%s'\n", argv[1]);
    if (Cipher::loadFreqMap(argv[1], freqMap) == false) {
        return -1;
    }

    printf("[+] Writing compressed n-grams to '%s'\n", argv[2]);
    Cipher::saveFreqMapBinary(argv[2], freqMap);

    return 0;
}
