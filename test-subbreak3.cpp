#include "subbreak3.h"

int main(int argc, char ** argv) {
    printf("Usage: %s n-gram.txt\n", argv[0]);
    if (argc < 2) {
        return -1;
    }

    Cipher::TFreqMap freqMap;

    if (Cipher::loadFreqMapBinary((std::string(argv[1]) + ".binary").c_str(), freqMap) == false) {
        if (Cipher::loadFreqMap(argv[1], freqMap) == false) {
            return -1;
        }

        Cipher::saveFreqMapBinary((std::string(argv[1]) + ".binary").c_str(), freqMap);
    }

    //if (Cipher::loadFreqMap(argv[1], freqMap) == false) {
    //    return -1;
    //}

    std::string plain = R"(
Dave found joy in the daily routine of life. He awoke at the same time, ate the same breakfast and drove the same commute. He worked at a job that never seemed to change and he got home at 6 pm sharp every night. It was who he had been for the last ten years and he had no idea that was all about to change.
    )";

    Cipher::TParameters params;
    Cipher::TResult result;

    Cipher::encryptExact(params, plain, result.clusters);
    Cipher::printEncoded(result.clusters);

    {
        std::map<TClusterId, char> answer;
        for (int i = 0; i < (int) result.clusters.size(); ++i) {
            answer[result.clusters[i]] = plain[i];
        }

        for (auto & p : answer) {
            printf("  answer: '%c' -> '%c'\n", Cipher::getEncodedChar(p.first), p.second);
        }
    }

    {
        auto tStart = std::chrono::high_resolution_clock::now();

        Cipher::beamSearch(params, freqMap, result);

        auto tEnd = std::chrono::high_resolution_clock::now();
        auto tDiff = std::chrono::duration_cast<std::chrono::milliseconds>(tEnd - tStart);
        printf("Time: %ld ms\n", tDiff.count());
    }

    return 0;
}
