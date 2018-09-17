/*! \file scale.cpp
 *  \brief Enter description here.
 *  \author Georgi Gerganov
 */

#include <fstream>
#include <cstdlib>
#include <vector>

bool g_terminate = false;

int main(int argc, char ** argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s input.kbd output.kbd nf0 nf1\n", argv[0]);
        return -127;
    }

    int nf0 = atoi(argv[3]);
    int nf1 = atoi(argv[4]);

    if (nf0%2 == 0) {
        printf("nf0 must be odd\n");
        return -1;
    }

    if (nf1%2 == 0) {
        printf("nf1 must be odd\n");
        return -1;
    }

    std::ifstream fin(argv[1], std::ios::binary);
    std::ofstream fout(argv[2], std::ios::binary);

    int bufferSize_frames = 1;
    fin.read((char *)(&bufferSize_frames), sizeof(bufferSize_frames));

    if (bufferSize_frames != nf0) {
        printf("Invalid number nf0 = %d.  Expected %d for this file\n", nf0, bufferSize_frames);
        return -1;
    }

    fout.write((char *)(&nf1), sizeof(nf1));

    int nadd = (nf1 - nf0)/2;
    printf("nadd = %d\n", nadd);

    while (true) {
        int key = -1;
        fin.read((char *)(&key), sizeof(key));
        if (fin.eof()) break;
        fout.write((char *)(&key), sizeof(key));

        std::vector<char> f0(4*1024);
        std::fill(f0.begin(), f0.end(), 0);
        std::vector<char> buf(nf0*4*1024);
        fin.read((char *)(buf.data()), buf.size());
        for (int i = 0; i < nadd; ++i) fout.write((char *)(f0.data()), f0.size());
        fout.write((char *)(buf.data()), buf.size());
        for (int i = 0; i < nadd; ++i) fout.write((char *)(f0.data()), f0.size());
    }

    fin.close();
    fout.close();

    return 0;
}
