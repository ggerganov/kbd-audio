/*! \file common.h
 *  \brief Common types and functions
 *  \author Georgi Gerganov
 */

#pragma once

#include <map>
#include <string>
#include <cstring>

static std::map<std::string, std::string> parseCmdArguments(int argc, char ** argv) {
    int last = argc;
    std::map<std::string, std::string> res;
    for (int i = 1; i < last; ++i) {
        if (argv[i][0] == '-') {
            if (strlen(argv[i]) > 1) {
                res[std::string(1, argv[i][1])] = strlen(argv[i]) > 2 ? argv[i] + 2 : "";
            }
        }
    }

    return res;
}

