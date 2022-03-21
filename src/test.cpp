#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "Peptides.h"

SCENARIO("normal test") {
    char *ptest[8];
    ptest[0] = strdup("foo.exe");
    ptest[1] = strdup("min.fasta");
    ptest[2] = strdup("-o");
    ptest[3] = strdup("min_mimic.fasta");
    ptest[4] = strdup("-m");
    ptest[5] = strdup("10");
    ptest[6] = strdup("-I");
    ptest[7] = nullptr;
    unsigned int minLength = 0u;
    unsigned int maxTries = 50u;
    Peptides peptides(minLength, {}, maxTries);
    peptides.parseOptions(7, ptest);
    peptides.run();

    std::vector<string> expected{
            ">mimic|Random_1|shuffle_1",
            "AICK",
            ">mimic|Random_1|shuffle_2",
            "ARIK",
            ">mimic|Random_1|shuffle_3",
            "ACPK",
            ">mimic|Random_1|shuffle_4",
            "APLK",
            ">mimic|Random_1|shuffle_5",
            "ACRK",
            ">mimic|Random_1|shuffle_6",
            "ACQK",
            ">mimic|Random_1|shuffle_7",
            "AARK",
            ">mimic|Random_1|shuffle_8",
            "ACYK",
            ">mimic|Random_1|shuffle_9",
            "AYIK",
            ">mimic|Random_1|shuffle_10",
            "AMIK",

    };

    string line;
    int i = 0;
    ifstream myfile("min_mimic.fasta");
    while (getline(myfile, line)) {
        INFO("shuffle: " + std::to_string(i/2+1))
        CHECK(line == expected[i]);
        ++i;
    }
    myfile.close();
}
