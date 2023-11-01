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
    Peptides peptides(minLength, {}, std::mt19937{}, maxTries);
    peptides.parseOptions(7, ptest);
    peptides.run();

    std::vector<string> expected{
            ">mimic|Random_1|shuffle_1",
            "ACIKDCE",
            ">mimic|Random_1|shuffle_2",
            "AIFKDQE",
            ">mimic|Random_1|shuffle_3",
            "AIQKDAE",
            ">mimic|Random_1|shuffle_4",
            "AIPKDVE",
            ">mimic|Random_1|shuffle_5",
            "AIAKDPE",
            ">mimic|Random_1|shuffle_6",
            "AIVKDFE",
            ">mimic|Random_1|shuffle_7",
            "AIDKDME",
            ">mimic|Random_1|shuffle_8",
            "AIHKDGE",
            ">mimic|Random_1|shuffle_9",
            "AIEKDWE",
            ">mimic|Random_1|shuffle_10",
            "AIYKDNE",
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


map<char, double> testDist(){
    map<char, double> dist;
    dist['A'] = 0.081;
    dist['I'] = 0.057;
    dist['R'] = 0.057;
    dist['W'] = 0.013;
    return dist;
}
void checkSum(double d){
    // allow one bit difference
    CHECK(d >= std::nextafter(1.0, 0.0));
    CHECK(d <= std::nextafter(1.0, 2.0));
}
SCENARIO("Test AA dist"){
    WHEN("default Dist and no IL removal is used") {
        AminoAcidDist defaultDist(false);
        auto distSum = [](const map<char, double> &dist) {
            double sum = std::accumulate(dist.cbegin(), dist.cend(), 0.0,
                                         [](const double &cum, const std::pair<char, double> &pair) {
                                             return cum + pair.second;
                                         });
            return sum;
        };
        CHECK(defaultDist.getDist().size() == 18);
        checkSum(distSum(defaultDist.getDist()));
    }
    WHEN("default Dist and IL removal is used") {
        AminoAcidDist defaultDist(true);
        auto distSum = [](const map<char, double> &dist) {
            double sum = std::accumulate(dist.cbegin(), dist.cend(), 0.0,
                                         [](const double &cum, const std::pair<char, double> &pair) {
                                             return cum + pair.second;
                                         });
            return sum;
        };
        CHECK(defaultDist.getDist().size() == 16);
        checkSum(distSum(defaultDist.getDist()));
    }
    WHEN("default new dist is set") {
        AminoAcidDist defaultDist(false);
        defaultDist.setDist(testDist(), false);
        auto distSum = [](const map<char, double> &dist) {
            double sum = std::accumulate(dist.cbegin(), dist.cend(), 0.0,
                                         [](const double &cum, const std::pair<char, double> &pair) {
                                             return cum + pair.second;
                                         });
            return sum;
        };
        CHECK(defaultDist.getDist().size() == 3);
        checkSum(distSum(defaultDist.getDist()));
    }
}

SCENARIO("Test Abs AA dist") {
    AbsAminoAcidDist dist;
    dist.add('A');
    dist.add('A');
    dist.add('K');
    dist.add('L');
    auto relDist = dist.getDist();
    CHECK(relDist.at('A') == 0.5);
    CHECK(relDist.at('K') == 0.25);
    CHECK(relDist.at('L') == 0.25);
}