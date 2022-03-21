#define CATCH_CONFIG_MAIN

#include "catch2/catch.hpp"
#include "Peptides.h"

SCENARIO("test") {
    char *ptest[5];
    ptest[0] = strdup("foo.exe");
    ptest[1] = strdup("min.fasta");
    ptest[2] = strdup("-o");
    ptest[3] = strdup("min_mimic.fasta");
    ptest[4] = nullptr;
    Peptides peptides;
    peptides.parseOptions(4, ptest);
    peptides.run();

    std::vector<string> expected{">mimic|Random_1",
                                 "MTPESSPPEDELQWFLVSDSEPQKLLDDAPQESMPIFLWPATNLGPDPMDELSDSVDPPE",
                                 "QLENRMPAVAAAPPASVPAPASPAPAPTSAPQPSLWSEPPKTYYQSGGFRLGGLFTHSAK",
                                 "SNVTTLSPCYAKMALQFCKTQVWDTLPPCVPGSPTRVRAAMYIKQTVHSMQVERRCHPHE",
                                 "RCAILSDDPLQGPSHRVGNLERVDLYEDRNTFRHMGMVTYVNNNSMPGSCVHYTECCDES",
                                 "GYSPVIPRRPLISGILISEGLTDNTLRNVSFERVGCPACRDRRTELNEERKKGGEHPTSP",
                                 "HELPKRALPQTSNSSNPPKKKPLQGTFYEDLIRGRERFEFMRELALLENEKDGAQAKEGS",
                                 "PGRASHSLHKSKKGSTQSRHKKLMFKTDEPDSG",};

    string line;
    int i = 0;
    ifstream myfile("min_mimic.fasta");
    if (myfile.is_open()) {
        while (getline(myfile, line)) {
            CHECK(line == expected[i]);
            ++i;
        }
        myfile.close();
    }
}