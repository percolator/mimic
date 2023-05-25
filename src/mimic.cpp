//
// Created by mgraber on 21/03/2022.
//
#include "Peptides.h"

int main(int argc, char **argv){
    Peptides pPeptides;
    int retVal = -1;
    if (pPeptides.parseOptions(argc, argv)) {
        retVal = pPeptides.run();
    }
    return retVal;
}