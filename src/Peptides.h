/*******************************************************************************
    Copyright 2006-2009 Lukas Käll <lukas.kall@cbr.su.se>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

 *******************************************************************************/
#ifndef PEPTIDES_H_
#define PEPTIDES_H_
#include <string>
#include <vector>
#include <map>
#include <set>
#include <random>
#include "AminoAcidDist.h"
using namespace std;

class Peptides
{
  public:
    Peptides(unsigned int minLen, set<string> usedPeptides, AminoAcidDist background, std::mt19937 rGen, unsigned int maxTries = 1000, bool replaceI=false);
    Peptides(unsigned int minLen, set<string> usedPeptides, std::mt19937 rGen, unsigned int maxTries = 1000, bool replaceI=false);
    Peptides();
    
    bool parseOptions(int argc, char **argv);
    int run();

    void readFasta(string& path, bool write, std::ostream& os);
    void printAll(const vector<string>& connectorStrings, 
                  const std::string& suffix,
                  std::ostream& os);
    
    void addPeptide(const string& peptide, unsigned int pepNo);
    void shuffle(const map<string,set<unsigned int> >& normalPep2ixs, std::ofstream &logger);
    void shuffle(const string& in, string& out);
    void mutate(const string& in, string& out);
    bool checkAndMarkUsedPeptide(const string& pep, bool force=false);
    
    void cleaveProtein(string seq, unsigned int & pepNo);
    double uniformDist(double a = 0.0, double b = 1.0);
  protected:
    static const unsigned int lineLen_ = 60;
    unsigned int maxTries_;
    
    map<string,set<unsigned int> > pep2ixs_;
    vector<string> connectorStrings_;
    set<string> usedPeptides_;
    
    // input arguments
    string inFile_, outFile_;
    unsigned int seed_;
    bool replaceI_;
    // N.B.: the shortest shuffled peptide will be minLen_ + 1, as we conserve 
    // the last AA of each peptide
    unsigned int minLen_;
    string proteinNamePrefix_;
    unsigned int multFactor_;
    double sharedPeptideRatio_;
    bool prependOriginal_;
    
    AminoAcidDist background_;
    std::mt19937 rGen_;
    bool inferAAFrequency_;
    AbsAminoAcidDist absBackground_;
    bool isVerbose_ = false;
    bool noDigest_ = false;
};

#endif /*PEPTIDES_H_*/
