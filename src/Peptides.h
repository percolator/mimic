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

class Peptides
{
  public:
    Peptides();
    virtual ~Peptides();
    
    bool parseOptions(int argc, char **argv);
    int run();
    
    void readFasta(string path);
    void printAll(const vector<string>& connectorStrings, 
                  const std::string& suffix = "");
    
    void addPeptide(const string& peptide, unsigned int pepNo);
    void shuffle(const map<string,set<unsigned int> >& normalPep2ixs);
    void shuffle(const string& in, string& out);
    void mutate(const string& in, string& out);
    bool checkAndMarkUsedPeptide(const string& pep, bool force=false);
    
    void cleaveProtein(string seq, unsigned int & pepNo);
  protected:
    map<string,set<unsigned int> > pep2ixs;
    vector<string> connectorStrings_;
    set<string> usedPeptides_;
    string inFile;
    bool replaceI;
    static const unsigned int lineLen = 60;
    static const unsigned int maxTries = 1000;
    static unsigned int minLen;
    static string proteinNamePrefix_;
    static unsigned int multFactor_;
    static double sharedPeptideRatio_;
    static AminoAcidDist background;
};

#endif /*PEPTIDES_H_*/
