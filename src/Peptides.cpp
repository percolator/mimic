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
#include <map>
#include <set>
#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <numeric>
#include <random>
#include <algorithm>
using namespace std;
#include "AminoAcidDist.h"
#include "Option.h"
#include "Peptides.h"

std::vector<size_t> getNotILIdxs(const string& sequence){
    std::vector<std::size_t> idxs;
    idxs.reserve(sequence.length());
    for (std::size_t i = 0; i < sequence.length(); ++i){
        if(sequence[i] != 'I' && sequence[i] != 'L') {
            idxs.push_back(i);
        }
    }
    return idxs;
}
std::vector<size_t> getIdxs(const string& sequence){
    std::vector<std::size_t> idxs;
    idxs.reserve(sequence.length());
    for (std::size_t i = 0; i < sequence.length(); ++i){
        idxs.push_back(i);
    }
    return idxs;
}

Peptides::Peptides(unsigned int minLength, set<string> usedPeptides, AminoAcidDist background, std::mt19937 rGen, unsigned int maxTries, bool replaceI) :
        maxTries_{maxTries}, usedPeptides_{std::move(usedPeptides)}, minLen_(minLength), replaceI_(replaceI), seed_(1u),
        multFactor_(1u), sharedPeptideRatio_(0.0),
        proteinNamePrefix_("mimic|Random_"), prependOriginal_(false),
        background_{std::move(background)}, inferAAFrequency_{false}, rGen_{rGen} {}

void Peptides::printAll(const vector<string>& connectorStrings,
    const std::string& suffix, std::ostream& os) {
  unsigned int count = 0;
  vector<string> outPep(connectorStrings.size());
  auto it = pep2ixs_.begin();
  for (; it != pep2ixs_.end(); it++) {
    auto ixIt = it->second.begin();
    for (; ixIt != it->second.end(); ixIt++) {
      assert(outPep[*ixIt].empty());
      outPep[*ixIt] = it->first;
      count++;
    }
  }
  assert(count == connectorStrings.size());
  ostringstream first("");
  for (size_t ix = 0; ix < connectorStrings.size(); ix++) {
    if (connectorStrings[ix][0] == '>') {
      string str = first.str();
      size_t p = 0;
      while (p < str.length()) {
        os << str.substr(p, lineLen_) << endl;
        p+=lineLen_;
      }
      first.str("");
      os << connectorStrings[ix] << suffix << endl;
    } else {
      first << connectorStrings[ix];
    }
    first << outPep[ix];
  }
  string str=first.str();
  size_t p=0;
  while(p<str.length()) {
    os << str.substr(p, lineLen_) << endl;
    p+=lineLen_;
  }
  first.str("");
}

void Peptides::addPeptide(const string& peptide, unsigned int pepNo) {
  for(char aa: peptide){
      absBackground_.add(aa);
  }
  assert(pep2ixs_[peptide].count(pepNo) == 0);
  pep2ixs_[peptide].insert(pepNo);
  assert(pepNo+1 == connectorStrings_.size());
}

void Peptides::cleaveProtein(string seq, unsigned int& pepNo) {
  size_t lastPos = 0, pos = 0, protLen = seq.length();
  if (protLen > 0 && seq[protLen-1]=='*') {
    --protLen;
  }
  for (; pos<protLen; pos++) {
    if (pos == 0 || seq[pos] == 'K' || seq[pos] == 'R' || pos==protLen-1) {
      // store peptide without C-terminal 'K/R'
        addPeptide(seq.substr(lastPos,pos-lastPos), pepNo++);

      // store single/multiple consecutive 'K/R' to leave unchanged in
      // shuffled protein sequence
      int len = 1;
      while (pos+len+1 < protLen && (seq[pos+len-1] == 'K' || seq[pos+len-1] == 'R')) {
        len++;
      }
      connectorStrings_.push_back(seq.substr(pos,len));

      // update position to jump over 'K/R' region
      lastPos = pos + len;
      pos += len - 1;
    }
  }

  // store protein C-terminal peptide
    addPeptide(seq.substr(lastPos,pos-lastPos), pepNo++);
}


void Peptides::readFasta(string& path, bool write, std::ostream& os) {
  ifstream fastaIn(path.c_str(), ios::in);
  string line, pep, seq;
  unsigned int pepNo = 0,proteinNo = 0;
  bool spillOver = false;
  while (getline(fastaIn, line)) {
    if (write) os << line << std::endl;

    if (line[0] == '>') { // id line
      if (spillOver) {
        assert(!seq.empty());
        // add peptides and KR-stretches to connectorStrings_
        cleaveProtein(seq, pepNo);
        seq = "";
      }
      ostringstream newProteinName("");
      newProteinName << ">" << proteinNamePrefix_ << ++proteinNo;
      connectorStrings_.push_back(newProteinName.str());
      assert(connectorStrings_.size() == pepNo+1);
    } else { // amino acid sequence
      seq += line;
      spillOver = true;
    }
  }
  if (spillOver) {
    cleaveProtein(seq, pepNo);
  }
  if(inferAAFrequency_){
      background_.setDist(absBackground_.getDist(), replaceI_);
  }
}

void Peptides::shuffle(const string& in,string& out) {
  out=in;
  std::shuffle(out.begin(), out.end(), rGen_);
}
void Peptides::mutate(const string& in, string& out) {
  out=in;
    std::vector<size_t> okIds;
  if(replaceI_) {
      okIds = getNotILIdxs(in);
  } else {
      okIds = getIdxs(in);
  }
  std::size_t i = okIds.size();
  if (i==0){
      return;
  }
  std::uniform_int_distribution<size_t> distrib(0, i-1);
  double d= uniformDist();
  auto mutateId = okIds[distrib(rGen_)];
  out[mutateId] = background_.generateAA(d);
}

bool Peptides::checkAndMarkUsedPeptide(const string& pep, bool force) {
  string checkPep=pep;
  if (replaceI_) {
    for(char & it : checkPep) {
      if (it=='I')
        it='L';
    }
  }
  bool used = usedPeptides_.count(checkPep) > 0;
  if (force || !used)
    usedPeptides_.insert(checkPep);
  return used;
}


void Peptides::shuffle(const map<string,set<unsigned int> >& normalPep2ixs, std::ofstream &logger) {
  auto it = normalPep2ixs.begin();
  bool force = true;
  for(; it != normalPep2ixs.end(); it++) {
    checkAndMarkUsedPeptide(it->first, force);
  }
  it = normalPep2ixs.begin();
  for (; it != normalPep2ixs.end(); it++) {
    double uniRand = uniformDist();
    string scrambledPeptide = it->first;
    bool wasMutated = false;
    if (uniRand >= sharedPeptideRatio_ && scrambledPeptide.length() > 0) {
      size_t tries = 0;
      bool peptideUsed;
      do {
        shuffle(it->first, scrambledPeptide);
        peptideUsed = checkAndMarkUsedPeptide(scrambledPeptide);
      } while ( peptideUsed &&
               (it->first).length() >= minLen_ &&
                ++tries < maxTries_);

      // if scrambling does not give a usable peptide, mutate some AAs
      if (tries == maxTries_) {
        tries = 0;
        scrambledPeptide = it->first;
        string mutatedPeptide;
        do {
          mutate(scrambledPeptide, mutatedPeptide);
          wasMutated = true;
          peptideUsed = checkAndMarkUsedPeptide(mutatedPeptide);
          scrambledPeptide = mutatedPeptide;
        } while ( peptideUsed && ++tries < maxTries_);
      }
      //if (tries == maxTries_) cerr << "Gave up on peptide " << it->first << endl;
    }
    if(isVerbose_ && scrambledPeptide.length() > minLen_){
        logger << it->first << ',' << scrambledPeptide << ',' << scrambledPeptide.length() << ',' << wasMutated << '\n';
    }
    pep2ixs_[scrambledPeptide].insert(it->second.begin(),it->second.end());
  }
}

int Peptides::run() {
  srand(seed_);
  rGen_.seed(seed_);

  std::ofstream outFileStream;
  if (!outFile_.empty()) {
    outFileStream.open(outFile_.c_str(), ios::out);
  }
  std::ostream &outStream = !outFile_.empty() ? outFileStream : std::cout;

  std::ofstream logger("shuffle.log", ios::trunc);
  logger << "original" << ',' << "scrambled" << ',' << "length" << ',' << "wasMutated" << '\n';

  cerr << "Reading fasta file and in-silico digesting proteins" << endl;
  readFasta(inFile_, prependOriginal_, outStream);
  background_.print(cerr);
  for (unsigned int m = 0; m < multFactor_; ++m) {
    Peptides entrapmentDB(minLen_, usedPeptides_, background_, rGen_, maxTries_, replaceI_);

    cerr << "Shuffling round: " << (m+1) << endl;
    entrapmentDB.shuffle(pep2ixs_, logger);

    ostringstream suffix("");
    if (multFactor_ > 1) {
      suffix << "|shuffle_" << (m+1);
    }

    entrapmentDB.printAll(connectorStrings_, suffix.str(), outStream);
    usedPeptides_ = entrapmentDB.usedPeptides_;
      rGen_ = entrapmentDB.rGen_;
  }
  return 0;
}

bool Peptides::parseOptions(int argc, char **argv){
  ostringstream intro;
  intro << "Usage:" << endl;
  intro << "   mimic <fasta-file>" << endl;
  CommandLineParser cmd(intro.str());

  cmd.defineOption("o",
      "out-file",
      "File to write the mimicked protein sequences to (Default: prints to stdout)",
      "string");
  cmd.defineOption("p",
      "prefix",
      "Prefix to mimic proteins (Default: \"mimic|Random_\")",
      "string");
  cmd.defineOption("m",
      "mult-factor",
      "Number of times the database should be multiplied (Default: 1)",
      "int");
  cmd.defineOption("s",
      "shared-pept-ratio",
      "Ratio of shared peptides that will stay preserved in the mimic database (Default: 0.0)",
      "int");
  cmd.defineOption("S",
      "seed",
      "Set seed of the random number generator. Default = 1",
      "int");
  cmd.defineOption("P",
      "prepend",
      "Prepend the original fasta file to the output",
      "",
      TRUE_IF_SET);
    cmd.defineOption("I",
                     "replaceI",
                     "New sequences are checked if they were already generated. If this flag is added all I will count as L for this check",
                     "",
                     TRUE_IF_SET);
    cmd.defineOption("e",
                     "empiric",
                     "Instead of assuming the AA frequency infer it from the given fasta file",
                     "",
                     TRUE_IF_SET);

  cmd.parseArgs(argc, argv);

  if (cmd.optionSet("o")) {
    outFile_ = cmd.options["o"];
  }
  if (cmd.optionSet("m")) {
    multFactor_ = cmd.getInt("m", 1, 1000);
  }
  if (cmd.optionSet("s")) {
    sharedPeptideRatio_ = cmd.getDouble("s", 0.0, 1.0);
  }
  if (cmd.optionSet("S")) {
    seed_ = cmd.getInt("S", 1, 2147483647);
  }
  if (cmd.optionSet("p")) {
    proteinNamePrefix_ = cmd.options["p"];
  }
  if (cmd.optionSet("P")) {
    prependOriginal_ = true;
  }
  if (cmd.optionSet("I")) {
    replaceI_ = true;
    background_ = AminoAcidDist(true);
  }
  if (cmd.optionSet("e")) {
      inferAAFrequency_ = true;
  }

  if (!cmd.arguments.empty()) {
    inFile_ = cmd.arguments[0];
  } else {
    std::cerr << "No fasta file specified" << std::endl;
    return false;
  }
  return true;
}

Peptides::Peptides(unsigned int minLen, set<string> usedPeptides, std::mt19937 rGen, unsigned int maxTries, bool replaceI):
Peptides(minLen, std::move(usedPeptides), AminoAcidDist{replaceI}, rGen, maxTries, replaceI) {

}

Peptides::Peptides() : Peptides(4, {}, AminoAcidDist{false}, {}, 1000, false) {

}

double Peptides::uniformDist(double a, double b) {
    auto dist = std::uniform_real_distribution(a,b);
    return dist(rGen_);
}





