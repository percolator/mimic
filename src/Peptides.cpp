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
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <assert.h>
using namespace std;
#include "AminoAcidDist.h"
#include "Option.h"
#include "Peptides.h"

const string Peptides::proteinNamePrefix = "Random_";

AminoAcidDist Peptides::background = AminoAcidDist();

Peptides::Peptides():replaceI(false)
{
}

Peptides::~Peptides()
{
}


void Peptides::printAll() {
  unsigned int count = 0;
  vector<string> outPep(between.size());
  for(map<string,set<unsigned int> >::const_iterator it=pep2ixs.begin();it!=pep2ixs.end();it++) {
    for(set<unsigned int>::const_iterator ixIt=it->second.begin();ixIt!=it->second.end();ixIt++) {
      assert(outPep[*ixIt]=="");
      outPep[*ixIt]=it->first;
      count++;
    }
  }
  assert(count==between.size());
  ostringstream first("");
  for(size_t ix=0;ix<between.size();ix++) {
    if (between[ix][0]=='>') {
      string str=first.str();
      size_t p=0;
      while(p<str.length()) {
        cout << str.substr(p,lineLen) << endl;
        p+=lineLen;
      }
      first.str("");
      cout << between[ix] << endl;
    } else {
      first << between[ix];
    }
    first << outPep[ix];
  }
  string str=first.str();
  size_t p=0;
  while(p<str.length()) {
    cout << str.substr(p,lineLen) << endl;
    p+=lineLen;
  }
  first.str("");
}

void Peptides::addPeptide(const string& peptide,unsigned int pepNo) {
  assert(pep2ixs[peptide].count(pepNo)==0);
  pep2ixs[peptide].insert(pepNo);
  assert(pepNo+1==between.size());
}

void Peptides::cleaveProtein(string seq,unsigned int & pepNo) {
  size_t lastPos = 0,pos=0,protLen=seq.length();
  if (protLen>0 && seq[protLen-1]=='*') {
    --protLen;
  }
  for (; pos<protLen; pos++) {
    if (seq[pos] == 'K' || seq[pos] == 'R' || pos == 0) {
      addPeptide(seq.substr(lastPos,pos-lastPos),pepNo++);
      int len = 1;
      while (pos+len+1<protLen && (seq[pos+len-1]=='K' || seq[pos+len-1] == 'R')) {
        len++;
      }
      between.push_back(seq.substr(pos,len));
      lastPos = pos + len;
      pos += len - 1;
    }
  }
  addPeptide(seq.substr(lastPos,pos-lastPos),pepNo++);

}


void Peptides::readFasta(string path) {
  ifstream fastaIn(path.c_str(),ios::in);
  string line(""),pep,seq("");
  unsigned int pepNo=0,proteinNo=0;
  bool spillOver=false;
  while (getline(fastaIn,line)) {
    if (line[0] == '>') {
      // id line
      if(spillOver) {
        assert(seq.size()>0);
        cleaveProtein(seq,pepNo);
        seq = "";
      }
      ostringstream newProteinName("");
      newProteinName << ">" << proteinNamePrefix << ++proteinNo;
      between.push_back(newProteinName.str());
      assert(between.size()==pepNo+1);
    } else {
      // amino acid sequence
      seq += line;
      spillOver = true;
    }
  } 
  if(spillOver) {
    cleaveProtein(seq,pepNo);
  }
}

void Peptides::shuffle(const string& in,string& out) {
  out=in;
  int i = in.length();
  if (i==0) return;
  // fisher_yates_shuffle
  for (; --i; ) {
    double d=((double)(i+1)*((double)rand()/((double)RAND_MAX+(double)1)));
    int j=(int)d;
    if (i!=j){
      char a = out[i];
      out[i] = out[j];
      out[j] = a;
    }
  }
}
void Peptides::mutate(const string& in,string& out) {
  out=in;
  int i = in.length();
  int j=(int)((double)(i+1)*((double)rand()/((double)RAND_MAX+(double)1)));
  double d=((double)rand()/((double)RAND_MAX+(double)1));
  out[j] = background.generateAA(d);
}

bool Peptides::recurringPeptide(const string& pep, bool force) {
  string checkPep=pep;
  if (replaceI) {
    for(string::iterator it=checkPep.begin();it!=checkPep.end();it++) {
      if (*it=='I')
        *it='L';
      if (*it=='N')
        *it='L';
      if (*it=='Q')
        *it='E';
      if (*it=='K')
        *it='E';
    }
  }
  bool used = usedPeptides.count(checkPep) > 0;
  if (force || !used)
    usedPeptides.insert(checkPep);
  return used;
}


void Peptides::shuffle(const Peptides& normals) {
  between=normals.between;
  for(map<string,set<unsigned int> >::const_iterator it = normals.pep2ixs.begin();it!=normals.pep2ixs.end();it++) {
    recurringPeptide(it->first,true);
  }
  for(map<string,set<unsigned int> >::const_iterator it = normals.pep2ixs.begin();it!=normals.pep2ixs.end();it++) {
    size_t tries=0;
    string scram;
    bool again;
    do {
      shuffle(it->first,scram);
      again = recurringPeptide(scram);
    } while ( again && 
             (it->first).length()>minLen &&  
             ++tries < maxTries);
    if (tries==maxTries) {
      tries = 0;
      scram = it->first;
      string mutated;
      do {
        mutate(scram,mutated);
        again = recurringPeptide(mutated);
        scram = mutated;
      } while ( again &&  ++tries < maxTries);
    }
    if (tries==maxTries) cerr << "Gave up on peptide " << it->first << endl;
    pep2ixs[scram].insert(it->second.begin(),it->second.end());
  }
}

int Peptides::run() {
  srand(time(0));
  readFasta(inFile);
  Peptides fake;
  fake.shuffle(*this);
  fake.printAll();
  return 0;
}

bool Peptides::parseOptions(int argc, char **argv){
  ostringstream intro;
  intro << "Usage:" << endl;
  intro << "   mimic <fasta-file> [<fasta-file 2> [<fasta-file 3>[...]]]" << endl;
  CommandLineParser cmd(intro.str());
  cmd.parseArgs(argc, argv);
  if (cmd.arguments.size()>0)
     inFile = cmd.arguments[0];
  else
     return false;
  return true;
}


int main(int argc, char **argv){
  Peptides *pPeptides=new Peptides();
  int retVal=-1;
  if(pPeptides->parseOptions(argc,argv))
  {
    retVal=pPeptides->run();
  }
  delete pPeptides;
  return retVal;
}   


