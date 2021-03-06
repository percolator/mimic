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
using namespace std;
#include "AminoAcidDist.h"
#include <assert.h>

AminoAcidDist::AminoAcidDist()
{
  dist['A'] = 0.081;
  dist['C'] = 0.015;
  dist['D'] = 0.054;
  dist['E'] = 0.061;
  dist['F'] = 0.040;
  dist['G'] = 0.068;
  dist['H'] = 0.022;
  dist['I'] = 0.057;
  dist['K'] = 0.056;
  dist['L'] = 0.093;
  dist['M'] = 0.025;
  dist['N'] = 0.045;
  dist['P'] = 0.049;
  dist['Q'] = 0.039;
  dist['R'] = 0.057;
  dist['S'] = 0.068;
  dist['T'] = 0.058;
  dist['V'] = 0.067;
  dist['W'] = 0.013;
  dist['Y'] = 0.032;
}

AminoAcidDist::~AminoAcidDist()
{
}

char AminoAcidDist::generateAA(double p) {
  assert(p >= 0.0 && p <= 1.0);
  double cumulative = 0.0;
  map<char,double>::iterator it;
  for(it=dist.begin();it!=dist.end();it++) {
    cumulative += it->second;
    if (p<=cumulative) {
      return it->first;
    }
  }
  assert(false);
}
