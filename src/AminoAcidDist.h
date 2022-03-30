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

#ifndef AMINOACIDDIST_H_
#define AMINOACIDDIST_H_
#include <map>
#include <iostream>

using namespace std;

class AminoAcidDist
{
public:
	explicit AminoAcidDist(bool removeIl);
    char generateAA(double p);
    static void normalize(map<char,double> &dist);
protected:
    map<char,double> dist_;
public:
    [[nodiscard]] const map<char, double> &getDist() const;

public:
    void setDist(map<char, double> dist, bool removeIL);
    void print(std::ostream& os);
};

class AbsAminoAcidDist
{
public:
    explicit AbsAminoAcidDist() = default;
    void add(char AA);
    inline static const string validAAs = "ACDEFGHIKLMNOPQRSTUVWY";
protected:
    map<char,int> dist_;
public:
    [[nodiscard]] map<char, double> getDist() const;
};

#endif /*AMINOACIDDIST_H_*/
