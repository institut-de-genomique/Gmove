/*******************************************************************************
+
+  GffRecord.cpp
+
+  Copyright (c) 2017 Genoscope, CEA, CNS, Evry, France
+  Author : Marion Dubarry, jmaury@genoscope.cns.fr
+
 *******************************************************************************/

#include "GffRecord.h"

using namespace std;

GffRecord::GffRecord(string line, s8 datatype) {
  string tmp;
  stringstream ss,ss2;
  vector<string> tmpV;
  ss<< line ;
  while(std::getline(ss, tmp, '\t')) {
    tmpV.push_back(tmp);
  }
  
  _seqname = tmpV[0];
  _method = tmpV[1];
  _type = tmpV[2];
  _start = atoi(tmpV[3].c_str());
  _end = atoi(tmpV[4].c_str());
  _score = atoi(tmpV[5].c_str());
  _strand = tmpV[6];
  _phase = tmpV[7];
  _attribute = tmpV[8];
  prepareAttribute(_attribute);
  _color=0;
  _datatype = datatype;
}

void GffRecord::prepareAttribute(string& attribute){
  std::replace( attribute.begin(), attribute.end(), '=', ' ');
  string delimiter;
  size_t pos = 0;
  if ((pos = attribute.find(' ')) != string::npos){
    delimiter='=';
    s32 start= pos + delimiter.length();
    attribute.erase(0, start);
  }
  else if((pos = attribute.find("Parent")) != string::npos){
    delimiter = "Parent";
    s32 start= pos + delimiter.length()+1;
    attribute.erase(0, start);
  }
  else if ((pos = attribute.find("ID")) != string::npos){
    delimiter="ID";
    s32 start= pos + delimiter.length()+1;
    attribute.erase(0, start);
  }
  
  //remove all ';' in string
  std::replace( attribute.begin(), attribute.end(), ';', ' ');
  //clean after if there is still some spaces
  if((pos = attribute.find(' ')) != string::npos){
    attribute.erase(pos,attribute.length());
  }
}

// to print the GffRecord
ostream& operator<<(ostream& ostr, const GffRecord& r) {
  ostr << r.getSeqName()<< "\t" << r.getMethod() << "\t";
  ostr << r.getType() << "\t" << r.getStart() << "\t" << r.getEnd() << "\t" << r.getScore() << "\t";
  ostr << r.getStrand() << "\t" << r.getPhase() << "\t" << r.getAttribute() ;
  return ostr;
}

