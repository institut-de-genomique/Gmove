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

GffRecord::GffRecord(string line) {
	  string tmp;
	  stringstream ss,ss2;
	  vector<string> tmpV;
	//  checkFormatGff(fileName );
	//  fstream file;
	//  file.open(fileName, ios::in);
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

 }



void GffRecord::prepareAttribute(string& attribute){
	//TODO prepare attribute. Extract ID different if gff2 or gff3
	// detecter si on est en gff2 ou gff3
	// detecter si tag ID Parent est présent sinon récupérer le premier Id

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

	std::replace( attribute.begin(), attribute.end(), ';', ' '); //remove all ';' in string
	//clean after if there is still some spaces
	if((pos = attribute.find(' ')) != string::npos){
		attribute.erase(pos,attribute.length());
	}
}

/*
void GffRecord::split(string s, string delimiter){
	size_t pos=0;
	list <string> tokens;
	pos = s.find(delimiter);
	//while ((pos = s.find(delimiter)) != string::npos) {
		s32 start= pos + delimiter.length();

				tokens.push_back(s.substr(start, s.length()));
	//			line.erase(0, pos + delimiter.length());
	//		}
	for(list<string>::iterator itList = tokens.begin();itList!= tokens.end();++itList){
		cout << "tokens " << *itList << endl;
	}
}
*/

// to print the GffRecord
ostream& operator<<(ostream& ostr, const GffRecord& r) {
    ostr << r.getSeqName()<< "\t" << r.getMethod() << "\t";
    ostr << r.getType() << "\t" << r.getStart() << "\t" << r.getEnd() << "\t" << r.getScore() << "\t";
    ostr << r.getStrand() << "\t" << r.getPhase() << "\t" << r.getAttribute() ;
  return ostr;
}

