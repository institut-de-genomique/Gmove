/*******************************************************************************
+
+  SSRContigLists.cpp
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
 *******************************************************************************/

#include "SSRContigLists.h"

using namespace std;

// class constructor 
SSRContigLists::SSRContigLists(GffRecordList listRecord, map<string, string>& sequences, map<string, s32>& badintrons) {
  string seqname, tagPrevious, exonIdPrevious;
  s32 start, end, tmpEnd;
  map<string,SSRContig*> seenCovtigs;
  f8 coverage;
  s32 startPrevious = 0, endPrevious = 0;
  Tstrand strand;
  map<s32,TSSRList*> mapTranscrit;
  s32 removed = 0;
  
  GffRecordL* allRecords = listRecord.getRecords();
  
  for(GffRecordL::iterator itRecord = allRecords->begin() ; itRecord != allRecords->end(); ++itRecord) {
    
    GffRecord* record = *itRecord;
    seqname = record->getSeqName();
    string currentId = record->getAttribute();
    
    if(seqname.empty()) { break; }
    
    start = record->getStart(); //+ 1;
    end = record->getEnd();
    coverage = 0;
    if(record->getStrand() == "+") strand = SSRContig::FORWARD; 
    else if(record->getStrand() == "-") strand = SSRContig::REVERSE;
    else  strand = SSRContig::UNKNOWN;
    map<string, string>::iterator itSeq = sequences.find(seqname);
    if (itSeq == sequences.end()) {
      cerr << "[SSRContigList] Can't find sequence " << seqname << " in map structure." << endl;
      exit(2);
    }
    string* seq = &itSeq->second;
    
    SSRContig* ctg = new SSRContig(seqname, start, end, coverage, seq, strand);
    s32 color = record->getColor();
    ctg->setIdTranscrit(color);
    ostringstream oss;
    oss << seqname << "@"<< start << "@"<<end << "@"<<strand;
    string key = oss.str();
    map<string,SSRContig*>::iterator itSeenCovtigs = seenCovtigs.find(key);
    if( itSeenCovtigs == seenCovtigs.end()){// The covtig doesnt exist
      if(end-start+1 > SSRContig::MINSIZEEXON) {
	TcontigLists::iterator itContigs = _contigs.find(seqname);
	if (itContigs == _contigs.end()) {
	  SSRContigList* liste = new SSRContigList();
	  itContigs = _contigs.insert(make_pair(seqname, liste)).first;
	}
	(itContigs->second)->push_back(ctg); // it->second = TSSRList _contigs
	map<s32,list<string> >::iterator itTranscrit = _transcrit.find(color);
	if(itTranscrit != _transcrit.end()){
	  itTranscrit->second.push_back(key);
	}
	else{
	  list<string> tmpList;
	  tmpList.push_back(key);
	  _transcrit.insert(make_pair(color,tmpList));
	}
	
	seenCovtigs.insert(make_pair(key,((itContigs->second)->back())));
      }
      else delete(ctg);      
    }
    else {
      (itSeenCovtigs->second)->setIdTranscrit(color);
      map<s32,list<string> >::iterator itTranscrit = _transcrit.find(color);
      if(itTranscrit != _transcrit.end()){
	itTranscrit->second.push_back(key);
      }
      else{
	list<string> tmpList;
	tmpList.push_back(key);
	_transcrit.insert(make_pair(color,tmpList));
      }
      delete(ctg);
    }
    /*
     * Load Junctions
     */
    
    map<string, s32>::iterator itJS; // Known Junctions
    tmpEnd = end;
    if(record->getAttribute() == exonIdPrevious){
      if(strand == SSRContig::UNKNOWN) {
	removed += addJunction(seqname, startPrevious, endPrevious, start, tmpEnd, 1, record->getDatatype(), badintrons);
	removed += addJunction(seqname, startPrevious, endPrevious, start, tmpEnd, -1, record->getDatatype(), badintrons);
      } else {
	removed += addJunction(seqname, startPrevious, endPrevious, start, tmpEnd, (s32)strand, record->getDatatype(), badintrons);
      }
    }
    //update previous value before relooping
    startPrevious = start;
    endPrevious = tmpEnd;
    exonIdPrevious = record->getAttribute();
  }
  cout << removed << " introns from GFF files were removed thanks to the exclude_introns file." << endl << endl;
  for( TcontigLists::iterator  it=_contigs.begin(); it != _contigs.end(); it++ ) (it->second)->startposSort();
}

s32 SSRContigLists::addJunction(string seqname, s32 start1, s32 end1, s32 start2, 
				s32 end2, s32 strand, s8 datatype, map<string, s32>& badintrons) {
  ostringstream ossintron, oss;
  ossintron << seqname << "@" << (end1 + 1) << "@" << (start2 - 1) << "@" << strand;
  string keyI = ossintron.str();
  map<string, s32>::iterator it = badintrons.find(keyI);
  if (it == badintrons.end()) {
    oss << seqname << "@" << start1 << "@" << end1 << "@" << start2 << "@" << end2 << "@" << strand;
    string key = oss.str();
    it = _kwJunctions.find(key);
    if (it == _kwJunctions.end()) _kwJunctions.insert(make_pair(key, datatype));
    else it->second += 1;
    return 0;
  } else return 1;
}

// to test junctions using the dictionary and the junctions file given as input
list<NetEx*> SSRContigLists::buildGraph() {
  list<NetEx*> NetList;
  for( TcontigLists::iterator it =_contigs.begin(); it != _contigs.end(); it++ ){//loop on transcrits by scaffold
    NetList.push_back( (it->second)->buildGraph( _kwJunctions) );
  }
  return NetList;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}


void SSRContigLists::cleanJunctions(){
	//Delete some jonctions based on length and coverage
	map<string,pair<s32,list<string> > >::iterator itJ;
	string key;
	map<string,pair<s32,list<string> > > newJunctions;

for(map<string,s32>::iterator itKJ = _kwJunctions.begin(); itKJ != _kwJunctions.end(); ++itKJ){
	vector<string> blop = split(itKJ->first,'@');
	s32 startJunction =  atoi(blop[2].c_str());
	s32 endJunction =  atoi(blop[3].c_str());
	s32 sizeJunction = endJunction-startJunction+1;

	if(sizeJunction > 5000 ){
		key = blop[2];
		key += "@";
		key += blop[3];
		 itJ = newJunctions.find(key);
		if(itJ == newJunctions.end()){
			list<string> tmpList;
			tmpList.push_back(itKJ->first);
			pair<s32,list<string> > tmpP =  make_pair(itKJ->second,tmpList);
			pair<string,pair<s32,list<string> > > tmpPP = make_pair(key,tmpP);
			newJunctions.insert(tmpPP);
		}
		else{
			itJ->second.first += itKJ->second;
			itJ->second.second.push_back(itKJ->first);
		}
	}
}
	for(map<string,pair<s32,list<string> > >::iterator itNJ =newJunctions.begin();  itNJ != newJunctions.end(); ++itNJ){
		vector<string> blop = split(itNJ->first,'@');
		if(itNJ->second.first < 2 ){
			for(list<string>::iterator itList = itNJ->second.second.begin(); itList != itNJ->second.second.end(); ++itList){
				_kwJunctions.erase(*itList);
			}
		}
	}
}

// to print the covtigs
ostream& operator<<(ostream& ostr, const SSRContigLists& d) {
  TcontigLists::const_iterator it;
  for( it=(d._contigs).begin(); it != (d._contigs).end(); it++ ) ostr << *(it->second);
  return ostr;
}
