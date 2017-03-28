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
using namespace seqan;

// class constructor 
/*
SSRContigLists::SSRContigLists(char* covfilename,char* juncfilename, map<string, string>& sequences) {
  string seqname,strd;
  s32 start, end, startPrec, endPrec,startNext,endNext,startPrevious,endPrevious;
  f8 coverage;
  s32 id=0;
  Tstrand tmpStrand, strand;
  fstream fstrm;
  fstrm.open(covfilename, ios_base::in | ios_base::binary);
  while(!fstrm.eof()) {
    seqname.erase();
    fstrm>>seqname;
    if(seqname.empty()) { break; }
    fstrm>>start;
    fstrm>>end;
    fstrm>>coverage;
    fstrm>>tmpStrand;

    if(tmpStrand == '+') strand = SSRContig::FORWARD;
    else if(tmpStrand == '-')  	strand= SSRContig::REVERSE;
    else strand=SSRContig::UNKNOWN;

    map<string, string>::iterator itSeq = sequences.find(seqname);
    if (itSeq == sequences.end()) {
      continue;
    }
    string* seq = &itSeq->second;
    TcontigLists::iterator it = _contigs.find(seqname);
    if (it == _contigs.end()) {
      SSRContigList* liste = new SSRContigList();
      it = _contigs.insert(make_pair(seqname, liste)).first;
      startPrec=0;
      endPrec=0;
    }
    SSRContig* ctg = new SSRContig(seqname, start, end, coverage, seq,strand);
    ctg->setID(id);
    cout << *ctg<< endl;
    _nbContigs++;
    (it->second)->push_back(ctg); // it->second = TSSRList _contigs
  }
  for( TcontigLists::iterator  it=_contigs.begin(); it != _contigs.end(); it++ ) (it->second)->startposSort();
  fstrm.close();
  
  //* Known Junctions	*//*
  map<string, s32>::iterator itJS;
  cout << "load junction " <<endl;
  fstrm.open(juncfilename, ios_base::in | ios_base::binary);
  int nbKnownJunctions = 0;
  char ch[2048];
  while(!fstrm.eof()) {
	  cout << "while " << endl;
    fstrm>>seqname;
    if(seqname.empty()) { break; }
    fstrm>>start;
    fstrm>>end;
    fstrm>>strd;
    fstrm.getline(ch, 2048);
    ostringstream oss,oss2;
    if( strd == "0"){
    	oss << seqname << "@" << start << "@" << end << "@" << '1';
		string key = oss.str();
		itJS = _kwJunctions.find(key);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key, 1) );
		  nbKnownJunctions++;
		}
		oss2 << seqname << "@" << start << "@" << end << "@" << '-1';
		string key2 = oss2.str();
		itJS = _kwJunctions.find(key2);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key2, 1) );
		  nbKnownJunctions++;
		}
    }
    else{
    	oss << seqname << "@" << start << "@" << end << "@" << strd;
    	string key = oss.str();
		itJS = _kwJunctions.find(key);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key, 1) );
		  nbKnownJunctions++;
		}
    }
  }
  if(nbKnownJunctions != 0 && SSRContigList::VERBOSE) { cerr << nbKnownJunctions << " different junction(s) loaded from file [" << juncfilename << "]." << endl; }
  fstrm.close();
}
*/
// class constructor
SSRContigLists::SSRContigLists(list<pair< GffRecord, string> > listRecord, map<string, string>& sequences) {
  string seqname, tagPrevious, exonIdPrevious;
  s32 start, end, tmpEnd;
  map<string,SSRContig*> seenCovtigs;
  _nbContigs = 0;
  f8 coverage;
  s32 id=0; //TODO why the id is always 0 ?
  s32 idTranscrit = 0;
  s32 startPrevious = 0, endPrevious = 0;
  Tstrand tmpStrand, strand;
  for(list<pair<GffRecord,string > >::iterator itRecord = listRecord.begin() ; itRecord != listRecord.end(); ++itRecord){
	seqname.erase();
    assign(seqname,itRecord->first.ref);
    if(seqname.empty()) { break; }
    start = itRecord->first.beginPos + 1;
    end = itRecord->first.endPos;
    coverage = 7;//TODO change coverage of the exon ?
    tmpStrand = itRecord->first.strand;
    if(tmpStrand == '+') {strand = SSRContig::FORWARD; }
    else if(tmpStrand == '-') { 	strand= SSRContig::REVERSE;}
    else { strand=SSRContig::UNKNOWN;}
    map<string, string>::iterator itSeq = sequences.find(seqname);
    if (itSeq == sequences.end()) {
      cerr << "[SSRContigList] Can't find sequence " << seqname << " in map structure." << endl;
      exit(2);
    }
    string* seq = &itSeq->second;
    TcontigLists::iterator itContigs = _contigs.find(seqname);
    if (itContigs == _contigs.end()) {
      SSRContigList* liste = new SSRContigList();
      itContigs = _contigs.insert(make_pair(seqname, liste)).first;
    }
    //TODO else increase coverage exon

    SSRContig* ctg = new SSRContig(seqname, start, end, coverage, seq,strand);
    ctg->setID(id);
    ctg->setIdTranscrit(idTranscrit);
    ostringstream oss;
    oss << seqname << "@"<< start << "@"<<end << "@"<<strand;
    string key = oss.str();
    if(seenCovtigs.find(key) == seenCovtigs.end()){// The covtig already exist ?
    	_nbContigs++;
    	if(end-start+1 > SSRContig::MINSIZEEXON){
    		(itContigs->second)->push_back(ctg); // it->second = TSSRList _contigs
    	}
    	seenCovtigs.insert(make_pair(key,((itContigs->second)->back())));
    }
    else
    	(seenCovtigs.find(key)->second)->setIdTranscrit(idTranscrit);

  /*
   * Load Junctions
   */

    map<string, s32>::iterator itJS; //* Known Junctions	*//
  int nbKnownJunctions = 0;
  tmpEnd = end;
  if(itRecord->first.tagValues[0] == exonIdPrevious){
    ostringstream oss,oss2;

    //Different junctions if exon at the beginnning or end
    if( strand == SSRContig::UNKNOWN){//If we do not know the strand : junctions in both strand
    	oss << seqname << "@" << startPrevious<<"@"<<endPrevious << "@"<<start<<"@" << tmpEnd << "@" << 1;
		string key = oss.str();
		itJS = _kwJunctions.find(key);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key, 1) ); // key, coverage
		  nbKnownJunctions++; //XXX nbKnownJunction, what is it for ?
		}
		else
			itJS->second += 1 ;
		oss2 << seqname << "@"  << startPrevious<<"@"<<endPrevious << "@"<<start<<"@" << tmpEnd<<"@" << -1;
		string key2 = oss2.str();
		itJS = _kwJunctions.find(key2);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key2, 1) );
		  nbKnownJunctions++;
		}
		else
			itJS->second += 1 ;
    }
    else{
    	int tmp = strand;
    	oss << seqname << "@" << startPrevious<<"@"<<endPrevious << "@"<<start<<"@" << tmpEnd<<"@"<< tmp;
    	string key = oss.str();

		itJS = _kwJunctions.find(key);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key, 1) );
		  nbKnownJunctions++;
		}
		else
			itJS->second += 1 ;
    }
  }
  else if (!exonIdPrevious.empty()){
	  ++idTranscrit;
  }
  //update previous value before relooping
	startPrevious = start;
	endPrevious = tmpEnd;
	tagPrevious = itRecord->second;
	assign(exonIdPrevious,itRecord->first.tagValues[0]);
  }

  for( TcontigLists::iterator  it=_contigs.begin(); it != _contigs.end(); it++ ) (it->second)->startposSort();
  listRecord.clear();
 }


// to test junctions using the dictionary and the junctions file given as input
list<NetEx*> SSRContigLists::testJunctions(DnaDictionary& dict) {
  list<NetEx*> NetList;
  for( TcontigLists::iterator it =_contigs.begin(); it != _contigs.end(); it++ ){
    NetList.push_back( (it->second)->testJunctions(dict, _kwJunctions) );
  }
  return NetList;
}

// to extend covtigs using the dictionary
void SSRContigLists::extendCovtigs(DnaDictionary& dict) {
  TcontigLists::iterator it;
  for( it=_contigs.begin(); it != _contigs.end(); it++ ) (it->second)->extendCovtigs(dict);
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
