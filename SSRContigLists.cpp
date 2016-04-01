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
      //cerr << "[SSRContigList] Can't find sequence " << seqname << " in map structure." << endl;
      //exit(2);
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
  
  //* Known Junctions	*//
  map<string, s32>::iterator itJS;
  
  fstrm.open(juncfilename, ios_base::in | ios_base::binary);
  int nbKnownJunctions = 0;
  char ch[2048];
  while(!fstrm.eof()) {
    fstrm>>seqname;
    if(seqname.empty()) { break; }
    map<string, string>::iterator itSeq = sequences.find(seqname);
       if (itSeq == sequences.end()) {
         continue;
         //cerr << "[SSRContigList] Can't find sequence " << seqname << " in map structure." << endl;
         //exit(2);
       }
    fstrm>>startPrevious;
    fstrm>>endPrevious;
    fstrm>>startNext;
    fstrm>>endNext;
    fstrm>>strd;
    fstrm.getline(ch, 2048);
    ostringstream oss,oss2;
    if( strd == "0"){
    	oss << seqname << "@" << startPrevious<<"@"<<endPrevious << "@"<<startNext<<"@" << endNext << "@" << '1';
		string key = oss.str();
		itJS = _kwJunctions.find(key);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key, 1) );
		  nbKnownJunctions++;
		}
		oss2 << seqname << "@"  << startPrevious<<"@"<<endPrevious << "@"<<startNext<<"@" << endNext<<"@" << '-1';
		string key2 = oss2.str();
		itJS = _kwJunctions.find(key2);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key2, 1) );
		  nbKnownJunctions++;
		}
    }
    else{
    	oss << seqname << "@" << startPrevious<<"@"<<endPrevious << "@"<<startNext<<"@" << endNext<<"@"<< strd;
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

// class constructor
SSRContigLists::SSRContigLists(list<pair< GffRecord, string> > listRecord, map<string, string>& sequences) {
  string seqname, tagPrevious, exonIdPrevious;
  s32 start, end, startPrec, endPrec,tmpEnd;
  map<string,bool> seenCovtigs;
  f8 coverage;
  s32 id=0; //TODO why the id is always 0 ?
  s32 startPrevious = 0, endPrevious = 0;
  Tstrand tmpStrand, strand;
  for(list<pair<GffRecord,string > >::iterator itRecord = listRecord.begin() ; itRecord != listRecord.end(); ++itRecord){
	seqname.erase();
    assign(seqname,itRecord->first.ref);
    if(seqname.empty()) { break; }
    start = itRecord->first.beginPos + 1;
    end = itRecord->first.endPos;
    coverage = 7;
    tmpStrand = itRecord->first.strand;
    if(tmpStrand == '+') strand = SSRContig::FORWARD;
    else if(tmpStrand == '-')  	strand= SSRContig::REVERSE;
    else strand=SSRContig::UNKNOWN;
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
      startPrec=0;
      endPrec=0;
    }
    SSRContig* ctg = new SSRContig(seqname, start, end, coverage, seq,strand);
    ctg->setID(id);

    ctg->setTag(itRecord->second);//FIXME update the tag if we find the key but the previous tag is d or f and the new tag is i
    ostringstream oss;
    oss << seqname << "@"<< start << "@"<<end << "@"<<strand;
    string key = oss.str();
    if(seenCovtigs.find(key) == seenCovtigs.end()){
    	seenCovtigs.insert(make_pair(key,true));
    	_nbContigs++;
    	(itContigs->second)->push_back(ctg); // it->second = TSSRList _contigs
    }
  //  cout << *ctg <<" " <<itRecord->second<< endl;

  /*
   * Load Junctions
   */

    map<string, s32>::iterator itJS; //* Known Junctions	*//
  int nbKnownJunctions = 0;
  tmpEnd = end;
  if(itRecord->first.tagValues[0] == exonIdPrevious){
	//  cout << "id " << itRecord->first.tagValues[0] << " "<< exonIdPrevious << endl;
    ostringstream oss,oss2;

    //Different junctions if exon at the beginnning or end
/*    if(tagPrevious == "d" || tagPrevious == "f"){
    	cout << " enter in IF "<< endl;
    	startPrevious =-1;
    	tmpEnd = -1;
    }
 */   if( strand == SSRContig::UNKNOWN){//If we do not know the strand : junctions in both strand
    	oss << seqname << "@" << startPrevious<<"@"<<endPrevious << "@"<<start<<"@" << tmpEnd << "@" << 1;
		string key = oss.str();
		itJS = _kwJunctions.find(key);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key, 1) );
		  nbKnownJunctions++; //XXX nbKnownJunction, what is it for ?
		}
		oss2 << seqname << "@"  << startPrevious<<"@"<<endPrevious << "@"<<start<<"@" << tmpEnd<<"@" << -1;
		string key2 = oss2.str();
		itJS = _kwJunctions.find(key2);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key2, 1) );
		  nbKnownJunctions++;
		}
    }
    else{
    	int tmp = strand;
    	oss << seqname << "@" << startPrevious<<"@"<<endPrevious << "@"<<start<<"@" << tmpEnd<<"@"<< tmp;
    	string key = oss.str();

  //  	cout << " jonction " << key << endl;
		itJS = _kwJunctions.find(key);
		if (itJS == _kwJunctions.end()) {
		  _kwJunctions.insert( make_pair(key, 1) );
		  nbKnownJunctions++;

		}
    }
  }
	//update previous value before relooping
	startPrevious = start;
	endPrevious = tmpEnd;
	tagPrevious = itRecord->second;
	assign(exonIdPrevious,itRecord->first.tagValues[0]);
 // if(nbKnownJunctions != 0 && SSRContigList::VERBOSE) { cerr << nbKnownJunctions << " different junction(s) loaded from file [" << juncfilename << "]." << endl; }
//  ++id;

 //   cout << " Nombre de contigs " << itContigs->second->size() << " "<< _nbContigs << " " << "nombre de jonctions "<< _kwJunctions.size() << endl;
  }
  for( TcontigLists::iterator  it=_contigs.begin(); it != _contigs.end(); it++ ) (it->second)->startposSort();
  listRecord.clear();
}


// to test junctions using the dictionary and the junctions file given as input
list<NetEx*> SSRContigLists::testJunctions(DnaDictionary& dict, ofstream& out) {

  list<NetEx*> NetList;

  for( TcontigLists::iterator it =_contigs.begin(); it != _contigs.end(); it++ ){
    NetList.push_back( (it->second)->testJunctions(dict, _kwJunctions, out) );
  }
  return NetList;
}

// to extend covtigs using the dictionary
void SSRContigLists::extendCovtigs(DnaDictionary& dict) {
  TcontigLists::iterator it;
  for( it=_contigs.begin(); it != _contigs.end(); it++ ) (it->second)->extendCovtigs(dict);
}

// to print the covtigs
ostream& operator<<(ostream& ostr, const SSRContigLists& d) {
  TcontigLists::const_iterator it;
  for( it=(d._contigs).begin(); it != (d._contigs).end(); it++ ) ostr << *(it->second);
  return ostr;
}
