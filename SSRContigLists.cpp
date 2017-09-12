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
//using namespace seqan;

// class constructor 
//TODO simplifier pas besoin de pair<GffRecord,string>
SSRContigLists::SSRContigLists(GffRecordList listRecord, map<string, string>& sequences) {
  string seqname, tagPrevious, exonIdPrevious;
  s32 start, end, tmpEnd;
  map<string,SSRContig*> seenCovtigs;
  f8 coverage;
  s32 idTranscrit = 0;
  s32 startPrevious = 0, endPrevious = 0;
  Tstrand strand;
  map<s32,TSSRList*> mapTranscrit;

  GffRecordL* allRecords = listRecord.getRecords();
  for(GffRecordL::iterator itRecord = allRecords->begin() ; itRecord != allRecords->end(); ++itRecord){

	  GffRecord* record = *itRecord;
//	seqname.erase();
    //assign(seqname,itRecord->first.ref);
	  seqname = record->getSeqName();
    string currentId = record->getAttribute();
    if(seqname.empty()) { break; }
    // assign(currentId,itRecord->first.tagValues[0]);
     if(exonIdPrevious != currentId && !exonIdPrevious.empty()){
     		  ++idTranscrit;
     	}
  //   cout <<"blop exonIdPrevious "<<  exonIdPrevious << endl;
    start = record->getStart(); //+ 1;
    end = record->getEnd();
    coverage = 0;//TODO change coverage of the exon ?
    if(record->getStrand() == "+") {strand = SSRContig::FORWARD; }
    else if(record->getStrand() == "-") { 	strand= SSRContig::REVERSE;}
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
    ctg->setIdTranscrit(idTranscrit);
    ostringstream oss;
    oss << seqname << "@"<< start << "@"<<end << "@"<<strand;
    string key = oss.str();
    if(seenCovtigs.find(key) == seenCovtigs.end()){// The covtig doesnt exist
  //  	_nbContigs++;
    	if(end-start+1 > SSRContig::MINSIZEEXON){
    		(itContigs->second)->push_back(ctg); // it->second = TSSRList _contigs
    		map<s32,list<string> >::iterator itTranscrit = _transcrit.find(idTranscrit);
    		if(itTranscrit != _transcrit.end()){
    			itTranscrit->second.push_back(key);
    		}
    		else{
    			list<string> tmpList;
    			tmpList.push_back(key);
    			_transcrit.insert(make_pair(idTranscrit,tmpList));
    		}
    	}
    	seenCovtigs.insert(make_pair(key,((itContigs->second)->back())));
    }
    else{
    	(seenCovtigs.find(key)->second)->setIdTranscrit(idTranscrit);
    	map<s32,list<string> >::iterator itTranscrit = _transcrit.find(idTranscrit);
			if(itTranscrit != _transcrit.end()){
				itTranscrit->second.push_back(key);
			}
			else{
				list<string> tmpList;
				tmpList.push_back(key);
				_transcrit.insert(make_pair(idTranscrit,tmpList));
			}
    }
  /*
   * Load Junctions
   */

    map<string, s32>::iterator itJS; // Known Junctions
  int nbKnownJunctions = 0;
  tmpEnd = end;
 // cout << record->getAttribute() << " exonIdPrevious "<< exonIdPrevious<<endl;
  if(record->getAttribute() == exonIdPrevious){
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
	//  cout << "exonIdPrevious add +1 to idTranscrit "<< idTranscrit<< endl;
	  ++idTranscrit;
  }

  //update previous value before relooping
	startPrevious = start;
	endPrevious = tmpEnd;
	exonIdPrevious = record->getAttribute();
//	tagPrevious = itRecord->second;
//	assign(exonIdPrevious,itRecord->first.tagValues[0]);

	//cout << "put id Transcrit in constructor " << idTranscrit << endl;
  }
	//XXX TEST
/*	for(map<s32,list<string> >::iterator itT = _transcrit.begin(); itT != _transcrit.end(); ++itT){
		cout << itT->first << " ";
		for(list<string>::iterator itTT = itT->second.begin() ; itTT!= itT->second.end(); ++itTT){
			cout << *itTT << " ";
		}
		cout << endl;
	}

*/
  for( TcontigLists::iterator  it=_contigs.begin(); it != _contigs.end(); it++ ) (it->second)->startposSort();
//  listRecord.clear();
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
