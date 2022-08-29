/*******************************************************************************
+
+  SSRContigList.cpp
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
 *******************************************************************************/

#include "SSRContigList.h"

using namespace std;

s32 SSRContigList::MINOVERLAPJUNCTION;
s32 SSRContigList::NBNEIGHBOUR;
s32 SSRContigList::MINSIZEINTRON;
s32 SSRContigList::MAXSIZEINTRON;
s32 SSRContigList::MINCOVWORD2ORIENTATE;
//s32 SSRContigList::REDUCEEXONLIST;
//s32 SSRContigList::VERBOSE;


// to build a list of potential exons from a covtig
void SSRContigList::populateExons(SSRContig* ctg, map<string, bool>& _seenExons) {
  ctg->populateExons();

  TSSRList* Lexons = ctg->getExons();
  for(TSSRList::iterator itEx = Lexons->begin(); itEx != Lexons->end(); itEx++) {
    SSRContig* exon = *itEx;
    ostringstream oss;
    oss << exon->seqName() << "@" << exon->start() << "@" << exon->end() << "@" << exon->strand();
    string key = oss.str();
    map<string, bool>::iterator it = _seenExons.find(key);
    if (it != _seenExons.end()) {
      itEx = Lexons->erase(itEx);
      delete exon;
      itEx--;
    } else {
      _seenExons.insert(make_pair(key, true));
    }
  }
}

map<s32,list<SSRContig*> > SSRContigList::mapCovtigs(TSSRList cov){
  //time_t  before = time(NULL);
  map<s32,list<SSRContig*> > mCov;
  map<s32,list<SSRContig*> >::iterator itMCov;
  for(TSSRList::iterator itCov = cov.begin(); itCov != cov.end();++itCov){
    itMCov = mCov.find((*itCov)->start());
    
    if(itMCov == mCov.end()){
      list<SSRContig*> tmpList;
      tmpList.push_back(*itCov);
      mCov.insert(make_pair((*itCov)->start(),tmpList));
    }
    else{
      itMCov->second.push_back(*itCov);
    }
  }
  //time_t  after = time(NULL);
  return mCov;
}

map<pair<string, s32>,list<string> > SSRContigList::mapJunctions(map<string,s32> junction){
  map<pair<string, s32>,list<string> > mJunction;
  map<pair<string,s32>,list<string> >::iterator itMJunction;
  
  for(map<string,s32>::iterator itJunction = junction.begin(); itJunction != junction.end();++itJunction){
    vector<string> splitJunction = split(itJunction->first,'@');
    s32 start = atoi(splitJunction[2].c_str());
    string name = splitJunction[0];
    itMJunction = mJunction.find(make_pair(name,start));
    if(itMJunction == mJunction.end()){
      list<string> tmpList;
      tmpList.push_back(itJunction->first);
      pair<string,s32> tmpPair = make_pair(name,start);
      mJunction.insert(make_pair(tmpPair,tmpList));
    }
    else itMJunction->second.push_back(itJunction->first);
  }
  return mJunction;
}

NetEx* SSRContigList::buildGraph( map<string, s32>& KJunctions) {
  TSSRList::iterator itPrec;
  map<string, bool> _seenExons;
  
  for( itPrec = _contigs.begin(); itPrec != _contigs.end(); itPrec++) {
    SSRContig* covtigA = *itPrec;
    if(!covtigA->isPopulate())  this->populateExons(covtigA, _seenExons);
  }
  SSRContig* covtigA = *(_contigs.begin());
  
  string seqname = covtigA->seqName();
  TSSRList* vertices = exons2vertices();
  map<Edge,s32> weight;
  list<Edge>* edges = junctions2edges(vertices, KJunctions, weight);
  NetEx* ExNtk = new NetEx(vertices, edges, seqname, weight);
  return ExNtk;
}

TSSRList* SSRContigList::exons2vertices() {
  TSSRList* vertices = new TSSRList();
  s32 i = 0;
  for(TSSRList::iterator it = _contigs.begin(); it != _contigs.end(); it++) {
    SSRContig* ctg = *it;
    TSSRList* exons = ctg->getExons();

    for(TSSRList::iterator itEx = exons->begin(); itEx != exons->end(); itEx++) {
      SSRContig* exon = *itEx;
      exon->setID(i);
      vertices->push_back( exon );
      ++i;
    }
  }
  return vertices;
}

// to add junctions between vertices in the list of graph edges
list<Edge>* SSRContigList::junctions2edges(TSSRList* vertices, map<string,s32> KJunctions,  map<Edge,s32>& weight) {
  list<Edge>* edges = new list<Edge>;
  
  for(TSSRList::iterator itV = vertices->begin(); itV != vertices->end(); itV++) {
    SSRContig* v = *itV;
    for(TSSRList::iterator itW = vertices->begin(); itW != vertices->end(); itW++) {
      
      s32 coverage = 0;
      SSRContig* exonB = *itW;
      if( !distantFrom(v, exonB, SSRContigList::MAXSIZEINTRON) && distantFrom(v, exonB, SSRContigList::MINSIZEINTRON)) {
	ostringstream oss;
	oss << v->seqName() << "@" << v->start()<<"@"<<v->end() << "@"<<exonB->start()<<"@" << exonB->end() << "@" << v->strand();
	string keyJunction = oss.str();
	map<string,s32>::iterator itKJ = KJunctions.find(keyJunction);
	if(itKJ != KJunctions.end()) {
	  coverage = itKJ->second;
	  Edge pairEdge =  make_pair(v->getID(), exonB->getID());
	  edges->push_back( pairEdge);
	  map<Edge,s32>::iterator itW = weight.find(pairEdge);
	  if(itW == weight.end())
	    weight.insert(make_pair(pairEdge,coverage));
	  else{
	    if(itW->second != coverage){
	      cout << " junctions2edges else junction found  "<< itW->first.first << "-"<<itW->first.second  << "  "<< itW->second<< " coverage " << coverage << endl;
	      cout << "Error itW->second != coverage " << endl;
	      exit(1);
	    }
	  }
	}
      }
    }
  }
  return edges;
}

// to print the covtig list
ostream& operator<<(ostream& ostr, const SSRContigList& d) {
  TSSRList::const_iterator itList;
  for( itList=d._contigs.begin(); itList != d._contigs.end(); itList++ ) {
    if(itList != d._contigs.begin()) { ostr << endl; }
    SSRContig* ctg = *itList;
    ostr << *ctg;
  }
  return ostr;
}
