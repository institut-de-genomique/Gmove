/*******************************************************************************
+
+  SSRContig.cpp
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
*******************************************************************************/

#include "SSRContig.h"

using namespace std;

s32 SSRContig::EXTEND;
s32 SSRContig::MINSIZEEXON;
s32 SSRContig::MINCOVSPLICESITES;
s32 SSRContig::MINNBWORD;


// to build a list of potential exons from the enlarged sequence of the covtig
void SSRContig::populateExons() {
  string sequence = getEnlargeSeq(SSRContig::EXTEND);
  s32 from = max(0, _startpos - 1 - SSRContig::EXTEND);
  if(this->strand() == SSRContig::UNKNOWN){
	  _populate(sequence, from, splices3for, splices5for, SSRContig::FORWARD);
	  _populate(sequence, from, splices5rev, splices3rev, SSRContig::REVERSE);
  }
  else
	  _populate(sequence, from, splices3for, splices5for, this->strand());
}

// to build a list of potential exons from an oriented sequence
void SSRContig::_populate(string& sequence, s32 from, vector<string>& splicesLeft, vector<string>& splicesRight, Tstrand strand) {
vector<string>::iterator it;
  vector<s32> sitesL, sitesR;
  /* found 5' splice sites */
  for ( it = splicesLeft.begin() ; it < splicesLeft.end(); it++ ) {
    string::size_type loc = sequence.find(*it, 0);
    while( loc != string::npos && (s32)loc <= 2*SSRContig::EXTEND) {
      sitesL.push_back((s32)loc+3);
      loc = sequence.find(*it, loc+1);
    }
  }
  /* found 3' splice sites */
  for ( it = splicesRight.begin() ; it < splicesRight.end(); it++ ) {
    string::size_type loc = sequence.find(*it, sequence.size()-2*SSRContig::EXTEND);
    while( loc != string::npos ) {
      sitesR.push_back((s32)loc);
      loc = sequence.find(*it, loc+1);
    }
  }

  /* derive all potential exons from 5' and 3' splice sites */
  vector<s32>::iterator itLSites, itRSites;
  for ( itLSites = sitesL.begin() ; itLSites < sitesL.end() ; itLSites++ ) {
    for ( itRSites = sitesR.begin() ; itRSites < sitesR.end() ; itRSites++ ) {
      s32 start = (*itLSites)+from, stop = (*itRSites)+from;
      if((stop - start + 1) >= SSRContig::MINSIZEEXON && 
	 start <= _endpos-SSRContig::MINSIZEEXON &&
	 stop >= _startpos+SSRContig::MINSIZEEXON) {
    	  SSRContig* exon = new SSRContig(_seqname, start, stop, _coverage, _wholeseq,strand);//FIXME memory leak
    	  exon->setMaster(this);
    	  _exons.push_back(exon);
      }
    }
  }

  /* covtig adds itself to its exons list */
  if((_endpos - _startpos + 1) >= SSRContig::MINSIZEEXON &&
		  _startpos <= _endpos-SSRContig::MINSIZEEXON &&
	 _endpos >= _startpos+SSRContig::MINSIZEEXON) {
   SSRContig* exonSelf = new SSRContig(_seqname, _startpos, _endpos, _coverage, _wholeseq,strand); //FIXME memory leak
  exonSelf->setMaster(this);
  exonSelf->setTag(_tag);
  _exons.push_back(exonSelf);
  }
  /* sorts all exons according to distance with covtig */
  _exons.sort(SortByDistWithMaster());
}

// to get the enlarged sequence of the covtig
string SSRContig::getEnlargeSeq(s32 extend) const { 
  s32 start = max(0, _startpos-1-extend);
  s32 size = min((s32)_wholeseq->length(), _endpos-start+extend);
  return _wholeseq->substr(start, size);
}

// to get the enlarged sequence of the covtig on its 5' side
string SSRContig::getLeftSeq(s32 fromStart) const { 
  s32 size = min((s32)_wholeseq->length(), fromStart);
  return _wholeseq->substr(_startpos-1, size);
}

// to get the enlarged sequence of the covtig on its 3' side
string SSRContig::getRightSeq(s32 fromEnd) const { 
  s32 start = max(0, _endpos-fromEnd);
  s32 size = min((s32)_wholeseq->length(), fromEnd);
  return _wholeseq->substr(start, size);
}

// to know if 2 covtigs are distant from mindist
bool distantFrom(SSRContig* ctg1, SSRContig* ctg2, s32 mindist) {
  if(ctg1->start() < ctg2->start()) { return ctg1->end() + mindist < ctg2->start(); }
  else { return ctg2->end() + mindist < ctg1->start(); }
}

// to know if 2 covtigs overlap
bool overlap(SSRContig* ctg1, SSRContig* ctg2) {
  return !distantFrom(ctg1, ctg2, 0);
}

void leftRight(SSRContig* ctg, s32& left,s32& right){
	if(ctg->strand() == SSRContig::REVERSE){
		 left = ctg->end();
		 right = ctg->start();
	 }
	 else{
		 left = ctg->start();
		 right = ctg->end();
	 }
}
bool sameEnd(SSRContig* ctg1, SSRContig* ctg2){
	return ctg1->end() == ctg2->end();
}

bool sameStart(SSRContig* ctg1, SSRContig* ctg2){
	return ctg1->start() == ctg2->start();
}


bool sameStrand(SSRContig* ctg1,SSRContig* ctg2){
	return (ctg1->strand() == ctg2->strand()) || (ctg1->strand() == SSRContig::UNKNOWN) || (ctg2->strand()== SSRContig::UNKNOWN);
}
// to order 2 covtigs based on their starting position
bool SSRContig::startposSort(const SSRContig* ctg1, const SSRContig* ctg2) { return (ctg1->start() < ctg2->start()); }

// to order 2 covtigs based on their ending position
bool SSRContig::endposSort(const SSRContig* ctg1, const SSRContig* ctg2) { return (ctg1->end() < ctg2->end()); }


bool SSRContig::startposSortAndLength(const SSRContig* ctg1, const SSRContig* ctg2) { return (ctg1->start() < ctg2->start() && ctg1->size() > ctg2->size()); }

bool SSRContig::endposSortAndLength(const SSRContig* ctg1, const SSRContig* ctg2) { return (ctg1->end() < ctg2->end() && ctg1->size() > ctg2->size()); }

// to calculate the number of times each junction between 2 covtigs is seen in the dictionary
s32 nbJunctions(SSRContig* exon1, SSRContig* exon2, DnaDictionary& dict, s32 min_overlap, string& listeJ) {
  s32 nbTot = 0, first = 1, nbWords = 0;//, nbOffset = 0, nbThreshold = 0;
  ostringstream oss;
  for(s32 i = dict.getWordSize()-min_overlap ; i >= min_overlap ; i--) {
    string word, wordR;
    word = exon1->getRightSeq(i);
    wordR = exon2->getLeftSeq(dict.getWordSize()-i);
    word.append(wordR);
    s32 nb = dict.nbOccWord(word);
    nbTot += nb;
    if(nb > 0) nbWords++;
    if(!first) oss << ",";
    first = 0;
    oss << nb;
  }
  // If the minimal number of different words detected is not higher than the treshold, return 0
  if(nbWords < SSRContig::MINNBWORD) nbTot = 0;
  listeJ = oss.str();
  return nbTot;
}

// to print a covtig's information
ostream& operator<<(ostream& ostr, const SSRContig& d) {
  return ostr<< d.getID()<<" " << d.seqName() << " " << d.start() << " " << d.end() << " " << d.strand();
}

// to get start-end in a string form
string SSRContig::getName() {
  string name_exon = convert(start(), end());
  return name_exon;
  
}

//** exon labelling with boundaries **//
string convert(s32 start, s32 end) {
  std::stringstream ss;
  ss << start << "-" << end;
  return ss.str();
}



