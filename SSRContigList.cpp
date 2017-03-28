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
s32 SSRContigList::REDUCEEXONLIST;
s32 SSRContigList::VERBOSE;


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
//	cout << " mapCovtig "<<endl;
	time_t  before = time(NULL);
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
	 time_t  after = time(NULL);
	// cout << " time SSRContigList::mapCovtigs " << difftime(after,before)<<endl;
	return mCov;
}

map<pair<string, s32>,list<string> > SSRContigList::mapJunctions(map<string,s32> junction){
	time_t  before = time(NULL);
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
		else{
			itMJunction->second.push_back(itJunction->first);
		}
	}

	 time_t  after = time(NULL);
//	 cout << "time SSRContigList::mapJunctions " << difftime(after,before)<<endl;
	return mJunction;
}

NetEx* SSRContigList::testJunctions(DnaDictionary& dict, map<string, s32>& KJunctions) {
	TSSRList::iterator itPrec, itNext;
	map<string, s32> _seenJunctions;
	map<string, bool> _seenExons;
	s32 nbCovtigs = 0, posnegCount = 0;
	string seqname;

	map<s32,list<SSRContig*> > mCov;
	map<pair<string,s32>,list<string> > mJunction;
	mCov = mapCovtigs(_contigs);
	mJunction = mapJunctions(KJunctions);
	if(SSRContigList::VERBOSE) cerr << "validation of junctions :" << endl;

    for( itPrec = _contigs.begin(); itPrec != _contigs.end(); itPrec++) {
		if(SSRContigList::VERBOSE) cerr << "\r  -> " << ++nbCovtigs << " of " << _contigs.size() << " name " << (*itPrec)->seqName() << flush;
		SSRContig* covtigA = *itPrec;
		if(seqname.empty()) seqname = covtigA->seqName();
				 if(!covtigA->isPopulate()) { this->populateExons(covtigA, _seenExons); }
		map<pair<string,s32>,list<string> >::iterator itMJunction;
		pair<string,s32> tmpPair = make_pair(covtigA->seqName(),covtigA->end());
		itMJunction = mJunction.find(tmpPair);
		if(itMJunction == mJunction.end()){
			continue;
		}
		for(list<string>::iterator ititMJunction = itMJunction->second.begin() ; ititMJunction != itMJunction->second.end(); ++ititMJunction ){
			vector<string> splitKey = split(*ititMJunction,'@');
			map<s32,list<SSRContig*> >::iterator itMCov;
			itMCov = mCov.find(atoi(splitKey[3].c_str()));
			if(itMCov == mCov.end()){
				continue;
			}
			for(list<SSRContig*>::iterator ititMCov = itMCov->second.begin(); ititMCov != itMCov->second.end(); ++ititMCov){
				 SSRContig* covtigB = *ititMCov;
				 if (!overlap(covtigA, covtigB) && sameStrand(covtigA,covtigB)){
						  _testJunctionsBetweenCovtigs(covtigA, covtigB, dict, KJunctions, _seenJunctions, _seenExons);
				 }
			}
		}
    }
  cleanExons(posnegCount);
  TSSRList* vertices = exons2vertices();
  map<Edge,s32> weight;
  list<Edge>* edges = junctions2edges(vertices,KJunctions, weight);
  NetEx* ExNtk = new NetEx(vertices, edges, seqname,weight);
  return ExNtk;
}


// to test junctions between 2 covtigs
void SSRContigList::_testJunctionsBetweenCovtigs(SSRContig* covtigA, SSRContig* covtigB, DnaDictionary& dict, map<string, s32>& KJunctions, map<string, s32>& seenJunctions,  map<string, bool>& _seenExons) {
  TSSRList* LexonsA = covtigA->getExons();
  if(!covtigB->isPopulate()) { this->populateExons(covtigB, _seenExons); }
  /*XXX Why we populate covtigB here and covtigA in SSRContigList::testJunctions
   * We need to put one more argument in SSRContigList::_testJunctionsBetweenCovtigs : map<string, s32>& seenJunctions   *
   */
  TSSRList* LexonsB = covtigB->getExons();
  s32 idA=-1, idB=-1;
  // For each exonA
  for(TSSRList::iterator itExA = LexonsA->begin(); itExA != LexonsA->end(); itExA++) {
    SSRContig* exonA = *itExA;
    idA++; idB=-1;
    if(covtigA->isKnownStrand() && covtigA->strand() != exonA->strand()) { continue;}
    
    // For each exonB
    for(TSSRList::iterator itExB = LexonsB->begin(); itExB != LexonsB->end(); itExB++) {
      SSRContig* exonB = *itExB;
      idB++;
      if(covtigB->isKnownStrand() && covtigB->strand() != exonB->strand()) { continue;}
      _testJunctionsBetweenExons(exonA, exonB, dict, KJunctions, seenJunctions);
      if(idB >= SSRContigList::REDUCEEXONLIST) { itExB = LexonsB->end(); itExB--; }
    }
    if(idA >= SSRContigList::REDUCEEXONLIST) { itExA = LexonsA->end(); itExA--; }
  }

}

// to test junctions between 2 exons
s32 SSRContigList::_testJunctionsBetweenExons(SSRContig* exonA, SSRContig* exonB, DnaDictionary& dict, map<string, s32>& KJunctions, map<string, s32>& seenJunctions) {
  ostringstream oss;
  s32 nb = -1;
  SSRContig *ctg1 = exonA->getMaster(), *ctg2 = exonB->getMaster();
  if (exonA->strand() != exonB->strand()) return nb; // TO PREVENT EXONS WITH CONTRADICTORY JUNCTIONS
   
  oss << exonA->seqName() << "@" << exonA->start()<<"@"<<exonA->end() << "@" << exonB->start()<<"@"<<exonB->end() << "@" << exonA->strand();

  string key = oss.str();
  string sA, sB, siteA, siteB;
  // 1. We look if junction between those 2 exons is already validated ?
  map<string, s32>::iterator itJunctions = seenJunctions.find(key);
  if (itJunctions != seenJunctions.end())
	  nb = itJunctions->second; //We save the nb of time this junction has been validated

  // present in known junctions (input file)
  if(nb < 0) { // The junction was never validated
	  if( !distantFrom(exonA, exonB, SSRContigList::MAXSIZEINTRON) && distantFrom(exonA, exonB, SSRContigList::MINSIZEINTRON)) {
		  itJunctions = KJunctions.find(key); // we look if the junction exist from the junctionFile

		  if(itJunctions != KJunctions.end()) {
			  nb = itJunctions->second;
			  if(nb >= 1) {
				sA = exonA->getEnlargeSeq(2);
				sB = exonB->getEnlargeSeq(2);
				siteA = sA.substr(sA.size()-2, 2);
				siteB = sB.substr(0, 2);
				transform(siteA.begin(), siteA.end(), siteA.begin(), (int(*)(int)) toupper);
				transform(siteB.begin(), siteB.end(), siteB.begin(), (int(*)(int)) toupper);
			  }
		  }
	  }
  }
  // need to check in the word dictionary
  if(nb < 0) {
    nb = 0;
    if(exonA->strand() == exonB->strand() && !distantFrom(exonA, exonB, SSRContigList::MAXSIZEINTRON) && distantFrom(exonA, exonB, SSRContigList::MINSIZEINTRON)) {
    	string listeJ;
      nb = nbJunctions(exonA, exonB, dict, SSRContigList::MINOVERLAPJUNCTION, listeJ);
      if(nb >= 1) {
		string sA = exonA->getEnlargeSeq(2);
		string sB = exonB->getEnlargeSeq(2);
		string siteA = sA.substr(sA.size()-2, 2);
		string siteB = sB.substr(0, 2);
		transform(siteA.begin(), siteA.end(), siteA.begin(), (int(*)(int)) toupper);
		transform(siteB.begin(), siteB.end(), siteB.begin(), (int(*)(int)) toupper);
      }
      seenJunctions.insert(make_pair(key, nb));
    }

  }

  if(nb > 0) {
	  (exonA->getLinkedWith())->push_back(exonB); // FIXME I change sens of graph
    if(nb >= SSRContigList::MINCOVWORD2ORIENTATE) {
      ctg1->setStrand(exonA->strand());
      ctg2->setStrand(exonB->strand());
    }
    ctg1->set_isAsingleexon(false);
    ctg2->set_isAsingleexon(false);

     exonB->setValstate(1); // exonB could be 0 or 1 but we set it at 1
    if (exonA->getValstate() == 0) exonA->setValstate(2);
    if (exonA->getValstate() == 1) {	
      exonA->setValstate(3);
      ctg1->set_isdblval(true);
    }
  }
  return nb;
}

// to check if potential exons are validated as being vertice in the graph
TSSRList* SSRContigList::exons2vertices() {
  TSSRList* vertices = new TSSRList();
  for(TSSRList::iterator it = _contigs.begin(); it != _contigs.end(); it++) {
    SSRContig* ctg = *it;
    TSSRList* exons = ctg->getExons();

    for(TSSRList::iterator itEx = exons->begin(); itEx != exons->end(); itEx++) {
      SSRContig* exon = *itEx;
      if(exon->isAFutureVertice()) { vertices->push_back( exon ); }
    }
  }
  return vertices;
}

// to add junctions between vertices in the list of graph edges
list<Edge>* SSRContigList::junctions2edges(TSSRList* vertices, map<string,s32> KJunctions,  map<Edge,s32>& weight) {
  list<Edge>* edges = new list<Edge>;

  for(TSSRList::iterator itV = vertices->begin(); itV != vertices->end(); itV++) {
    SSRContig* v = *itV;
    if(v->isAFutureVertice()) {
      TSSRList* linkedWith = v->getLinkedWith();
      for(TSSRList::iterator itEx = linkedWith->begin(); itEx != linkedWith->end(); itEx++) {
    	  SSRContig* exonB = *itEx;
    	  s32 coverage = 0;
    	  if(exonB->isAFutureVertice()){
    		  ostringstream oss;
    		  oss << v->seqName() << "@" << v->start()<<"@"<<v->end() << "@"<<exonB->start()<<"@" << exonB->end() << "@" << v->strand();
    		  string keyJunction = oss.str();
    		  map<string,s32>::iterator itKJ = KJunctions.find(keyJunction);
    		  if(itKJ == KJunctions.end())
    			  cerr << " error junction not found : SSRContigList::junctions2edges " <<keyJunction << " exonA " << *v << " exonB " << *exonB << endl; //FIXME probleme des fois on ne trouve pas la jonction alors que les exons sont liés. Verifié comment est fait  la liste linkedWith
    		  else
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

// to clean exons that are not validated or have discrepancies (present in junctions on strand + and -)
void SSRContigList::cleanExons(s32& posnegCount) {
  s32 n = 0, tmp_n = 0;
  map<s32, bool> rightSites, leftSites;
  for(TSSRList::iterator it = _contigs.begin(); it != _contigs.end(); it++){

      tmp_n = n;
      n = cleanExon(*it, n, rightSites, leftSites);
      if(n-tmp_n==2) { //XXX What does that mean ?
    	  SSRContig* tmp_it = *it;
    	  tmp_it->setPosNegID(++posnegCount);
      }
  }
}

// to clean exons that are not validated
s32 SSRContigList::cleanExon(SSRContig* ctg, s32 n, map<s32, bool>& right, map<s32, bool>& left){
  s32 witness3 = 0, witness5 = 0, witness35 = 0;
  TSSRList* exons = ctg->getExons();
  for(TSSRList::iterator it = exons->begin(); it != exons->end(); ) {

    SSRContig* exon = *it;
    bool keep = false;
    // exon is double validated - 3prim and 5prim
    if(exon->getValstate() == 3) { keep = true; exon->setID(n++); exon->isAFutureVertice(true); }

    // exon is last exon, move exon end position to covtig end position and keep only one ID for all last exons
    if(!ctg->isdbval() && exon->getValstate() == 1) {
      if(!witness3) { witness3 = n+1; n++; exon->isAFutureVertice(true); }
      exon->setEnd(ctg->end());
      exon->setID(witness3-1);
      keep = true;
      right.insert(make_pair(exon->start(), true));
    }

    // exon is first exon, move exon start position to covtig start position and keep only one ID for all first exons
    if(!ctg->isdbval() && exon->getValstate() == 2) {
      if(!witness5) { witness5 = n+1; n++; exon->isAFutureVertice(true);  }
      exon->setStart(ctg->start());
      exon->setID(witness5-1);
      keep = true; 
      left.insert(make_pair(exon->end(), true));
    }

    // exon seems to be a single-exon model and is longer than the minimal exon size
    if(ctg->isAsingleexon() && exon->getValstate() == 0 && ctg->size() > SSRContig::MINSIZEEXON) {
		witness35 = n+1;
		exon->setID(witness35-1);
		n++;
		exon->setStart(ctg->start());
		exon->setEnd(ctg->end());
		exon->isAFutureVertice(true);
		keep = true;
    }
    // have to delete the current exon
    if(!keep) {
      it = exons->erase(it);
      delete exon;
    } else it++;
  }

  return n;
}

// to extend the covtigs using the word dictionary
void SSRContigList::extendCovtigs(DnaDictionary& dict) {
	TSSRList::iterator itPrec, itNext;
	SSRContig *previous=0, *current=0, *next=0;
	for(itPrec = _contigs.begin(); itPrec != _contigs.end(); ) {
		current = *itPrec;
		bool enlargeLeft=1, enlargeRight=1;
		bool joinPrevious=0, joinNext=0;
		while(enlargeLeft || enlargeRight) {
			itNext = itPrec;
			itNext++;
			next = (itNext == _contigs.end()) ? 0 : *itNext;
			s32 offset = (current->size() < dict.getWordSize()) ? dict.getWordSize() - (s32)current->size() : 1;
			string covtig = current->getEnlargeSeq(offset);
			if((s32)covtig.length() < dict.getWordSize()) { enlargeLeft=enlargeRight=0; }
			if(current->start() == 1) { enlargeLeft=0; }
			if(current->end() == current->getSeqSize()) { enlargeRight=0; }
			if(enlargeLeft) {
				string word = covtig.substr(0, dict.getWordSize());
				s32 nb = dict.nbOccWord(word);
				if(nb>0) {
					current->setStart(current->start()-1);
					if(previous && current->start() < previous->end()+1) {
						current->setStart(previous->start());
						enlargeLeft = 0;
						joinPrevious = 1;
					}
				} else
				enlargeLeft = 0;
			}
			if(enlargeRight) {
				string word = covtig.substr((current->size()+2*offset)-dict.getWordSize()-1, dict.getWordSize());
				s32 nb = dict.nbOccWord(word);
				if(nb>0) {
					current->setEnd(current->end()+1);
					if(next && current->end() > next->start()-1) {
						current->setEnd(next->end());
						enlargeRight = 0;
						joinNext=1;
					}
				}
				else
					enlargeRight = 0;
			}
		}
		if(joinPrevious || joinNext) {
			if(joinPrevious) {
				previous->setEnd(current->end());
				current = previous;
			}
			if(joinNext)
				next->setStart(current->start());
			itPrec = _contigs.erase(itPrec);
		}
		else { itPrec++; }
			previous = current;
	}
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
