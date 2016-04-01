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

// to test junctions of the dictionary and the junctions file given as input
NetEx* SSRContigList::testJunctions(DnaDictionary& dict, map<string, s32>& KJunctions, ofstream& out) {
  TSSRList::iterator itPrec, itNext;
  map<string, s32> _seenJunctions;
  map<string, bool> _seenExons;
  s32 nbCovtigs = 0, posnegCount = 0;
  string seqname;

  if(SSRContigList::VERBOSE) cerr << "validation of junctions :" << endl;
  for( itPrec = _contigs.begin(); itPrec != _contigs.end(); itPrec++) {
    if(SSRContigList::VERBOSE) cerr << "\r  -> " << ++nbCovtigs << " of " << _contigs.size() << flush;
    SSRContig* covtigA = *itPrec;
    if(seqname.empty()) seqname = covtigA->seqName();
    if(!covtigA->isPopulate()) { this->populateExons(covtigA, _seenExons); }
    s32 i = 1;
    for( itNext = itPrec, itNext++; (itNext != _contigs.end() && i <= SSRContigList::NBNEIGHBOUR); i++, itNext++ ) {
      SSRContig* covtigB = *itNext;
 //     cout << "\n covtigB " << *covtigB << endl;
      if (!overlap(covtigA, covtigB) && sameStrand(covtigA,covtigB)){
   // 	  cout << *covtigA << " covtigB "<< *covtigB << endl;
    	  _testJunctionsBetweenCovtigs(covtigA, covtigB, dict, KJunctions, _seenJunctions, _seenExons, out);
      }
      else i--;
    }
  }
  cleanExons(posnegCount);
  TSSRList* vertices = exons2vertices();
  list<Edge>* edges = junctions2edges(vertices);
  NetEx* ExNtk = new NetEx(vertices, edges, seqname);
  return ExNtk;
}

// to test junctions between 2 covtigs
void SSRContigList::_testJunctionsBetweenCovtigs(SSRContig* covtigA, SSRContig* covtigB, DnaDictionary& dict, map<string, s32>& KJunctions, map<string, s32>& seenJunctions,  map<string, bool>& _seenExons, ofstream& out) {
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
  //  cout <<"exonA->seqName() " << exonA->seqName()<<" " <<exonA->start()<< " "<< exonA->end()<< " strand "<< exonA->strand() ;
    idA++; idB=-1;
    if(covtigA->isKnownStrand() && covtigA->strand() != exonA->strand()) { continue;}
    
    // For each exonB
    for(TSSRList::iterator itExB = LexonsB->begin(); itExB != LexonsB->end(); itExB++) {
      SSRContig* exonB = *itExB;
      idB++;
      if(covtigB->isKnownStrand() && covtigB->strand() != exonB->strand()) { continue;}
   //   cout <<" exonB->seqName() " << exonB->seqName()<<" " <<exonB->start()<< " "<< exonB->end()<< " strand "<< exonB->strand() << endl;
      _testJunctionsBetweenExons(exonA, exonB, dict, KJunctions, seenJunctions, out);
      if(idB >= SSRContigList::REDUCEEXONLIST) { itExB = LexonsB->end(); itExB--; }
    }
    if(idA >= SSRContigList::REDUCEEXONLIST) { itExA = LexonsA->end(); itExA--; }
  }

}

// to test junctions between 2 exons
s32 SSRContigList::_testJunctionsBetweenExons(SSRContig* exonA, SSRContig* exonB, DnaDictionary& dict, map<string, s32>& KJunctions, map<string, s32>& seenJunctions, ofstream& out) {
  ostringstream oss;
  s32 nb = -1;
  SSRContig *ctg1 = exonA->getMaster(), *ctg2 = exonB->getMaster();

  if (exonA->strand() != exonB->strand()) return nb; // TO PREVENT EXONS WITH CONTRADICTORY JUNCTIONS
   
  oss << exonA->seqName() << "@" << exonA->start()<<"@"<<exonA->end() << "@" << exonB->start()<<"@"<<exonB->end() << "@" << exonA->strand();

  string key = oss.str();
  string sA, sB, siteA, siteB;
  // 1. We look if junction between those 2 exons is already validated ?
  map<string, s32>::iterator itJunctions = seenJunctions.find(key);
  if (itJunctions != seenJunctions.end()){
	  nb = itJunctions->second; //We save the nb of time this junction has been validated
  }
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

		out << ctg1->seqName() <<" "<<  exonA->end() << " " << exonB->start() << " " << exonA->strand()
			<< " siteA= " << siteA << " siteB= " << siteB
			<< " coverageA= " << exonA->coverage() << " coverageB= " << exonB->coverage()
			<< " covtigA= " << ctg1->start() << "-" << ctg1->end() << " covtigB= "
			<< ctg2->start() << "-" << ctg2->end() << endl;
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

		out << ctg1->seqName() <<" "<<  exonA->end() << " " << exonB->start() << " " << exonA->strand()
			<< " TotW= " << nb << " WordList= " << listeJ << " siteA= " << siteA << " siteB= " << siteB
			<< " coverageA= " << exonA->coverage() << " coverageB= " << exonB->coverage()
			<< " covtigA= " << ctg1->start() << "-" << ctg1->end() << " covtigB= "
			<< ctg2->start() << "-" << ctg2->end() << endl;
      }
      seenJunctions.insert(make_pair(key, nb));
    }

  }

  if(nb > 0) {
//	  cout << " jonction is validated between "<< exonA->start() << " " << exonA->end() << " exonB " << exonB->start() << " " << exonB->end() << endl;
 //   if(exonA->strand() == SSRContig::FORWARD) (exonA->getLinkedWith())->push_back(exonB); //XXX Apres avoir vérifié la jonction, on met le sens du lien
 //   else (exonB->getLinkedWith())->push_back(exonA);
//	  cout << " jonction validates exonA start " << exonA->start() << " exonBstart " << exonB->start() << endl;
	  (exonA->getLinkedWith())->push_back(exonB); // FIXME I change sens of graph
    if(nb >= SSRContigList::MINCOVWORD2ORIENTATE) {
      ctg1->setStrand(exonA->strand());
      ctg2->setStrand(exonB->strand());
    }
    ctg1->set_isAsingleexon(false);
    ctg2->set_isAsingleexon(false);

    // MODIF JM - 04/08/2015 - modif apportee en r28, retour en r27
    // Le changement de statut des exons ne se fait que dans un sens, on ne doit pas changer directement 
    // le statut des exons ayant d�ja un statut !!!
    //    if(exonA->strand() == SSRContig::FORWARD) {
    exonB->setValstate(1); // exonB could be 0 or 1 but we set it at 1
    if (exonA->getValstate() == 0) exonA->setValstate(2);
    if (exonA->getValstate() == 1) {	
      exonA->setValstate(3);
      ctg1->set_isdblval(true);
    }
      //    }
      //    else {
      //exonA->setValstate(1);
      //if (exonB->getValstate() == 0) exonB->setValstate(2);
      //if (exonB->getValstate() == 1) {	
      //exonB->setValstate(3);
      //ctg2->set_isdblval(true);
      //}
    //}
    // FIN MODIF JM - 04/08/2015 - modif apportee en r28, retour en r27
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
 // if(SSRContigList::VERBOSE) cerr << "Number of vertice(s): " << vertices->size() << endl;
  return vertices;
}

// to add junctions between vertice in the list of graph edges
list<Edge>* SSRContigList::junctions2edges(TSSRList* vertices) {
  list<Edge>* edges = new list<Edge>;
  
  for(TSSRList::iterator itV = vertices->begin(); itV != vertices->end(); itV++) {
    SSRContig* v = *itV;
    if(v->isAFutureVertice()) {
      TSSRList* linkedWith = v->getLinkedWith();
      for(TSSRList::iterator itEx = linkedWith->begin(); itEx != linkedWith->end(); itEx++) {
    	  SSRContig* exonB = *itEx;
    	  if(exonB->isAFutureVertice()) edges->push_back( make_pair(v->getID(), exonB->getID()) );
      } 
    }
  }
 // if(SSRContigList::VERBOSE) cerr << "Number of edge(s): " << edges->size() << endl;
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
    //  if(right.find(exon->start()) != right.end()) { cout << " we find right" << endl;break; } //XXX This is a test : if the end exon already exist, we do not validate the other !
      if(!witness3) { witness3 = n+1; n++; exon->isAFutureVertice(true); }
      exon->setEnd(ctg->end());
      exon->setID(witness3-1);
      keep = true;
      right.insert(make_pair(exon->start(), true));
    }

    // exon is first exon, move exon start position to covtig start position and keep only one ID for all first exons
    if(!ctg->isdbval() && exon->getValstate() == 2) {
   //   if(left.find(exon->end()) != left.end()) { break; }
      if(!witness5) { witness5 = n+1; n++; exon->isAFutureVertice(true);  }
      exon->setStart(ctg->start());
      exon->setID(witness5-1);
      keep = true; 
      left.insert(make_pair(exon->end(), true));
    }

    // exon seems to be a single-exon model and is longer than the minimal exon size
    if(ctg->isAsingleexon() && exon->getValstate() == 0 && ctg->size() > SSRContig::MINSIZEEXON) {
      if(!witness35) {
		witness35 = n+1;
		exon->setID(witness35-1);
		n++;
		exon->setStart(ctg->start());
		exon->setEnd(ctg->end());
		exon->isAFutureVertice(true);
		keep = true;
      }
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
	} else {
	  enlargeLeft = 0;
	}
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
	} else { 
	  enlargeRight = 0;
	}
      }
    }
    if(joinPrevious || joinNext) {
      if(joinPrevious) {
	previous->setEnd(current->end());
	current = previous;
      }
      if(joinNext) {
	next->setStart(current->start());
      }
      itPrec = _contigs.erase(itPrec);
    } else { itPrec++; }
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
