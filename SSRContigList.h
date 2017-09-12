/*******************************************************************************
+
+  SSRContigList.h
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
*******************************************************************************/

#ifndef JM_SSR_CONTIG_LIST_H
#define JM_SSR_CONTIG_LIST_H

#include "LocalType.h"
#include "SSRContig.h"
#include "NetEx.h"

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <list>

using namespace std;

typedef list<SSRContig*> TSSRList;
typedef pair<s32,s32> Edge;

class SSRContigList {
  friend ostream& operator<<(ostream&, const SSRContigList&);
  
 protected:
  TSSRList _contigs;


 public:
  static s32 MINOVERLAPJUNCTION;
  static s32 NBNEIGHBOUR;
  static s32 MINSIZEINTRON;
  static s32 MAXSIZEINTRON;
  static s32 MINCOVWORD2ORIENTATE;
  static s32 REDUCEEXONLIST;
  static s32 VERBOSE;
  
  /* Constructors and Destructors*/
  SSRContigList() {}
  ~SSRContigList() {
	  for(TSSRList::iterator itCtg = _contigs.begin();itCtg!=_contigs.end();++itCtg){
		  delete *itCtg;
	  }
  }
  
  void push_back(SSRContig* ctg) { return _contigs.push_back(ctg); }
  s32 size() { return _contigs.size(); }
  SSRContig* back() { return _contigs.back();}
  void print(){//TODO delete, just for tests
	  for(TSSRList::iterator itCtg = _contigs.begin();itCtg!=_contigs.end();++itCtg){
	 		cout << **itCtg <<endl;
	 	  }
  }

  void populateExons(SSRContig*, map<string, bool>&);
  
  map<s32,list<SSRContig*> > mapCovtigs(TSSRList);
  map<pair<string, s32>,list<string> > mapJunctions(map<string,s32> junction);

  NetEx* buildGraph(map<string, s32>&);

//  void extendCovtigs();
  TSSRList* exons2vertices();
  list<Edge>* junctions2edges(TSSRList*, map<string,s32>, map<Edge,s32>&);

  void startposSort() { _contigs.sort(SSRContig::startposSort); }
  void endposSort() { _contigs.sort(SSRContig::endposSort); }



};
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) ;
std::vector<std::string> split(const std::string &s, char delim);

#endif
