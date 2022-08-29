/*******************************************************************************
+
+  SSRContigLists.h
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
 *******************************************************************************/

#ifndef JM_SSR_CONTIG_LISTS_H
#define JM_SSR_CONTIG_LISTS_H

#include "LocalType.h"
#include "SSRContig.h"
#include "SSRContigList.h"
#include "NetEx.h"
#include "GffRecordList.h"

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <string>
#include <map>
#include <list>

using namespace std;
using namespace boost;


typedef map<string, SSRContigList*> TcontigLists;

class SSRContigLists {
  
  friend ostream& operator<<(ostream&, const SSRContigLists&);
  
 protected:
  TcontigLists _contigs;
  // map of Junctions from the junctionFile
  map<string, s32> _kwJunctions;
  map<s32,list<string> > _transcrit;
  
 public:
  /* Constructors and Destructors*/
  SSRContigLists() { /*_nbContigs=0;*/ }
  SSRContigLists(GffRecordList listRecord, map<string, string>& sequences, map<string, s32>& badintrons);

  ~SSRContigLists() {
    for(TcontigLists::iterator itCtg =  _contigs.begin(); itCtg != _contigs.end();++itCtg)
      delete itCtg->second;
  }
  
  s32 addJunction(string, s32, s32, s32, s32, s32, s8, map<string, s32>&);
  list<NetEx*> buildGraph();
  void cleanJunctions();
};


#endif
