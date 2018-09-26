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
//#include <seqan/basic.h>
//#include <seqan/gff_io.h>
//#include <seqan/stream.h>

using namespace std;
using namespace boost;
//using namespace seqan;

typedef map<string, SSRContigList*> TcontigLists; // map<scaffold, listExons>


class SSRContigLists {

	friend ostream& operator<<(ostream&, const SSRContigLists&);

protected:
	TcontigLists _contigs;
//	s32 _nbContigs;
	map<string, s32> _kwJunctions; // map of Junctions from the junctionFile
	map<s32,list<string> > _transcrit;

public:
	/* Constructors and Destructors*/
	SSRContigLists() { /*_nbContigs=0;*/ }
	SSRContigLists(GffRecordList listRecord, map<string, string>& sequences) ;
		~SSRContigLists() {
		for(TcontigLists::iterator itCtg =  _contigs.begin(); itCtg != _contigs.end();++itCtg){
			delete itCtg->second;
		}
	}

	list<NetEx*> buildGraph();
//	void extendCovtigs(DnaDictionary&);
	void cleanJunctions();
};


#endif
