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
#include "DnaDictionary.h"
//#include "PETags.h"

#include <boost/shared_ptr.hpp>
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <seqan/basic.h>
#include <seqan/gff_io.h>
#include <seqan/stream.h>

using namespace std;
using namespace boost;
using namespace seqan;

typedef map<string, SSRContigList*> TcontigLists;


class SSRContigLists {

	friend ostream& operator<<(ostream&, const SSRContigLists&);

protected:
	TcontigLists _contigs;
	s32 _nbContigs;
	map<string, s32> _kwJunctions; // map of Junctions from the junctionFile

public:
	/* Constructors and Destructors*/
	SSRContigLists() { _nbContigs=0; }
//	SSRContigLists(char*, char*, map<string, string>&);
	SSRContigLists(list<pair< GffRecord, string> > , map<string, string>& ) ;
		~SSRContigLists() {
		for(TcontigLists::iterator itCtg =  _contigs.begin(); itCtg != _contigs.end();++itCtg){
			delete itCtg->second;
		}
	}

	list<NetEx*> testJunctions(DnaDictionary&);
	void extendCovtigs(DnaDictionary&);
	void cleanJunctions();
};


#endif
