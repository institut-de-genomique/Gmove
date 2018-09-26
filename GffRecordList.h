/*******************************************************************************
+
+  GffRecordList.h
+
+  Copyright (c) 2017 Genoscope, CEA, CNS, Evry, France
+  Author : Marion Dubarry, jmaury@genoscope.cns.fr
+
*******************************************************************************/

#ifndef GFFRECORD_H_
#define GFFRECORD_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>     /* exit, EXIT_FAILURE */
#include <sstream>
#include <list>
#include <map>

#include "GffRecord.h"
#include "LocalType.h"

using namespace std;

typedef list<GffRecord*> GffRecordL;

class GffRecordList {

	protected:

	GffRecordL _records;

 public:

  /* Constructors and Destructors*/
	//Lsta_scaffold10 Gmove   mRNA    100764  101283  .       +       .       ID=mRNA.Lsta_scaffold10.1.1;Name=mRNA.Lsta_scaffold10.1.1;start=1;stop=1;cds_size=165;model_size=520
	GffRecordList(){};
	GffRecordList(char * fileName,map<string,string>fastaDB, string typeData);

  ~GffRecordList() {
/*	  for(GffRecordL::iterator itGffRecord = _records.begin(); itGffRecord != _records.end(); ++itGffRecord){
		  delete *itGffRecord;
	 }
*/
  }
  void copyGffRecordList(GffRecordList gffRecord);

  void printGffRecord();
  map<string,GffRecordL> extractMono();
  void checkFormatGff(char* line);
  void loadGff(map<string, map<s32,s32> >& cds);
  void loadAnnotation( map<string, map<s32,s32> >& cds, bool useCds);
  bool empty(){return _records.empty();};
  GffRecordL* getRecords(){return &_records;}
  void insertMap(map<string,GffRecordL>& mapRecord, GffRecordL listMono);
  void fusionMonoExons(GffRecordL& listMono);
  void cleanMono();
  void intron();
};

void checkUniqId(map<string,bool>& mapAttribute, string previousAttribute);
void cdsPhase(s32& phase, list< pair< s32,s32> > pos,  string name, map<string,map<s32,s32> >& cds);
bool sortMonoGff (GffRecord* record1 ,  GffRecord* record2);
bool sortRecordGff (GffRecord*  record1 , GffRecord* record2);




#endif
