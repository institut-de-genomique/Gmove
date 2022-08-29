/*******************************************************************************
+
+  GffRecordList.h
+
+  Copyright (c) 2017 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury & Marion Dubarry, 
+           jmaury@genoscope.cns.fr
+
*******************************************************************************/

#ifndef GFFRECORD_H_
#define GFFRECORD_H_

#include <iostream>
#include <fstream>
#include <stdlib.h>
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
  GffRecordList(){};
  GffRecordList(char*, map<string,string>, s8);
  
  ~GffRecordList() {
  }

  void clean() {
    for(GffRecordL::iterator itGffRecord = _records.begin(); itGffRecord != _records.end(); ++itGffRecord)
      delete *itGffRecord;
  }
  void copyGffRecordList(GffRecordList);
  
  void printGffRecord();
  map<string,GffRecordL> extractMono();
  void checkFormatGff(char*);
  void loadGff(map<string, map<s32,s32> >&);
  void loadAnnotation( map<string, map<s32,s32> >&, bool);
  bool empty() { return _records.empty(); };
  GffRecordL* getRecords() { return &_records; }
  void insertMap(map<string,GffRecordL>&, GffRecordL);
  void fusionMonoExons(GffRecordL&);
  void cleanMono();
  void intron();
};

void checkUniqId(map<string,bool>&, string);
void cdsPhase(s32&, list< pair< s32,s32> >, string, map<string,map<s32,s32> >&);
bool sortMonoGff (GffRecord*, GffRecord*);
bool sortRecordGff (GffRecord*, GffRecord*);

#endif
