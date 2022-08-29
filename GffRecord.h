/*******************************************************************************
+
+  GffRecord.h
+
+  Copyright (c) 2017 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury & Marion Dubarry, 
+           jmaury@genoscope.cns.fr
+
*******************************************************************************/

#ifndef GFFRECORD_H
#define GFFRECORD_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <list>
#include <algorithm>
#include <string>
#include <vector>

#include "LocalType.h"
using namespace std;

class GffRecord {
  friend ostream& operator<<(ostream&, const GffRecord&);

 protected:
  string _seqname;
  string _method;
  string _type;
  s32 _start;
  s32 _end;
  f4 _score;
  string _strand;
  string _phase;
  string _attribute;
  s32 _color;
  s8 _datatype;
  
 public:
  
  /* Constructors and Destructors*/
  GffRecord(string fileName, s8 datatype);
  ~GffRecord() {
  }
  
  void prepareAttribute(string& attribute);
  void split(string s, string delimiter);
  
  /* Accessors */
  string getSeqName() const { return _seqname; }
  string getMethod() const { return _method; }
  string getType() const { return _type; }
  s32 getStart() const { return _start; }
  s32 getEnd() const { return _end; }
  f4 getScore() const { return _score; }
  string getStrand() const { return _strand; }
  string getPhase() const { return _phase; }
  string getAttribute() const { return _attribute; }
  s32 getColor() const {return _color;}
  s8 getDatatype() const { return _datatype; }
  
  void setType(string s){_type = s;}
  void setStart(s32 s){_start = s;}
  void setEnd(s32 e ){_end = e;}
  void setColor(s32 c){_color = c;}
};

#endif
