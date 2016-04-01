/*******************************************************************************
+
+  ReadFile.h
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
*******************************************************************************/

#ifndef JM_READ_FILE_H
#define JM_READ_FILE_H

#include "LocalType.h"
#include "DnaDictionary.h"
#include <iostream>
#include <fstream>
#include <string>
#include <time.h>
#include <map>

using namespace std;

class ReadFile {
  
 private:
  u32 _file_format;
  char* _filename;
  fstream _fstrm;
  
  void _checkFormat();
  
 public:
  static s32 VERBOSE;
  
  /* Constructors and Destructors*/
  ReadFile(char* filename) { 
    _filename = filename;
    _file_format = 0;
    _checkFormat();
  }
  ~ReadFile() {}
  
  /* Accessors */
   
  /* Methods */
  int next(string&, s32&);
  s32 loadAndCount(DnaDictionary&);
};

#endif
