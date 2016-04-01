/*******************************************************************************
+
+  ReadFile.cpp
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
 *******************************************************************************/

#include <malloc.h>
#include "ReadFile.h"
#include <time.h>

s32 ReadFile::VERBOSE;

using namespace std;

// to check what file format is given as input, accepts fasta, customized fasta and fastq
void ReadFile::_checkFormat() {
  string s1,s2,s3,s4;
  char ch[1000];
  _fstrm.open(_filename, ios_base::in | ios_base::binary);
  _fstrm>>s1;
  _fstrm.getline(ch, 1000);
  _fstrm>>s2;
  _fstrm.getline(ch, 1000);
  string fasta_form = "AaCcGgTtNn>";
  
  if('>' == s1[0]){
    _fstrm>>s3;
    if(fasta_form.find(s3[0]) != string::npos) {
      _file_format = 1;} //genome or read fasta file
    else _file_format = 3;} //read fasta file with number of occurences
  else if('@' == s1[0]) {
    _fstrm>>s3;
    _fstrm.getline(ch, 1000);
    _fstrm>>s4;
    _fstrm.getline(ch, 1000);
    _file_format = 2; //read fastq file
    if(s2.size() != s4.size()) {
      cerr<<"fatal error: fq format, sequence length not equal to quality length\n";
      exit(1);
    }
  }
  else {
    cerr<<"fatal error: unrecognizable format of reads file.\n";
    exit(1);
  }
  _fstrm.seekg(0);
}

// to fill the dictionary with the DNA sequence (DnaDictionary splits the sequence into words)
s32 ReadFile::loadAndCount(DnaDictionary& dict) {
  string sSeq;
  s32 nbSeq = 0;
  s32 nbOCcur = 0;
  
  s32 cpt=0;
  
  if (_file_format==3){
    while( next(sSeq, nbOCcur) != 0 ) {
      dict.countWords(sSeq, nbOCcur);
      nbSeq=nbSeq+nbOCcur;
    }
  }
  else{
    while( next(sSeq,nbOCcur) != 0 ) {
      dict.countWords(sSeq, 1);
      nbSeq++;
      cpt++;
      if(cpt==10000) { if(ReadFile::VERBOSE) cerr << "\r  -> " << nbSeq << " processed" << flush; cpt=0; }
    }
  }
  if(ReadFile::VERBOSE) cerr << "\r  -> " << nbSeq << " processed" << endl;
  return nbSeq;
}

// to get the next sequence in file
int ReadFile::next(string& seq, s32& nbRead) {
  char ch[1000];
  char c;
  string name, qual, tmp_str;
  
  _fstrm>>c;
  if(_fstrm.eof()) return 0;
  _fstrm>>name;
  _fstrm>>seq;
  if (_file_format==3) _fstrm>>nbRead; //read fasta file with nb of occurences
  if (_file_format==1) { //genome or read fasta file
    _fstrm>>c;
    while (('>' != c) && (!_fstrm.eof())){
      seq += c;
      _fstrm>>tmp_str;
      seq += tmp_str;
      tmp_str.erase();
      _fstrm>>c;
    }
  }
  if(_file_format==2) {//read fastq file
    _fstrm>>ch;
    _fstrm.getline(ch, 1000);
    _fstrm>>qual;
  }
  return 1;
}
