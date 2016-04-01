/*******************************************************************************
+
+  DnaDictionary.cpp
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
 *******************************************************************************/
#include <malloc.h>
#include "DnaDictionary.h"

using namespace std;

// to insert a new word in the dictionary
void DnaDictionary::_insert(TDnaWord& w, s32 multiplicator) {
  Tdictionary::iterator it = _dictionary.find(w);
  if (it == _dictionary.end()) {
    it = _dictionary.insert(make_pair(w,0)).first;
    _nbDiffWords++;
  }
  it->second+=multiplicator;
  _nbWords++;
}

// to shift in position each time a "N" is found (begin a new word) to only add determined words to the dictionary
bool DnaDictionary::beginWord(const string& s, s32& pos, TDnaWord& fwd, TDnaWord& rev) {
  fwd = 0;
  rev = 0;
  s32 i = 0, N = s.length();

  if ((N - pos) < _wordlen) { return false; }
  for(i = 0; i < _wordlen - 1;  i++) {
    if(_packedFwdWord(fwd, s[pos])) {
      _packedRevWord(rev, s[pos]);
      pos++;
    }
    else {
      pos++;
      return beginWord(s, pos, fwd, rev);
    }
  }
  return true;
}

// to extract all words present in a string and insert them into the dictionary
void DnaDictionary::countWords(const string& s, s32 multiplicator) {
  TDnaWord fwd = 0, rev = 0;
  s32 i = 0, N = s.length();
  if(!beginWord(s, i, fwd, rev)) {return;}
  
  Tdictionary::iterator it;
  while (i < N) {
    if(!_packedFwdWord(fwd, s[i])) { 
      i++;
      if(!beginWord(s, i, fwd, rev)) {
	return;}
      continue;
    }
    _packedRevWord(rev, s[i]);
    
    string word = s.substr(i-_wordlen+1, _wordlen);
    if(!isLowComplexity(word))
      if (fwd < rev) { _insert(fwd, multiplicator); }
      else { _insert(rev, multiplicator);}    
    else _discardedW++;
    i++;
  }
  return;
}

// to test existence of a word (TDnaWord form)
bool DnaDictionary::existWord(const TDnaWord& w) {
  Tdictionary::iterator it = _dictionary.find(w);
  if (it == _dictionary.end()) return FALSE;
  return TRUE;
}

// to test existence of a word (string form)
bool DnaDictionary::existWord(const string& s) {
  if((s32)s.length() != _wordlen) return FALSE;
  TDnaWord fwd=0, rev=0;
  
  s32 oldBadchar = _badchar;
  for  (s32 i = 0;  i < _wordlen;  i++) {
    _packedFwdWord(fwd, s[i]);
    _packedRevWord(rev, s[i]);
  }
  if(oldBadchar != _badchar){ 
    _badchar = oldBadchar; //_badchar represents bad chars present in the dictionary, should not change when testing words
    return FALSE;
  }
  //return (fwd < rev) ? existWord(fwd) : existWord(rev);
  Tdictionary::iterator it = _dictionary.find(fwd);
  if (it != _dictionary.end()) return TRUE;
  it = _dictionary.find(rev);
  if (it != _dictionary.end()) return TRUE;
  return FALSE;
}

// to find the number of occurences of a word (TDnaWord form)
s32 DnaDictionary::nbOccWord(const TDnaWord& w) {
  Tdictionary::iterator it = _dictionary.find(w);
  if (it == _dictionary.end()) return 0;
  return it->second;
}

// to find the number of occurences of a word (string form)
s32 DnaDictionary::nbOccWord(const string& s) {
  if((s32)s.length() != _wordlen) return 0;
  TDnaWord fwd=0, rev=0;
  for  (s32 i = 0;  i < _wordlen;  i++) {
    _packedFwdWord(fwd, s[i]);
    _packedRevWord(rev, s[i]);
  }
  Tdictionary::iterator it = _dictionary.find(fwd);
  if (it != _dictionary.end()) return it->second; 
  it = _dictionary.find(rev);
  if (it != _dictionary.end()) return it->second; 
  return 0;
}

ostream& operator<<(ostream& ostr, DnaDictionary& d) {
  Tdictionary::iterator it;
  for( it=(d._dictionary).begin(); it != (d._dictionary).end(); it++ ) {
    if(it != (d._dictionary).begin()) { ostr << endl; }
    string word;
    TDnaWord w = it->first;
    d.bin2String(w, word);
    ostr << word << " " << it->second;
  }
  return ostr;
}

// to test a word for low complexity
 bool DnaDictionary::isLowComplexity(string& mot) { //MTF algo
  if(!_testlowcomplexity) return false;
  u8 i,j;
  char c;
  char dico[4];
  dico[0] = 'a';
  dico[1] = 't';
  dico[2] = 'g';
  dico[3] = 'c';
  s32 score = 0;

  for (i=0; i < mot.length(); i++){
    c = tolower(mot[i]);	    
    for (j=0; j<4; j++){
      if (dico[j]==c) break;	    	
    }
    score=score+j;
    while (j>0){
      dico[j]=dico[j-1];
      j--;
    }
    dico[0]=c;	   
  }
  if (score <= ((s32)mot.length() / 2)) return true;
  else return false;
}

// to  clean unique words
void DnaDictionary::cleanUM (){
  Tdictionary::iterator it, itPrec;
  for( it = _dictionary.begin(); it != _dictionary.end(); ) {    
    if(it->second == 1) { itPrec=it; it++; _dictionary.erase(itPrec); }
    else it++;
  }
}



