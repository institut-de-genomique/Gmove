/*******************************************************************************
+
+  DnaDictionary.h
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
*******************************************************************************/

#ifndef JM_DNA_DICTIONARY_H
#define JM_DNA_DICTIONARY_H

#include "LocalType.h"
#include <time.h>
#include <iostream>
#include <string>
#include <ext/hash_map>

using namespace std;
namespace std { using namespace __gnu_cxx; }

typedef u64 TDnaWord;
struct eqwd {bool operator()(TDnaWord s1, TDnaWord s2) const{return s1 == s2;}};
typedef hash_map<TDnaWord, u16, hash<TDnaWord>, eqwd> Tdictionary;


class DnaDictionary {
  friend ostream& operator<<(ostream&, DnaDictionary&);

  static const s32 REPLACE_N_BY_A = 0;

 //protected:
 public:
  Tdictionary _dictionary;


 private:
  s32 _wordlen;
  s64 _nbWords;
  s64 _nbDiffWords;
  s32 _badchar;
  s32 _discardedW;
  bool _testlowcomplexity;

  TDnaWord _fwdMask;
  char* bintoascii;

  static void _initWord(TDnaWord& w) { w=0; }

  void _insert(TDnaWord&, s32 multiplicator=1);

  s32 _packedFwdWord(TDnaWord& w, char c) {
    w &= _fwdMask;
    w <<= 2;
    s32 old_badchar = _badchar;
    w |= char2Int(c);
    if(!REPLACE_N_BY_A && old_badchar != _badchar) { return 0; }
    return 1;
  }
  s32 _packedRevWord(TDnaWord& w, char c) {
    w >>= 2;
    s32 old_badchar = _badchar;
    w |= ((u64) (3 ^ char2Int(c)) << (2 * _wordlen - 2));
    if(!REPLACE_N_BY_A && old_badchar != _badchar) { return 0; }
    return 1;
  }
  s32 char2Int(char ch) {
    switch  (tolower (ch))
      {
      case  'a' : return  0;
      case  'c' : return  1;
      case  'g' : return  2;
      case  't' : return  3;
      default: _badchar++; return 0;
      };
    return  0;
  }

 public:
  /* Constructors and Destructors*/
  DnaDictionary(u32 word_len) { 
    _dictionary = Tdictionary();
    _fwdMask = ((u64) 1 << (2 * word_len - 2)) - 1;
    bintoascii = "ACGT";
    _wordlen = word_len;
    _nbWords = 0;
    _nbDiffWords = 0;
    _badchar = 0;
    _discardedW = 0;
    _testlowcomplexity = true;
  }
  ~DnaDictionary() {}

  /* Accessors */
  s32 getWordSize() const { return _wordlen; }
  s32 getNbWords() const { return _nbWords; }
  s32 getNbDiffWords() const { return _nbDiffWords; }
  s32 getNbBadChars() const { return _badchar; }
  Tdictionary* getDictionary() { return &_dictionary; }
  s32 getNbDiscardedWords() {return _discardedW;}
  bool empty() const { return !_nbWords; }
  bool testLowComplexity(bool val) { _testlowcomplexity = val; return _testlowcomplexity; }
  bool testLowComplexity() const { return _testlowcomplexity; }

  /* Methods */
  void bin2String(TDnaWord& w, string& s) {
    s.erase();
    s.resize(_wordlen);

    for (s32 i = 0; i < _wordlen; i++) {
      s[_wordlen-i-1] = bintoascii[w & 0x3];
      w >>= 2;
    }
  }
  bool beginWord(const string&, s32&, TDnaWord&, TDnaWord&);
  void countWords (const string&, s32 multiplicator=1);
  bool existWord(const TDnaWord&);
  bool existWord(const string&);
  s32 nbOccWord(const TDnaWord&);
  s32 nbOccWord(const string&);
  bool isLowComplexity(string&);
  void cleanUM();
  
};

#endif
