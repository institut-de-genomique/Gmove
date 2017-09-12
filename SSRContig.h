/*******************************************************************************
+
+  SSRContig.h
+
+  Copyright (c) 2002 Genoscope, CEA, CNS, Evry, France
+  Author : Jean-Marc Aury, jmaury@genoscope.cns.fr
+ 
 *******************************************************************************/

#ifndef JM_SSR_CONTIG_H
#define JM_SSR_CONTIG_H

#include "LocalType.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <map>

using namespace std;

typedef s8 Tstrand;
typedef u8 Valstate;
   
class SSRContig {
  friend ostream& operator<<(ostream&, const SSRContig&);
  friend bool overlap(SSRContig*, SSRContig*);
  friend void leftRight(SSRContig*, s32&,s32&);
  friend bool sameEnd(SSRContig* ctg1, SSRContig* ctg2);
  friend bool sameStart(SSRContig* ctg1, SSRContig* ctg2);
  friend bool sameStrand(SSRContig*, SSRContig*);// return True if same strand or if for one the strand is unknown, not enough info
  friend bool distantFrom(SSRContig*, SSRContig*, s32);
 // friend s32 nbJunctions(SSRContig*, SSRContig*, s32, string&);
  friend bool operator < (const SSRContig& e1, const SSRContig& e2) { return e1._startpos < e2._startpos; }

 protected:
  string _seqname;
  s32 _startpos;
  s32 _endpos;
  f8 _coverage;
  Tstrand _strand;
  s32 _id;
  std::map<s32,bool> _idTranscrit;
  Valstate _valstate; //diff validation states
  bool _dblval;
  bool _monoex;
  bool _isAFutureVertice;
  s32 _posnegID;
  
  string* _wholeseq;
  list<SSRContig*> _exons;
  list<SSRContig*> _linkedWith;
  SSRContig* _master;
  string _splice5;
  string _splice3;
 // string _tag; // define if the exon is a start, interne, end or mono exons

  vector<string> splices5for, splices3for, splices5rev, splices3rev;
  
  void _populate(string&, s32, vector<string>&, vector<string>&, Tstrand);
  
 public:
  
  const static Tstrand FORWARD =  1;
  const static Tstrand REVERSE = -1;
  const static Tstrand UNKNOWN =  0;
  
  const static Valstate noV = 0;
  const static Valstate begV = 1;
  const static Valstate endV = 2;
  const static Valstate bothV = 3;

  static s32 EXTEND;
  static s32 MINSIZEEXON;
  static s32 MINCOVSPLICESITES;
  static s32 MINCOVWORD2ORIENTATE;
  static s32 MINNBWORD;
  
  static bool startposSort(const SSRContig*, const SSRContig*);
  static bool endposSort(const SSRContig*, const SSRContig*);
  static bool startposSortAndLength(const SSRContig*, const SSRContig*);
  static bool endposSortAndLength(const SSRContig*, const SSRContig*);

  /* Constructors and Destructors*/  
  SSRContig(string& seqname, s32 st, s32 en, f8 coverage, string* seq, Tstrand strand) {
    _seqname = seqname;
    _startpos = st;
    _endpos = en;
    _coverage = coverage;
    _wholeseq = seq;
    _strand = strand;
    _master = NULL;
    _id = -1;
    _valstate = 0;
    _dblval = false;
    _monoex = true;
    _isAFutureVertice = false;
    _posnegID = 0;
 //   _idTranscrit.insert(make_pair(idTranscrit,true));
    splices5for.push_back("gt");
    splices5for.push_back("gc");
    splices3for.push_back("ag");
    splices5rev.push_back("ac");
    splices5rev.push_back("gc");
    splices3rev.push_back("ct");

  }
  
  
  ~SSRContig() {  
    for( list<SSRContig*>::iterator itExon = _exons.begin(); itExon!=_exons.end();++itExon){
    	delete *itExon;
    }
  }

  /* Accessors */
  s32 start() const { return _startpos; }
  s32 end() const { return _endpos; }
  f8 coverage() const { return _coverage; }
  s32 strand() const { return _strand; }
//  string tag() const { return _tag;}
  bool isKnownStrand() const { return (_strand != UNKNOWN); }
  string seqName() const { return _seqname; }
  string getSeq() const { return getEnlargeSeq(0); }
  s32 getSeqSize() const { return _wholeseq->length(); }
  list<SSRContig*>* getExons() { return &_exons; }
  list<SSRContig*>* getLinkedWith() { return &_linkedWith; }
  bool isPopulate() { return !(_exons.empty()); }
  SSRContig* getMaster() const { return _master; }
  string* getWholeSeq() const { return _wholeseq; }
  s32 getValstate() const { return _valstate;}
  bool isdbval() const { return _dblval; }
  bool isAsingleexon() const { return _monoex; }
  bool isAFutureVertice() const { return _isAFutureVertice; }
  bool isAFutureVertice(bool val) { _isAFutureVertice = val; return _isAFutureVertice; }
  s32 getID() const { return _id; }
  s32 size() const { return _endpos-_startpos+1; }
  s32 getPosNegID() const { return _posnegID; }
  
//  s32 distWithMaster() const { return (abs(_startpos-_master->start()) + abs(_endpos-_master->end())); }

  void setStrand(Tstrand strand) { _strand = strand; }
  void setMaster(SSRContig* master) { _master = master; }
  void setID(s32 ident){_id = ident;}
  void setIdTranscrit(s32 id){_idTranscrit.insert(make_pair(id,true));}
  void setValstate (s32 val){_valstate = val;}
  void set_isdblval(bool val) {_dblval = val;}
  void set_isAsingleexon(bool val) {_monoex = val;}
  void setExons (list<SSRContig*>* listEx) {_exons = *listEx;}
  void setStart(s32 start){_startpos = start;}
  void setEnd(s32 end){_endpos = end;}
  void setSplice5(string s) { _splice5 = s; }
  void setSplice3(string s) { _splice3 = s; }
  void setPosNegID(s32 posnegID) { _posnegID = posnegID; }
//  void setTag(string t){_tag = t;}
  
  /* Methods */
  string getEnlargeSeq(s32) const;
  string getLeftSeq(s32) const;
  string getRightSeq(s32) const;
  map<s32,bool> getIdTranscrit()const{return _idTranscrit;};
  void populateExons();
  void clearExons() {_exons.clear();}
  string getName();
};

string convert(s32 start, s32 end);

// Definition du foncteur servant a trier nos objets selon le nombre
/*struct SortByDistWithMaster
{ 
   bool operator() (const SSRContig* e1, const SSRContig* e2) const 
   { 
      return e1->distWithMaster() < e2->distWithMaster();
   } 
}; 
*/
#endif
