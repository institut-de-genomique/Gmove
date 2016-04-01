/*******************************************************************************
+
+  GeneModel.h
+
+  
+  
+ 
 *******************************************************************************/

#ifndef GENEMODEL_H
#define GENEMODEL_H

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <map>
#include "SSRContig.h"
#include "LocalType.h"
#include "DnaDictionary.h"
#include <iterator>
#include <ctime>

typedef list<SSRContig*> TSSRList;
typedef map<string, s32, less<string>, allocator<pair<const string, s32> > > DefaultMap;

extern bool PRINTTIME;

class GeneModel {
  s32 _strand;
  string _seqname;
  s32 _indMod;
  string _filename;
  list<s32> _path;
  TSSRList _exons;
  TSSRList* _whole_exons;
  
  pair<s32,s32> _cds;
  bool _start_found;
  bool _stop_found;
  s32 _cds_size;
  s32 _window;
  s32 _3pXL;
  s32 _5pXL;
  s32 _ref_score;
  map<string, s32> _protPhaseCoord;
  s32 _model_size;
  
 public:
  static bool UNORIENTED;
  
  /* Constructors and Destructors*/
  GeneModel(list<s32>& path, TSSRList* listExons, string fname, s32 iter, s32 window, map<string, list<s32> >& probExons) { //FIXME why we give all of exon from _vertices
    _filename = fname;
    _path = path;  
    _whole_exons = listExons;
    _indMod = iter;
    _start_found = 0;
    _stop_found = 0;
    _cds_size = 0;
    _model_size = 0 ;
    _cds = make_pair(0, 0);
    _window = window;
    _3pXL = 0;
    _5pXL = 0;
    _ref_score = 0;
    _protPhaseCoord = DefaultMap();
    for(list<s32>::iterator itP = _path.begin() ; itP != _path.end() ; itP++){

    	_exons.push_back(this->exon(*itP));
    	_model_size += this->exon(*itP)->size();
    	if(this->exon(*itP)->getPosNegID()) {
    		ostringstream oss;
    		oss << this->exon(*itP)->seqName() << "@" << this->exon(*itP)->start() << "@" << this->exon(*itP)->end();
    		string key = oss.str();
    		cout << "key GeneModel " << key << endl;
    		map<string, list<s32> >::iterator itPE;
    		itPE = probExons.find(key);
    		if(itPE == probExons.end()) {
    			list<s32> listita;
    			listita.push_back(iter);
    			probExons.insert( make_pair(key,listita) );
    		}
    		else (itPE->second).push_back(iter);
    	}
    }
    _exons.sort(SSRContig::startposSort);
    if(_exons.size() > 0) {
    	_strand = (*_exons.begin())->strand();
    	_seqname = (*_exons.begin())->seqName();
    }
  }
  
  ~GeneModel() { }
  
  /* Accessors */
  TSSRList& getExons() { return _exons; }
  string getFilename() const { return _filename; }
  s32 getID() { return _indMod; }
  string getSeqname() const { return _seqname; }
  list<s32>& getPath() { return _path; }
  pair<s32,s32> getCDS() const { return _cds; } 
  bool containCDS() const { return (_cds.first && _cds.second); }
  s32 getModelSize()const {return _model_size;}
  bool getStartFound()const { return _start_found;}
  bool getStopFound() const {return _stop_found;}
  
  /* Methods */
  bool findORF(const map<string, s32>& protPhaseCoord = DefaultMap());
  s32 selectORF(const map<string, s32>& protPhaseCoord = DefaultMap());
  s32 getPhaseScore(pair<s32, s32> orf, s32 strand, s32 mrna_len, const map<string, s32>& modelPhaseCoord = DefaultMap());
  pair<s32,s32> longestORF(vector<s32> lStart, vector<s32> lStop,map<string,bool>&alreadyLookPos);
  pair<s32,s32> bestORF(vector<s32> lStart, vector<s32> lStop, s32 strand, s32 mrna_len, s32& score,map<string,bool>& alreadyLookPos, const map<string,s32>& modelPhaseCoord = DefaultMap());
  void mapORF();
  SSRContig* exon(s32 id);
  void formatGTF (ofstream &forf, s32 nbp);
  void printAnnot(ofstream &forf, s32 nbp, bool format);
  void formatGFF (ofstream &forf, s32 nbp);
  void keepModel(map<string,GeneModel>& mapGeneModel);

  bool overlapOrf(pair<s32,s32>);
};

char complem (char base);
void selectModel( map<string,GeneModel> mapGeneModel);
void selectModel(map<s32, list<GeneModel> >&);
void selectModel(list<GeneModel>& );

#endif
