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
#include <iterator>
#include <ctime>

typedef list<SSRContig*> TSSRList;
typedef map<string, s32, less<string>, allocator<pair<const string, s32> > > DefaultMap;

extern bool PRINTTIME;

class GeneModel {
  s32 _strand;
  string _seqname;
  s32 _indMod;
  s32 _num_cc;
  s32 _genetic_code;
 string _filename; //FIXME Why GeneModel contain filename ?
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
  bool _toDelete;
  s32 _cluster;

  f8 _score_prot;
  f8 _score_lencds;
  f8 _score_input;
  f8 _score_mrna;
  
 public:
  static bool UNORIENTED;
  
  /* Constructors and Destructors*/
  GeneModel(list<s32>& path, TSSRList* listExons, string fname, s32 iter, s32 window, s32 genetic_code, map<string, list<s32> >& probExons) { //FIXME why we give all of exon from _vertices
    _filename = fname;
    _cluster = 0;
    _path = path;  
    _whole_exons = listExons;
    _num_cc = iter;
    _genetic_code = genetic_code;
    _indMod = 1;
    _start_found = 0;
    _stop_found = 0;
    _cds_size = 0;
    _model_size = 0 ;
    _cds = make_pair(0, 0);
    _window = window;
    _3pXL = 0;
    _5pXL = 0;
    _ref_score = 0;
    _toDelete = false;
    _protPhaseCoord = DefaultMap();
    for(list<s32>::iterator itP = _path.begin() ; itP != _path.end() ; itP++){
      _exons.push_back(this->exon(*itP));
      _model_size += this->exon(*itP)->size();
      //	cout <<"exon " << this->exon(*itP)->start() << " "<< this->exon(*itP)->end()<<endl;
      if(this->exon(*itP)->getPosNegID()) {
	ostringstream oss;
	oss << this->exon(*itP)->seqName() << "@" << this->exon(*itP)->start() << "@" << this->exon(*itP)->end();
	string key = oss.str();
    	//	cout << "key GeneModel " << key << endl;
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
    _score_prot = 0.0;
    _score_lencds = 0.0;
    _score_input = 0.0;
    _score_mrna = 0.0;
  }
 
  ~GeneModel() { }
  
  friend ostream& operator<<(ostream&, const GeneModel&);
  /* Accessors */
  bool getToDelete() const {return _toDelete;}
  void setToDelete(bool b){_toDelete = b;}
  TSSRList& getExons() { return _exons; }
  s32 nbExons() { return _exons.size(); }
  string getFilename() const { return _filename; }
  s32 getNumCC() const { return _num_cc; }
  s32 getGeneticCode() const { return _genetic_code; }
  s32 getID() const { return _indMod; }
  void setID(s32 id){_indMod = id;}
  string getSeqname() const { return _seqname; }
  list<s32>& getPath() { return _path; }
  pair<s32,s32> getCDS() const { return _cds; } 
  bool containCDS() const { return (_cds.first && _cds.second); }
  s32 getModelSize()const {return _model_size;}
  bool getStartFound()const { return _start_found;}
  bool getStopFound() const {return _stop_found;}
  s32 get3pXL()const {return _3pXL;}
  s32 get5pXL() const {return _5pXL;}
  s32 getCdsSize() const {return _cds_size;}
  s32 getCluster()const {return _cluster;}
  s32 getStrand()const { return _strand;}
  void setCluster(s32 num){_cluster = num;}
  s32 getMrnaLen() {
    s32 len = 0;
    for(TSSRList::iterator it = _exons.begin() ; it != _exons.end(); it++) {
      len += ((*it)->end() - (*it)->start() + 1); 
    }
    return len;
  }
  
  f8 score_prot() { return _score_prot; }
  f8 score_lencds() { return _score_lencds; }
  f8 score_input() { return _score_input; }
  f8 score_mrna() { return _score_mrna; }
  void score_prot(f8 arg) { _score_prot = arg; }
  void score_lencds(f8 arg) { _score_lencds = arg; }
  void score_input(f8 arg) { _score_input = arg; }
  void score_mrna(f8 arg) { _score_mrna = arg; }
  f8 score() { 
    f8 min_val = 0.1;
    f8 s = (_score_input < min_val) ? 0.1 : _score_input;
    return ((f8)_cds_size * s);
  }


  /* Methods */
  pair<s32,s32> getTransBound();
  s32 getModelSizeonGeno();
  bool findORF(const map<string, s32>& protPhaseCoord = DefaultMap());
  s32 selectORF(const map<string, s32>& protPhaseCoord = DefaultMap());
  s32 getPhaseScore(pair<s32, s32> orf, s32 strand, s32 mrna_len, const map<string, s32>& modelPhaseCoord = DefaultMap());
  pair<s32,s32> longestORF(vector<s32> lStart, vector<s32> lStop);
  pair<s32,s32> bestORF(vector<s32> lStart, vector<s32> lStop, s32 strand, s32 mrna_len, s32& score, const map<string,s32>& modelPhaseCoord = DefaultMap());
  void mapORF();
  SSRContig* exon(s32 id);
  void formatGTF (ofstream &forf);
  void printAnnot(ofstream &forf, bool format);
  void formatGFF (ofstream &forf);

  void keepModel(map<string,GeneModel>& mapGeneModel);
  bool cdsIsMono(void);
  list<pair<s32,s32> > getExonsCds();

  bool overlapOrf(pair<s32,s32>);
  bool overlapOrf(GeneModel model2);
  bool overlapModel(GeneModel model2);



};

char complem (char base);
void selectModel(map<s32, list<GeneModel> >&);

#endif
