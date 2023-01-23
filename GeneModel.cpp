/*******************************************************************************
+
+  GeneModel.cpp
+
+  
+  
+ 
 *******************************************************************************/

#include "GeneModel.h"

bool GeneModel::UNORIENTED;
bool PRINTTIME = false;
// to print the annotation in the output file
void GeneModel::printAnnot(ofstream &forf, bool format) {
  if (format == true) this->formatGFF(forf);
  else this->formatGTF(forf);
}

pair<s32,s32> GeneModel::getTransBound(){
  s32 debTrans = (*_exons.begin())->start();
  s32 finTrans = (*(--_exons.end()))->end();

  s32 debCDS = _cds.second, finCDS = _cds.first;
  if (_cds.first < _cds.second) {
    debCDS = _cds.first;
    finCDS = _cds.second;
  }

  if (_strand >= 0) {
    if(_start_found == 0)
      debTrans = debCDS;
    else if (_stop_found == 0)
      finTrans = finCDS;

    debTrans -= _5pXL;
    finTrans += _3pXL;
  }
  else if (_strand < 0) {
    if(_start_found == 0)
      finTrans = finCDS;
    else if (_stop_found == 0)
      debTrans = debCDS;

    finTrans += _5pXL;
    debTrans -= _3pXL;
  }
  return make_pair(debTrans, finTrans);
}

s32 GeneModel::getModelSizeonGeno() {
  pair<s32,s32> transBound = this->getTransBound();
  return transBound.second - transBound.first + 1;
}

// GFF3 file outputing: format fitting our in-house needs
void GeneModel::formatGFF (ofstream &forf) {
  string strdSt = (_strand < 0) ? "-" : "+";
  
  s32 debCDS = _cds.second, finCDS = _cds.first;
  if (_cds.first < _cds.second) {
    debCDS = _cds.first;
    finCDS = _cds.second;
  } 
  /*	s32 debTrans = (*_exons.begin())->start();
	s32 finTrans = (*(--_exons.end()))->end();
	
	if (_strand >= 0) {
	if(_start_found == 0)
	debTrans = debCDS;
	else if (_stop_found == 0)
	finTrans = finCDS;
	
	debTrans -= _5pXL;
	finTrans += _3pXL;
	}
	else if (_strand < 0) {
	if(_start_found == 0)
	finTrans = finCDS;
	else if (_stop_found == 0)
	debTrans = debCDS;
	
	finTrans += _5pXL;
	debTrans -= _3pXL;
	}
  */
  pair<s32,s32> transBound = this->getTransBound();
  s32 debTrans = transBound.first;
  s32 finTrans = transBound.second;
  string argCDS;
  ostringstream oss;
  f8 score = this->score();

  oss << ";start="<<_start_found<< ";stop=" << _stop_found << ";cds_size=" << _cds_size<<";model_size="<<_model_size+_5pXL+_3pXL<<";exons="<<_exons.size();
  argCDS = oss.str();	  
  forf << _seqname << "\tGmove\tmRNA\t" << debTrans << "\t" << finTrans << "\t" << score << "\t" << strdSt << "\t.\t" << "ID=mRNA."
       << _seqname << "." << _num_cc << "." << _indMod
       << ";Name=mRNA." << _seqname << "." << _num_cc << "." << _indMod
       << argCDS << endl;
  
  for (TSSRList::iterator it = _exons.begin(); it != _exons.end(); it++){
    s32 debExon = (*it)->start();
    s32 finExon = (*it)->end()  ;
    if (_strand >=0) {
      if (it == _exons.begin()) {debExon -= _5pXL;}
      if (it == --_exons.end()) {finExon += _3pXL;}
    }
    else if(_strand < 0) {
      if (it == --_exons.end()) {finExon += _5pXL;}
      if (it == _exons.begin()) {debExon -= _3pXL;}
    }
    if( finCDS < debExon || debCDS > finExon ) {
      forf << _seqname << "\tGmove\tUTR\t" << debExon << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA."
	   << _seqname << "." << _num_cc << "." << _indMod << endl;
    }
    else {
      s32 end = (finCDS < finExon) ? finCDS : finExon;
      s32 begin = (debCDS < debExon) ? debExon : debCDS;
      if(begin != debExon){
	if( (_strand > 0 && _start_found) || (_strand < 0 && _stop_found ) ){
	  forf << _seqname << "\tGmove\tUTR\t" << debExon << "\t" << begin-1 << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA."
	       << _seqname << "." << _num_cc << "." << _indMod << endl;
	}
      }
      forf << _seqname << "\tGmove\tCDS\t" << begin << "\t" << end << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA."
	   << _seqname << "." << _num_cc << "." << _indMod << endl;
      if(end != finExon){
	if( (_strand > 0 && _stop_found) || (_strand < 0 && _start_found) ){
	  forf << _seqname << "\tGmove\tUTR\t" << end+1 << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA."
	       << _seqname << "." << _num_cc << "." << _indMod << endl;
	}
      }
    }
  }
}

//** GTF file outputing **//
void GeneModel::formatGTF (ofstream &forf){
  string strdSt = (_strand < 0) ? "-" : "+";

  s32 leftB = _cds.second, rightB = _cds.first;
  if (_cds.first < _cds.second){  
    leftB = _cds.first; 
    rightB = _cds.second;
  }
  
  s32 debTrans = (*_exons.begin())->start();
  s32 finTrans = (*(--_exons.end()))->end();
  if (_strand >= 0) {
    debTrans -= _5pXL;
    finTrans += _3pXL;
  }
  else if (_strand < 0) {
    finTrans += _5pXL;
    debTrans -= _3pXL;
  }
    

  string argCDS;
  if(this->containCDS()) {
    ostringstream oss;
    oss << " start " << _start_found << " ; stop " << _stop_found << " ; cds_size " << _cds_size;
    argCDS = oss.str();	  
  }  
  f8 score = this->score();

  forf << _seqname << "\tGmove\ttranscript\t" << debTrans << "\t" << finTrans << "\t" << score << "\t" << strdSt << "\t.\t" << "gene_id.modele \""
       << _seqname << "." << _num_cc << "." << _indMod << "\"; transcript_id \"" << _seqname << "." << _num_cc << "." << _indMod  << "\";" << argCDS << endl;

  forf << _seqname << "\tGmove\tstart_codon\t" << _cds.first << "\t" << _cds.first+2 << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \""
       << _seqname << "." << _num_cc << "." << _indMod  << "\"; transcript_id \"" << _seqname << "." << _num_cc << "." << _indMod  << "\";" <<endl;

  forf << _seqname << "\tGmove\tstop_codon\t" << _cds.second-2 << "\t" << _cds.second << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \""
       << _seqname << "." << _num_cc << "." << _indMod  << "\"; transcript_id \"" << _seqname << "." <<_num_cc << "." << _indMod  << "\";" << endl;
  

  for (TSSRList::iterator it = _exons.begin(); it != _exons.end(); it++){
    s32 debExon = (*it)->start();
    s32 finExon = (*it)->end();
    if (_strand >=0) {
      if (it == _exons.begin()) {debExon -= _5pXL;}
      if (it == --_exons.end()) {finExon += _3pXL;}
    }
    else if(_strand < 0) {
      if (it == --_exons.end()) {finExon += _5pXL;}
      if (it == _exons.begin()) {debExon -= _3pXL;}
    }
    
    forf << _seqname << "\tGmove\texon\t" << debExon << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "."
	 <<_num_cc << "." << _indMod  << "\"; transcript_id \"" << _seqname << "." << _num_cc << "." << _indMod  << "\";" << endl;
    
    if(this->containCDS()) {
      if ((finExon >= leftB) && (debExon <= leftB) && (finExon <= rightB))
	forf << _seqname << "\tGmove\tCDS\t" << leftB << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "."
	     << _num_cc << "." << _indMod << "\"; transcript_id \"" << _seqname << "." << _num_cc << "." << _indMod  <<"\";" << endl;
      else if ((finExon >= leftB) && (debExon <= leftB))
	forf << _seqname << "\tGmove\tCDS\t" << leftB << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "."
	     << _num_cc << "." << _indMod << "\"; transcript_id \"" << _seqname << "." <<_num_cc << "." << _indMod << "\";" << endl;
      else if ((finExon <= rightB) && (debExon >= leftB))
	forf << _seqname << "\tGmove\tCDS\t" << debExon << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "."
	     << _num_cc << "." << _indMod << "\"; transcript_id \"" << _seqname << "." << _num_cc << "." << _indMod  << "\";" << endl;
      else if ((finExon >= rightB) && (debExon <= rightB))
	forf << _seqname << "\tGmove\tCDS\t" << debExon<< "\t" << rightB << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "."
	     << _num_cc << "." << _indMod  << "\"; transcript_id \"" << _seqname << "." << _num_cc << "." << _indMod  << "\";" << endl;
    }
  }
}

// to launch the ORF selection
bool GeneModel::findORF(const map<string, s32>& protPhaseCoord) {
  bool oriented = true;
  /* first pass to find an ORF using the proteic mapping phase information (if available) */
  if( (_exons.size() == 1) || GeneModel::UNORIENTED ){
    if(_exons.front()->strand() == SSRContig::UNKNOWN){
      oriented = false;
      (_exons.front())->setStrand(SSRContig::FORWARD);
    }
  }
  s32 phaseScore_fwd = this->selectORF(protPhaseCoord); // selectORF update _cds
  //cout << " select orf based in protPhaseCoord " << phaseScore_fwd << endl;
  if( (_exons.size() == 1) || GeneModel::UNORIENTED ){
    //	  cout << "if( (_exons.size() == 1) || GeneModel::UNORIENTED )"<<endl;
    s32 _3pXL_fwd =_3pXL, _5pXL_fwd=_5pXL;
    //  _3pXL = 0;//XXX Why ??
    //  _5pXL = 0;
    pair<s32, s32> cds_forward = _cds;
    s32 cds_forward_len = _cds_size, start_found = _start_found, stop_found = _stop_found;
    
    s32 phaseScore_rev = 0;
    if(_exons.front()->strand() == SSRContig::UNKNOWN || oriented == false){
      //  	cout << "if(_exons.front()->strand() == SSRContig::UNKNOWN || oriented == false)"<<endl;
      (_exons.front())->setStrand(SSRContig::REVERSE);
      _strand = (_exons.front())->strand();
      phaseScore_rev = this->selectORF(protPhaseCoord); // comparison between fwd and rev models to choose the best one
      
      if((_ref_score && phaseScore_rev < phaseScore_fwd) || (!_ref_score && _cds_size <= cds_forward_len)) {
	(_exons.front())->setStrand(SSRContig::FORWARD);
	_strand = (_exons.front())->strand();
	_cds = cds_forward;
	_cds_size = cds_forward_len;
	_start_found = start_found;
	_stop_found = stop_found;
	_3pXL = _3pXL_fwd;
	_5pXL = _5pXL_fwd;		  
	//	  cout << "cds " << _cds.first << " " << _cds.second << " cds_size " << _cds_size << endl;
	
      }
    }
  } 
  /* second pass if no ORF was found using the proteic mapping phase information --> information not used anymore */
  if(_ref_score && _cds_size==0){
    //  cout << " if(_ref_score && _cds_size==0)"<<endl;
    _ref_score = 0;
    if( (_exons.size() == 1) || GeneModel::UNORIENTED )
      if(_exons.front()->strand() == SSRContig::UNKNOWN)
	(_exons.front())->setStrand(SSRContig::FORWARD);
    this->selectORF();
    if( (_exons.size() == 1) || GeneModel::UNORIENTED ){
      //   	cout << " if( (_exons.size() == 1) || GeneModel::UNORIENTED )"<<endl;
      s32 _3pXL_fwd =_3pXL, _5pXL_fwd=_5pXL;
      //      _3pXL = 0; //XXX Why ?
      //      _5pXL = 0;
      pair<s32, s32> cds_forward = _cds;
      s32 cds_forward_len = _cds_size, start_found = _start_found, stop_found = _stop_found;
      if(_exons.front()->strand() == SSRContig::UNKNOWN){
	//   	  cout <<" if(_exons.front()->strand() == SSRContig::UNKNOWN)"<<endl;
	_3pXL = 0;
	_5pXL = 0;
	(_exons.front())->setStrand(SSRContig::REVERSE);
	_strand = (_exons.front())->strand();
	this->selectORF();  // comparison between fwd and rev models to choose the longest
      }
      
      if(_cds_size < cds_forward_len){
	//  	  cout<<" if(_cds_size < cds_forward_len)"<<endl;
	(_exons.front())->setStrand(SSRContig::FORWARD);
	_strand = (_exons.front())->strand();
	_cds = cds_forward;
	_cds_size = cds_forward_len;
	_start_found = start_found;
	_stop_found = stop_found;
	_3pXL = _3pXL_fwd;
	_5pXL = _5pXL_fwd;
	//		 cout << "cds " << _cds.first << " " << _cds.second << " cds_size " << _cds_size << endl;
      }
    }
  }
  
  if(_cds.first && _cds.second ) { // if a cds exist
    this->mapORF();
    if(_cds_size %3 != 0){
      cerr << "Error cds size is not %3 " << _cds_size<< " start found " << _start_found << " stop_found " << _stop_found << endl;
      exit(1);
    }
    _score_prot = (_cds_size != 0) ? 0 : (float)phaseScore_fwd / (float)_cds_size;
    _score_lencds = (this->getMrnaLen() != 0) ? 0 : (float)_cds_size / (float)this->getMrnaLen();
    return true;
  }
  return false;
}

// ORF finding 
s32 GeneModel::selectORF (const map<string, s32>& protPhaseCoord) {
  string exon_seq, mrna_seq, seq_name;	
  seq_name = this->getSeqname();
  s32 strand=0, refScore=0, mrna_coord=0;
  map<string, s32> modelPhaseCoord;
  /* recreate the mRNA sequence of the model and calculate the reference score for phase information */
  for(TSSRList::iterator it = _exons.begin() ; it != _exons.end(); it++) {
    if(!strand) { strand = (*it)->strand(); }
    //   cout << "selectORf pos "<< (*it)->start() << " "<< (*it)->end() << endl;
    exon_seq = (*it)->getSeq();
    mrna_seq += exon_seq;
    //   cout << "mra_seq length " << mrna_seq.size() << endl;
    if(!_ref_score) {
      for(s32 cursor=(*it)->start(); cursor<=(*it)->end();cursor++){
	mrna_coord++;
	ostringstream oss1, oss2;
	oss1 << seq_name << "@" << cursor; //absolute coordinate on genomic sequence
	oss2 << seq_name << "@" << mrna_coord; //relative coordinate on mrna sequence
	string key1 = oss1.str(), key2 = oss2.str();
	
	map<string, s32>::const_iterator itPhase = protPhaseCoord.find(key1);
	s32 t=0;
	if(itPhase != protPhaseCoord.end()) { // && itPhase->second){
	  refScore+=1;//alignment proteique exist ?
	  t++;
	  modelPhaseCoord.insert(make_pair(key2, itPhase->second)); //careful: strand<0 coordinates in the map are stored in forward
	}
	//cerr << "TEST key: " << key1 << " vu: " << t << " second: " << itPhase->second << " genomic: " << key1 << endl;
      }
    }
  }
  //cerr << "REFSCORE ==> " << refScore << endl;
  if(refScore) {
    _ref_score = refScore;
    _protPhaseCoord = modelPhaseCoord;
  }
  if(modelPhaseCoord != _protPhaseCoord) modelPhaseCoord = _protPhaseCoord;
  transform(mrna_seq.begin(), mrna_seq.end(), mrna_seq.begin(), (int(*)(int))toupper);
  /*reverse the sequence if the strand is - */
  if (strand < 0) { 
    reverse(mrna_seq.begin(),mrna_seq.end());
    transform(mrna_seq.begin(), mrna_seq.end(), mrna_seq.begin(), complem);
  }
  /*define the start and stop codons*/
  vector<string> startC, stopC;
  if(_genetic_code == 6) {
    //cerr << "genetic code 6" << endl;
    startC.push_back("ATG");
    stopC.push_back("TGA");
  }
  if(_genetic_code == 23) {
    //cerr << "genetic code 23" << endl;
    startC.push_back("ATG");
    startC.push_back("ATT");
    startC.push_back("GTG");
    stopC.push_back("TAA");
    stopC.push_back("TTA");
    stopC.push_back("TAG");
    stopC.push_back("TGA");
  }
  if(_genetic_code != 6 && _genetic_code != 23) { 
    //cerr << "genetic code 1" << endl;
    startC.push_back("ATG");
    stopC.push_back("TAA");
    stopC.push_back("TAG");
    stopC.push_back("TGA");
  }
  vector<string>::iterator it;
  vector<s32> start, stop;
  
  // search start codons
  for ( it=startC.begin() ; it < startC.end(); it++ ) {
    string::size_type loc = mrna_seq.find(*it, 0);
    while( loc != string::npos) {
      start.push_back((s32)loc);
      loc = mrna_seq.find(*it, loc+1);
    }
  }
  // search stop codons
  for ( it=stopC.begin() ; it < stopC.end(); it++ ) {
    string::size_type loc = mrna_seq.find(*it, 0);
    while( loc != string::npos) {
      stop.push_back((s32)loc);
      loc = mrna_seq.find(*it, loc+1);
    }
  }
  std::sort(start.begin(), start.end());
  std::sort(stop.begin(), stop.end());
  
  pair<s32,s32> cds;
  s32 cdslen = 0, phaseScore = 0;
  vector<s32>::iterator itStart, itStop;
  // try to find a full orf
  for (s32 i = 0 ; i < 3 ; i++) {
    //	  cout <<i << " full ORF ";
    bool findStop = false;
    for ( itStop = stop.begin() ; itStop < stop.end(); itStop++ ) {
      s32 pos = *itStop;
      if((pos-i)%3 == 0) {
	findStop = true;
	break; 
      }
    }
    if(!findStop) {
      //  	cout <<"no stop in phase " << endl;
      s32 start = i+1;
      s32 end = (s32)mrna_seq.length();
      while ((end-start +1)%3 != 0 ){
	end = end-1;
      }
      //   	cout << "start " << start << " end " << end << endl;
      pair<s32,s32> orf_begin_end = make_pair(start,end);
      s32 len_begin_end = orf_begin_end.second - orf_begin_end.first + 1;
      //     cout << "len_begin_end "<< len_begin_end << endl;
      if(_ref_score) {
	//   	  cout <<"enter in _ref_score " << endl;
	s32 score_begin_end = this->getPhaseScore(orf_begin_end, strand, mrna_seq.length(), modelPhaseCoord);
	if(score_begin_end > phaseScore) {
	  phaseScore = score_begin_end;
	  cds = orf_begin_end;
	  cdslen = len_begin_end;
	  _stop_found = 0;
	  _start_found = 0;
	  //cerr << "JM : found full ORF ; start: " << orf_begin_end.first << " stop: " << orf_begin_end.second << endl;
	}
      }
      else if(len_begin_end > cdslen) {
	cds = orf_begin_end;
	cdslen = len_begin_end;
	_stop_found = 0;
	_start_found = 0;
      }
    }
    //   cout <<"phase score " << phaseScore<< "cdslen "<< cdslen << endl;
  }
  
  
  // find the longest M->* orf
  pair<s32,s32> orf_M_stop;
  s32 len_M_stop;
  if(_ref_score) {
    //	  cout <<"M * orf ";
    s32 score_M_stop;
    orf_M_stop = this->bestORF(start, stop, strand, mrna_seq.length(), score_M_stop,modelPhaseCoord);
    len_M_stop = orf_M_stop.second - orf_M_stop.first + 1;
    if(score_M_stop > phaseScore) {
      phaseScore = score_M_stop;
      cds = orf_M_stop;
      cdslen = len_M_stop;
      _stop_found = 1;
      _start_found = 1;
      //cerr << "JM : found M -> * ; start: " << orf_M_stop.first << " stop: " << orf_M_stop.second << " score: " << score_M_stop << endl;
    }
    //   cout << "phase score "<< phaseScore << "cdslen "<< cdslen << endl;
  }  
  else {
    orf_M_stop = this->longestORF(start, stop);
    len_M_stop = orf_M_stop.second - orf_M_stop.first + 1;
    if(len_M_stop > cdslen) {
      cds = orf_M_stop;
      cdslen = len_M_stop;     
      _stop_found = 1;
      _start_found = 1;
    }  
  }
  // find the longest M->end orf
  vector<s32> allStop(stop);
  allStop.push_back(mrna_seq.length()-3); //FIXME Why 3, 4 , 5 ?
  allStop.push_back(mrna_seq.length()-4);
  allStop.push_back(mrna_seq.length()-5);
  std::sort(allStop.begin(), allStop.end());
  pair<s32,s32> orf_M_end;
  s32 len_M_end;
  
  if(_ref_score) {
    //	  cout << "M end "<<endl;
    s32 score_M_end;
    orf_M_end = this->bestORF(start, allStop, strand, mrna_seq.length(), score_M_end, modelPhaseCoord);
    len_M_end = orf_M_end.second - orf_M_end.first + 1;
    if(score_M_end > phaseScore) {
      phaseScore = score_M_end;
      cds = orf_M_end;
      cdslen = len_M_end;
      _start_found = 1;
      _stop_found = 0;
    }
    //    cout <<"phaseScore " << phaseScore << " cdslen " << cdslen << endl;
  }
  else {
    orf_M_end = this->longestORF(start, allStop);
    len_M_end = orf_M_end.second - orf_M_end.first + 1;
    if(len_M_end > cdslen) {
      cds = orf_M_end;
      cdslen = len_M_end;
      _start_found = 1;
      _stop_found = 0;
    }
  }
  // find the longest begin->* on a different frame that the existing longest cds, if the begin->* segment is at least 64bp longer
  vector<s32> allStart(start);
  for(s32 i = 0; i < 3; i++) if((i-cds.first+1)%3 !=0) allStart.push_back(i);
  std::sort(allStart.begin(), allStart.end());
  pair<s32,s32> orf_begin_stop;
  s32 len_begin_stop; 
  
  if(_ref_score) {
    //	  cout << "begin * " << endl;
    s32 score_begin_stop;
    orf_begin_stop = this->bestORF(allStart, allStop, strand, mrna_seq.length(), score_begin_stop, modelPhaseCoord);
    len_begin_stop = orf_begin_stop.second - orf_begin_stop.first + 1;
    if(score_begin_stop > phaseScore){
      phaseScore = score_begin_stop;
      cds = orf_begin_stop;
      cdslen = len_begin_stop;
      _start_found = 0;
      _stop_found = 1;
    }
    //  cout << "phaseScore " << phaseScore << " cdslen " << cdslen <<endl;
  } 
  else {
    orf_begin_stop = this->longestORF(allStart, stop);
    len_begin_stop = orf_begin_stop.second - orf_begin_stop.first + 1;
    if(len_begin_stop > cdslen + 64) {
      cds = orf_begin_stop;
      cdslen = len_begin_stop;
      _start_found = 0;
      _stop_found = 1;
    }
  }
  
  // SECOND STEP IF NO STOP WAS FOUND, TAKE THE FIRST ONE FOUND IN THE SAME FRAME, AT A MAX DISTANCE OF _WINDOW
  s32 frame;
  for (frame = 0; frame < 3; frame++){
    if ((cds.first-1-frame)%3==0)
      break;
    
  }
  string mrna_xlarge, mrna_3p, mrna_5p;
  string* whole_seq = (*_exons.begin())->getWholeSeq();
  s32 start_xlarge, size_xlarge; 
  s32 newStopPos = -1, newStartPos = -1;
  
  
  if (_stop_found == 0 && cds!=make_pair(0,0)) {
    //obtain 3' extended sequence
    if (strand >= 0) {
      start_xlarge =  min((*(--_exons.end()))->end(), (s32)whole_seq->length()-1);
      mrna_3p = whole_seq->substr(start_xlarge, _window);
      transform(mrna_3p.begin(), mrna_3p.end(), mrna_3p.begin(), (int(*)(int))toupper);
      mrna_xlarge = mrna_seq + mrna_3p;
    }
    
    else { // strand < 0
      start_xlarge = max(0,(*_exons.begin())->start() -1 -_window);
      size_xlarge = (*_exons.begin())->start() - start_xlarge -1;
      mrna_3p = whole_seq->substr(start_xlarge, size_xlarge);
      transform(mrna_3p.begin(), mrna_3p.end(), mrna_3p.begin(), (int(*)(int))toupper);
      reverse(mrna_3p.begin(),mrna_3p.end());
      transform(mrna_3p.begin(), mrna_3p.end(), mrna_3p.begin(), complem);
      mrna_xlarge = mrna_seq + mrna_3p;
      
    }
    
    // find STOP codon
    for ( it=stopC.begin() ; it < stopC.end(); it++ ) {
      string::size_type loc = mrna_xlarge.find(*it, 0);
      while( loc != string::npos) {
	if ((((s32)loc-frame)%3==0) && ((s32)loc+3-cds.first+1 >= cdslen)) {
	  if((s32)loc<newStopPos || newStopPos==-1) {newStopPos = (s32)loc;}
	  break;}
	loc = mrna_xlarge.find(*it, loc+1);
      }
    }
    if(newStopPos != -1) {
      cds = make_pair(cds.first, newStopPos + 3);
      cdslen = cds.second - cds.first + 1;
      _stop_found = 1;
      _3pXL = newStopPos + 3 - (s32)mrna_seq.length();
    }
  }
  
  
  // THIRD STEP IF NO START WAS FOUND, TAKE THE CLOSEST ONE TO THE BEGINNING OF THE CDS IN THE LIMIT OF (-_WINDOW,+_WINDOW)
  if(_start_found == 0 && cds!=make_pair(0,0)) {
    s32 distance = -1, distemp = 0;
    // obtain 5' extended sequence
    if (strand >= 0) {
      start_xlarge = max(0,(*_exons.begin())->start() -1 -_window);
      size_xlarge = (*_exons.begin())->start() - start_xlarge -1;
      mrna_5p = whole_seq->substr(start_xlarge, size_xlarge);
      transform(mrna_5p.begin(), mrna_5p.end(), mrna_5p.begin(), (int(*)(int))toupper);
      mrna_xlarge = mrna_5p + mrna_seq;
    }
    
    else { // strand < 0
      start_xlarge =  min((*(--_exons.end()))->end(), (s32)whole_seq->length()-1);
      mrna_5p = whole_seq->substr(start_xlarge, _window);
      transform(mrna_5p.begin(), mrna_5p.end(), mrna_5p.begin(), (int(*)(int))toupper);
      reverse(mrna_5p.begin(),mrna_5p.end());
      transform(mrna_5p.begin(), mrna_5p.end(), mrna_5p.begin(), complem);
      mrna_xlarge = mrna_5p + mrna_seq;
    }
    
    // find STOP codon in the enlarged sequence
    stop.clear();
    
    for ( it=stopC.begin() ; it < stopC.end(); it++ ) {
      string::size_type loc = mrna_xlarge.find(*it, 0);
      while( loc != string::npos && ((s32)loc < (s32)mrna_5p.length()+(s32)cds.first)) {
	if (((s32)loc-(s32)mrna_5p.length()-frame)%3==0) {stop.push_back((s32)loc); }
	loc = mrna_xlarge.find(*it, loc+1);
      }
    } 
    
    std::sort(stop.begin(), stop.end());
    s32 lastStop;
    if (stop.empty()) lastStop=0;
    else lastStop = stop.back()+3;
    
    // find START codon
    for ( it=startC.begin() ; it < startC.end(); it++ ) {
      string::size_type loc = mrna_xlarge.find(*it, lastStop);
      while( loc != string::npos) {
	// check this is the right ORF
	if((((s32)loc-(s32)mrna_5p.length()-frame)%3==0)){
	  // check the START is the closest to the initial cds_start
	  distemp = abs((int)((s32)mrna_5p.length()+cds.first-1-(s32)loc));
	  if(distemp<=_window && (distemp < distance || distance == -1)) {
	    distance = distemp;
	    newStartPos = (s32)loc; // we only keep the closest one
	  }
	}
	loc = mrna_xlarge.find(*it, loc+1);
      }
    }
    if (newStartPos != -1) {
      _start_found = 1;
      if(newStartPos+1 <= (s32)mrna_5p.length()) {  //the new START codon was found 5' of the old beginning of CDS
	_5pXL = (s32)mrna_5p.length()-newStartPos;
	cds = make_pair(1, cds.second + _5pXL);
      }
      else {
	_5pXL = 0;
	cds = make_pair(newStartPos-(s32)mrna_5p.length()+1,cds.second);
      }
      cdslen = cds.second - cds.first + 1;
    }
  }
  if(cdslen != 0) _cds = cds;
  else _cds = make_pair(0, 0);
  _cds_size = cdslen;
  //cerr << "PHASE SCORE: " << phaseScore << endl;
  return phaseScore;
}

// to calculate the score of an ORF based on the proteic mapping phase information
s32 GeneModel::getPhaseScore(pair<s32, s32> orf, s32 strand, s32 mrna_len, const map<string, s32>& modelPhaseCoord){
  s32 score=0, rev_j=0, nt=0;
  for(s32 j = orf.first; j <= orf.second; j++){
    ostringstream oss;
    if(strand > 0) { oss << this->getSeqname() << "@" << j; }
    if(strand < 0) {
      rev_j = mrna_len - j + 1;
      oss << this->getSeqname() << "@" << rev_j;
    }
    string key = oss.str();
    map<string, s32>::const_iterator it = modelPhaseCoord.find(key);
    nt++;
    if(it != modelPhaseCoord.end() && it->second){
      s32 phase_diff;
      if(strand > 0) { phase_diff = (j - orf.first)%3 +1 - it->second; }
      else { phase_diff = (j - orf.first)%3 +1 + it->second; }
      if(!phase_diff) score++;      
    }
  }
  //s32 debTrans = (*_exons.begin())->start();
  //s32 finTrans = (*(--_exons.end()))->end();
  
  //float s = (float)score / (float)nt;
  //float l = (float)nt / (float)mrna_len;
  //cerr << " Start= " << debTrans << " stop= " << mrna_len << " NT= " << nt << " score= " << score << " ratio= " << s << " ratiolen= " << l << " scorefinal= " << l*s << endl;
  return score;
}

// to find the longest ORF between all combinaisons of start and stop
pair<s32,s32> GeneModel::longestORF(vector<s32> lStart, vector<s32> lStop) {
//	std::clock_t c_start = std::clock();
  vector<s32>::iterator itStart, itStop;
  pair<s32,s32> orf;
  s32 max_len = 0;
  for ( itStart = lStart.begin() ; itStart < lStart.end(); itStart++ ) {
    for ( itStop = lStop.begin() ; itStop < lStop.end(); itStop++ ) {
      s32 pstart = *itStart, pstop = *itStop;
      s32 cds_len = pstop - pstart + 1 + 2;
      if(cds_len <=0 || cds_len%3 != 0) { continue; }
      else if(cds_len > max_len) {
	orf = make_pair(pstart + 1, pstop + 3);
	max_len = cds_len;
      }
      itStop = lStop.end();
    }
  }
  //	std::clock_t c_end = std::clock();
  //	 if(SSRContig::VERBOSE)   cout << " time longestORF " <<  1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << endl;
  return orf;
}

// to find the ORF between all combinaisons of start and stop that has the best phase score
pair<s32,s32> GeneModel::bestORF(vector<s32> lStart, vector<s32> lStop, s32 strand, s32 mrna_len, s32& score,const map<string, s32>& modelPhaseCoord ){
  //std::clock_t c_start = std::clock();
  vector<s32>::iterator itStart, itStop;
  pair<s32,s32> orf, tmp_orf;
  s32 max_score=0;
  for ( itStart = lStart.begin() ; itStart < lStart.end(); itStart++ ) {
    for ( itStop = lStop.begin() ; itStop < lStop.end(); itStop++ ) {
      s32 pstart = *itStart, pstop = *itStop;
      s32 cds_len = pstop - pstart + 1 + 2;
      if(cds_len <=0 || cds_len%3 != 0) { continue; }
      tmp_orf = make_pair(pstart + 1, pstop + 3);
      score = this->getPhaseScore(tmp_orf, strand, mrna_len, modelPhaseCoord);
      if(score > max_score) {
	orf = tmp_orf;
	max_score = score;
      }
      itStop = lStop.end();
    }
  }
  score = max_score;
  //std::clock_t c_end = std::clock();
  //if(SSRContig::VERBOSE) cout << " time bestORF " <<  1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << endl;
  return orf;
}

// to recalculate the absolute coordinate of the ORF on the genomic sequence
void GeneModel::mapORF() {
  s32 debOrf = _cds.first, finOrf = _cds.second;
  s32 posEp=0; // exonic sequence pos
  s32 posG=0; // genome position
  s32 orfgin=0, orfgout=0;
  if (_strand >= 0) { // if forward or unknown
    posEp=1;
    for(TSSRList::iterator it = _exons.begin() ; it != _exons.end(); it++) {
      s32 debExon = (*it)->start();
      s32 finExon = (*it)->end();
      if(it == --_exons.end()) finExon = min((finExon + _window), (s32)(*it)->getWholeSeq()->length());
      if(it == _exons.begin()) debExon-=_5pXL;
      for(posG = debExon ; posG <= finExon ; posG++){
    	  if (posEp == debOrf) orfgin=posG;
    	  if (posEp == finOrf) orfgout=posG;
    	  if (posEp > finOrf) break;
    	  posEp++;
      }
    }
    _cds = make_pair(orfgin,orfgout);
  }
  else { // _strand < 0
    posEp=1;
    s32 debExon = 0;
    for(TSSRList::reverse_iterator rit = _exons.rbegin() ; rit != _exons.rend(); rit++) {
       debExon = (*rit)->start();
      s32 finExon = (*rit)->end();
      if(rit == --_exons.rend()) debExon -= _3pXL;
      if(rit == _exons.rbegin()) finExon += _5pXL;
      for(posG = finExon ; posG >= debExon ; posG--){
    	  if (posEp == debOrf) orfgin=posG;
    	  if (posEp == finOrf) orfgout=posG;
    	  if (posEp > finOrf) break;
    	  posEp++;
      }
    }
 //   cout << "orfgout " << orfgout << " orfgin " << orfgin << endl;
      _cds = make_pair(orfgout,orfgin);
    if(!this->containCDS())
    	cerr << "error in mapORF _cds.first " << _cds.first << " _cds.second "<< _cds.second <<endl;//FIXME THAT APPEND SOMETIMES !!
    }
}

// to return the exon corresponding to the id (the exons list in GeneModel only has access to IDs)
SSRContig* GeneModel::exon(s32 id) {
  for (TSSRList::iterator iEx = _whole_exons->begin(); iEx != _whole_exons->end(); iEx++) if ((*iEx)->getID() == id) return *iEx;
  cerr << "[Error] exon not found id= " << id << endl;
  assert(false);
} 


// complementary strand function
char complem (char base){
  if      (base=='A') return 'T';
  else if (base=='T') return 'A';
  else if (base=='G') return 'C';
  else if (base=='C') return 'G';
  else { return base;}
}

void selectModel( map<string,GeneModel> mapGeneModel){
  for(map<string,GeneModel>::iterator itMap = mapGeneModel.begin(); itMap != mapGeneModel.end();++itMap) {
    for(map<string,GeneModel>::iterator itMap2 = mapGeneModel.begin(); itMap2 != mapGeneModel.end();++itMap2) {
      if(itMap == itMap2) continue;
      if(itMap->second.getStartFound() && itMap2->second.getStartFound() 
	 && itMap->second.getStopFound() && itMap2->second.getStopFound())continue;
    }
  }
}

bool GeneModel::cdsIsMono(){
  if(_exons.size()==1)
    return true;
  list<pair<s32,s32 > > exonsCds = getExonsCds();
  if(exonsCds.size() == 1 )
    return true;
  else
    return false;
}

list<pair<s32,s32> > GeneModel::getExonsCds() {
  //FIXME Error in getExonsCds ! when there is UTR and CDS on the same exon !
  list<pair<s32,s32> > exonCDS1;
  s32 debCDS1 = _cds.second, finCDS1 = _cds.first;
  if (_cds.first < _cds.second) {
    debCDS1 = _cds.first;
    finCDS1 = _cds.second;
  }
  for (TSSRList::iterator it = _exons.begin(); it != _exons.end(); ++it){
    if(_exons.size()==1 ){
      exonCDS1.push_back(_cds);
      break;
    }
    if(_cds.first >= (*it)->start() && _cds.second <= (*it)->end()){
      exonCDS1.push_back(_cds);
      break;
    }
    if( ( (*it)->start() >= _cds.first && (*it)->start() <=  _cds.second ) || ( (*it)->end()>= _cds.first &&  (*it)->end() <= _cds.second ) ){
      s32 debExon1 = (*it)->start();
      s32 finExon1 = (*it)->end();
      if( finCDS1 < debExon1 || debCDS1 > finExon1 )
	continue;
      else {
	s32 end = (finCDS1 < finExon1) ? finCDS1 : finExon1;
	s32 begin = (debCDS1 < debExon1) ? debExon1 : debCDS1;
	//	cout << "exon CDS " << begin << " " << end << endl;
	exonCDS1.push_back(make_pair(begin,end));
      }
    }
  }
  return exonCDS1;
}


void GeneModel::keepModel(map<string,GeneModel>& mapGeneModel){
  ostringstream oss;
  s32 debCDS = _cds.second, finCDS = _cds.first;
  if (_cds.first < _cds.second) {
    debCDS = _cds.first;
    finCDS = _cds.second;
  }
  oss << _seqname;
  for (TSSRList::iterator it = _exons.begin(); it != _exons.end(); it++){
    s32 debExon = (*it)->start();
    s32 finExon = (*it)->end();
    if (_strand >=0) {
      if (it == _exons.begin()) debExon -= _5pXL;
      if (it == --_exons.end()) finExon += _3pXL;
    }
    else if(_strand < 0) {
      if (it == --_exons.end()) finExon += _5pXL;
      if (it == _exons.begin()) debExon -= _3pXL;
    }
    if( finCDS < debExon || debCDS > finExon )
      continue;
    
    else {
      s32 end = (finCDS < finExon) ? finCDS : finExon;
      s32 begin = (debCDS < debExon) ? debExon : debCDS;
      oss <<"@"<< begin << "@" << end ;
    }
  }
  string key = oss.str();
  map<string,GeneModel>::iterator itMap = mapGeneModel.find(key);
  if( itMap == mapGeneModel.end()){
    pair<string,GeneModel> tmpPair = make_pair(key,*this);
    mapGeneModel.insert(tmpPair);
  }
  else{
    if(this->_model_size >= itMap->second._model_size)
      itMap->second= *this;
    //else do nothing
  }
}

bool GeneModel::overlapOrf(pair<s32,s32> cds2){
  return ((this->_cds.first <= cds2.first && this->_cds.second >= cds2.first) || 
	  (this->_cds.first >= cds2.first && this->_cds.first <= cds2.second));
}

bool GeneModel::overlapModel(GeneModel model2){
  s32 startModel1 =  (*(this->_exons.begin()))->start();
  s32 endModel1 =  (*(this->_exons.back())).end();
  s32 startModel2 =  (*(model2._exons.begin()))->start() ;
  s32 endModel2 = (*(model2._exons.back())).end();
  if (endModel1 >= startModel2  &&  endModel1 <= endModel2 )
    return true;
  else if( startModel1 >= startModel2 && startModel1 <= endModel2)
    return true;
  else if (startModel2 >= startModel1 && startModel2 <= endModel1)
    return true;
  else if(endModel2 >= startModel1 && endModel2 <= endModel1)
    return true;
  else
    return false;
}

bool GeneModel::overlapOrf(GeneModel model2){//overlap cds from model
  // overlap cluster location : look at all exons from both models and test if overlap 50 with a min contiguity of 10
  s32 startOverlap, endOverlap,sizeOverlap;
  s32 overlap = 0;
  list<pair<s32,s32> > listExons1 = this->getExonsCds();
  for(list<pair<s32,s32> >::iterator itListExons1 = listExons1.begin() ; itListExons1 != listExons1.end() ; ++itListExons1){
    s32 startExon1 = itListExons1->first;
    s32 endExon1 =  itListExons1->second;
    list<pair<s32,s32> > listExons2 = model2.getExonsCds();
    for(list<pair<s32,s32> >::iterator itListExons2 = listExons2.begin() ; itListExons2 != listExons2.end() ; ++itListExons2){
      s32 startExon2 = itListExons2->first;
      s32 endExon2 =  itListExons2->second;
      if( (startExon1 >= startExon2 && startExon1 <= endExon2 )|| (endExon1 >= startExon2 && endExon1 <= endExon2) || (startExon2 >= startExon1 && startExon2 <= endExon1) || (endExon2 >= startExon1 && endExon2 <= endExon1)){
	startOverlap = (startExon1 < startExon2) ? startExon2 : startExon1;
	endOverlap = (endExon1 < endExon2) ? endExon1 : endExon2;
	sizeOverlap = endOverlap - startOverlap +1;
	if(sizeOverlap > 10)
	  overlap += sizeOverlap;
	if(overlap > 50){
	  return true;
	  break;
	}
      }
    }
  }
  return false;
}


// to print a covtig's information
ostream& operator<<(ostream& ostr, const GeneModel& g) {
  return ostr<<g.getSeqname()<< " " << g.getNumCC() <<" " << g.getID() << " " << g.getCdsSize();
  // for (TSSRList::iterator it = g.getExons().begin(); it != g.getExons().end(); it++){
  //	ostr <<"exon " <<  (*it)->start() << " " << (*it)->end() <<endl;
  //}
  // return ostr;
}

