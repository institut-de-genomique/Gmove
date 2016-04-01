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
void GeneModel::printAnnot(ofstream &forf, s32 nbp, bool format) {
  if (format == true) this->formatGFF(forf, nbp);
  else this->formatGTF(forf, nbp);
}


// GFF3 file outputing: format fitting our in-house needs
void GeneModel::formatGFF (ofstream &forf, s32 nbp) {
  s32 idModel = getID();	
  string strdSt = (_strand < 0) ? "-" : "+";
  
  s32 debCDS = _cds.second, finCDS = _cds.first;
  if (_cds.first < _cds.second) {  
    debCDS = _cds.first; 
    finCDS = _cds.second;
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
    oss << ";start="<<_start_found<< ";stop=" << _stop_found << ";cds_size=" << _cds_size<<";model_size="<<_model_size+_5pXL+_3pXL;
    argCDS = oss.str();	  
  }

  forf << _seqname << "\tGmorse\tmRNA\t" << debTrans << "\t" << finTrans << "\t.\t" << strdSt << "\t.\t" << "ID=mRNA." 
       << _seqname << "." << idModel << "." << nbp 
       << ";Name=mRNA." << _seqname << "." << idModel << "." << nbp 
       << argCDS << endl;

  //cout << " idModel " << idModel << " nbp "<<nbp << " nb exons " << _exons.size() <<" this->containCDS() "<<this->containCDS() << endl;
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
    //    forf << _seqname << "\tGmorse\texon\t" << debExon << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA." 
    //	 << _seqname << "." << idModel << "." << nbp << endl; 	
    if(this->containCDS()) {
      if( finCDS < debExon || debCDS > finExon ) {
	forf << _seqname << "\tGmorse\tUTR\t" << debExon << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA." 
 	     << _seqname << "." << idModel << "." << nbp << endl;
      }

      else {
	s32 end = (finCDS < finExon) ? finCDS : finExon;
	s32 begin = (debCDS < debExon) ? debExon : debCDS;
	if(begin != debExon)
	  forf << _seqname << "\tGmorse\tUTR\t" << debExon << "\t" << begin-1 << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA." 
	       << _seqname << "." << idModel << "." << nbp << endl;
	forf << _seqname << "\tGmorse\tCDS\t" << begin << "\t" << end << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA." 
 	     << _seqname << "." << idModel << "." << nbp << endl;
	if(end != finExon)
	  forf << _seqname << "\tGmorse\tUTR\t" << end+1 << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "Parent=mRNA." 
 	       << _seqname << "." << idModel << "." << nbp << endl;
      }
    }
  }
}

//** GTF file outputing **//
void GeneModel::formatGTF (ofstream &forf, s32 nbp){
  s32 idModel = getID();	
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
  
  forf << _seqname << "\tGmorse\ttranscript\t" << debTrans << "\t" << finTrans << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" 
       << _seqname << "." << idModel << "." << nbp << "\"; transcript_id \"" << _seqname << "." << idModel << "." << nbp << "\";" << argCDS << endl; 

  forf << _seqname << "\tGmorse\tstart_codon\t" << _cds.first << "\t" << _cds.first+2 << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" 
       << _seqname << "." << idModel << "." << nbp << "\"; transcript_id \"" << _seqname << "." << idModel << "." << nbp << "\";" <<endl; 

  forf << _seqname << "\tGmorse\tstop_codon\t" << _cds.second-2 << "\t" << _cds.second << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" 
       << _seqname << "." << idModel << "." << nbp << "\"; transcript_id \"" << _seqname << "." << idModel << "." << nbp << "\";" << endl; 
  

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
    
    forf << _seqname << "\tGmorse\texon\t" << debExon << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "." 
	 <<idModel << "." << nbp << "\"; transcript_id \"" << _seqname << "." << idModel << "." << nbp << "\";" << endl; 	
    
    if(this->containCDS()) {
      if ((finExon >= leftB) && (debExon <= leftB) && (finExon <= rightB))
	forf << _seqname << "\tGmorse\tCDS\t" << leftB << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "." 
	     << idModel << "." << nbp << "\"; transcript_id \"" << _seqname << "." << idModel << "." << nbp <<"\";" << endl; 	
      else if ((finExon >= leftB) && (debExon <= leftB))
	forf << _seqname << "\tGmorse\tCDS\t" << leftB << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "."
	     << idModel << "." << nbp << "\"; transcript_id \"" << _seqname << "." << idModel << "." << nbp << "\";" << endl; 	
      else if ((finExon <= rightB) && (debExon >= leftB))
	forf << _seqname << "\tGmorse\tCDS\t" << debExon << "\t" << finExon << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "."
	     << idModel << "." << nbp << "\"; transcript_id \"" << _seqname << "." << idModel << "." << nbp << "\";" << endl; 	
      else if ((finExon >= rightB) && (debExon <= rightB))
	forf << _seqname << "\tGmorse\tCDS\t" << debExon<< "\t" << rightB << "\t.\t" << strdSt << "\t.\t" << "gene_id.modele \"" << _seqname << "."
	     << idModel << "." << nbp << "\"; transcript_id \"" << _seqname << "." << idModel << "." <<nbp << "\";" << endl; 	
    }
  }
}

// to launch the ORF selection
bool GeneModel::findORF(const map<string, s32>& protPhaseCoord) {
  /* first pass to find an ORF using the proteic mapping phase information (if available) */
//	for(TSSRList::iterator itE = _exons.begin(); itE != _exons.end();++itE){
	//	cout << " exons in findOrf " << *(*itE) << endl;
//	}
//	cout << " exon strand " << _exons.front()->strand()<<endl;
  if( (_exons.size() == 1) || GeneModel::UNORIENTED ){
	  if(_exons.front()->strand() == SSRContig::UNKNOWN){
		  (_exons.front())->setStrand(SSRContig::FORWARD);
	  }

  }

  s32 phaseScore_fwd = this->selectORF(protPhaseCoord); // selectORF update _cds
  if( (_exons.size() == 1) || GeneModel::UNORIENTED ){
//	  cout <<"if _exons.size() == 1 GeneModel::UNORIENTED "<<endl;
    s32 _3pXL_fwd =_3pXL, _5pXL_fwd=_5pXL;

 //   cout << " _3pXL " <<_3pXL<< " _3pXL_fwd "<<_3pXL_fwd<<endl;
  //  _3pXL = 0;//XXX Why ??
  //  _5pXL = 0;
    pair<s32, s32> cds_forward = _cds;
 //   cout << " put _3pXL at 0 , cds.first " << _cds.first << " cds.second " << _cds.second << endl;
    s32 cds_forward_len = _cds_size, start_found = _start_found, stop_found = _stop_found;
    
    s32 phaseScore_rev = 0;
    if(_exons.front()->strand() == SSRContig::UNKNOWN){
    	(_exons.front())->setStrand(SSRContig::REVERSE);
    	     _strand = (_exons.front())->strand();
    	    phaseScore_rev = this->selectORF(protPhaseCoord); // comparison between fwd and rev models to choose the best one

   // 	cout << " ref score " << _ref_score << endl;
 //   	    cout << " _cds_size " << _cds.first << " " << _cds.second << endl;
		if((_ref_score && phaseScore_rev < phaseScore_fwd) || (!_ref_score && _cds_size <= cds_forward_len)) {
		  (_exons.front())->setStrand(SSRContig::FORWARD);
		  _strand = (_exons.front())->strand();
		  _cds = cds_forward;
		  _cds_size = cds_forward_len;
		  _start_found = start_found;
		  _stop_found = stop_found;
		  _3pXL = _3pXL_fwd;
		  _5pXL = _5pXL_fwd;
//		  cout<<" rechange _3pXL " << _3pXL << endl;
		}
    }
  } 
  /* second pass if no ORF was found using the proteic mapping phase information --> information not used anymore */
  if(_ref_score && _cds_size==0){
//	  cout << "if _ref_score && _cds_size==0 "<<endl;
    _ref_score = 0;
    if( (_exons.size() == 1) || GeneModel::UNORIENTED )
    	  if(_exons.front()->strand() == SSRContig::UNKNOWN)
    		  (_exons.front())->setStrand(SSRContig::FORWARD);
   // cout << " select orf " << endl;
    this->selectORF();
    if( (_exons.size() == 1) || GeneModel::UNORIENTED ){
      s32 _3pXL_fwd =_3pXL, _5pXL_fwd=_5pXL;
 //     cout << " _3pXL " <<_3pXL<< " _3pXL_fwd "<<_3pXL_fwd<<endl;
//      _3pXL = 0;
//      _5pXL = 0;
      pair<s32, s32> cds_forward = _cds;
      s32 cds_forward_len = _cds_size, start_found = _start_found, stop_found = _stop_found;
      if(_exons.front()->strand() == SSRContig::UNKNOWN){
    	  _3pXL = 0;
          _5pXL = 0;
    	  (_exons.front())->setStrand(SSRContig::REVERSE);
    	  _strand = (_exons.front())->strand();
    	  this->selectORF();  // comparison between fwd and rev models to choose the longest
      }

      if(_cds_size < cds_forward_len){
		(_exons.front())->setStrand(SSRContig::FORWARD);
		_strand = (_exons.front())->strand();
		_cds = cds_forward;
		_cds_size = cds_forward_len;
		_start_found = start_found;
		_stop_found = stop_found;
		_3pXL = _3pXL_fwd;
		_5pXL = _5pXL_fwd;
      }
    }
  }
//  cout << " _cds.first " << _cds.first << " cds.second " << _cds.second<<" "<< this->containCDS() <<endl;
  if(_cds.first && _cds.second ) { // if a cds exist
//	  cout << "if _cds.first && _cds.second  "<<endl;
    this->mapORF();
//    cout << " _cds.first " << _cds.first << " cds.second " << _cds.second<<" "<< this->containCDS() <<endl;
    return true;
  }
 // cout << " return false " << endl;
  return false;
}

// ORF finding 
s32 GeneModel::selectORF (const map<string, s32>& protPhaseCoord) {
  string exon_seq, mrna_seq, seq_name;	
  seq_name = this->getSeqname();
  s32 strand=0, refScore=0, mrna_coord=0;
  map<string, s32> modelPhaseCoord;
  map<string,bool> alreadyLookPos;
  /* recreate the mRNA sequence of the model and calculate the reference score for phase information */
  for(TSSRList::iterator it = _exons.begin() ; it != _exons.end(); it++) {
    if(!strand) { strand = (*it)->strand(); }
    exon_seq = (*it)->getSeq();
    mrna_seq += exon_seq;
    if(!_ref_score) {
      for(s32 cursor=(*it)->start(); cursor<=(*it)->end();cursor++){
    	  mrna_coord++;
    	  ostringstream oss1, oss2;
    	  oss1 << seq_name << "@" << cursor; //absolute coordinate on genomic sequence
    	  oss2 << seq_name << "@" << mrna_coord; //relative coordinate on mrna sequence
    	  string key1 = oss1.str(), key2 = oss2.str();
    	  map<string, s32>::const_iterator itPhase = protPhaseCoord.find(key1);
    	  if(itPhase != protPhaseCoord.end() && itPhase->second){
    		  refScore+=1;//alignment proteique exist ?
    		  modelPhaseCoord.insert(make_pair(key2, itPhase->second)); //careful: strand<0 coordinates in the map are stored in forward
    	  }
      }
    }
  }
  if(refScore) {
    _ref_score = refScore;
    _protPhaseCoord = modelPhaseCoord;
  }
  if(modelPhaseCoord != _protPhaseCoord) modelPhaseCoord = _protPhaseCoord;
  transform(mrna_seq.begin(), mrna_seq.end(), mrna_seq.begin(), (int(*)(int))toupper);
  /*reverse the sequence if the strand is - */
  if (strand<0) { 
    reverse(mrna_seq.begin(),mrna_seq.end());
    transform(mrna_seq.begin(), mrna_seq.end(), mrna_seq.begin(), complem);
  }
  /*define the start and stop codons*/
  vector<string> startC, stopC;
  startC.push_back("ATG");
  stopC.push_back("TAA");
  stopC.push_back("TAG");
  stopC.push_back("TGA");
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
    bool findStop = false;
    for ( itStop = stop.begin() ; itStop < stop.end(); itStop++ ) {
      s32 pos = *itStop;
      if((pos-i)%3 == 0) { findStop = true; break; }
    }
    if(!findStop) {
      pair<s32,s32> orf_begin_end = make_pair(i+1, (s32)mrna_seq.length());
      s32 len_begin_end = orf_begin_end.second - orf_begin_end.first + 1;
      if(_ref_score) {
    	  s32 score_begin_end = this->getPhaseScore(orf_begin_end, strand, mrna_seq.length(), modelPhaseCoord);
    	  if(score_begin_end > phaseScore) {
    		  phaseScore = score_begin_end;
    		  cds = orf_begin_end;
    		  cdslen = len_begin_end;
    		  _stop_found = 0;
    		  _start_found = 0;
 //   		  cout<<" _stop_found=0 et start_found =0 "<< cds.first << " cds.second "<<cds.second<<endl;
    	  }
      }
      else if(len_begin_end > cdslen) {
    	  cds = orf_begin_end;
    	  cdslen = len_begin_end;
    	  _stop_found = 0;
    	  _start_found = 0;
 //   	  cout<<" else len_begin_end > cdslen "<< cds.first << " cds.second "<<cds.second <<endl;
      }
    }
  }

  // find the longest M->* orf
  pair<s32,s32> orf_M_stop;
  s32 len_M_stop;
  if(_ref_score) {
//	  cout << "score_M_stop"<<endl;
    s32 score_M_stop;
    orf_M_stop = this->bestORF(start, stop, strand, mrna_seq.length(), score_M_stop,alreadyLookPos, modelPhaseCoord);
 //   cout << " start orf orf_M_stop " << orf_M_stop.first << " stop Orf " << orf_M_stop.second << " score "<< score_M_stop<<endl;
    len_M_stop = orf_M_stop.second - orf_M_stop.first + 1;
//    cout << "score_M_stop "<<score_M_stop<<" phaseScore "<<phaseScore<<endl;
    if(score_M_stop > phaseScore) {
      phaseScore = score_M_stop;
      cds = orf_M_stop;
      cdslen = len_M_stop;
      _stop_found = 1;
      _start_found = 1;
 //     cout << "_stop_found = 1 et _start_found = 1 "<<cds.first <<" " <<cds.second<<endl;
    }
  }
  
  else{
//	  cout<<"else score_M_stop"<<endl;
    orf_M_stop = this->longestORF(start, stop,alreadyLookPos);
    len_M_stop = orf_M_stop.second - orf_M_stop.first + 1;
 //   cout<<"len_M_stop "<<len_M_stop<<" cdslen "<<cdslen<<endl;
    if(len_M_stop > cdslen) {
      cds = orf_M_stop;
      cdslen = len_M_stop;
      _stop_found = 1;
      _start_found = 1;
 //    cout << "_stop_found = 1 et _start_found = 1 "<<cds.first <<" " <<cds.second<<endl;
    }  
  }
  
  // find the longest M->end orf
  vector<s32> allStop(stop);
  allStop.push_back(mrna_seq.length()-3);
  allStop.push_back(mrna_seq.length()-4);
  allStop.push_back(mrna_seq.length()-5);
  std::sort(allStop.begin(), allStop.end());
  for(vector<s32>::iterator itAllStop = allStop.begin();itAllStop!=allStop.end();++itAllStop){
//	  cout <<"itAllStop " << *itAllStop<<endl;
  }
  pair<s32,s32> orf_M_end;
  s32 len_M_end;

  if(_ref_score) {
//	  cout <<"score_M_end"<<endl;
    s32 score_M_end;
    orf_M_end = this->bestORF(start, allStop, strand, mrna_seq.length(), score_M_end,alreadyLookPos, modelPhaseCoord);
  //  cout << " start orf orf_M_end " << orf_M_end.first << " stop Orf " << orf_M_end.second << " score "<< score_M_end<<endl;
    len_M_end = orf_M_end.second - orf_M_end.first + 1;
//    cout <<"score_M_end "<<score_M_end<<" phaseScore "<< phaseScore<<endl;
    if(score_M_end > phaseScore) {
      phaseScore = score_M_end;
      cds = orf_M_end;
      cdslen = len_M_end;
      _start_found = 1;
      _stop_found = 0;
 //     cout << "_start_found = 1 et _stop_found = 0 "<<cds.first << " "<<cds.second<<endl;
    }
  }
  else {
//	  cout<<"else score_M_end"<<endl;
    orf_M_end = this->longestORF(start, allStop,alreadyLookPos);
    len_M_end = orf_M_end.second - orf_M_end.first + 1;
 //   cout<<"len_M_end "<<len_M_end<<" cdslen "<<cdslen <<endl;
    if(len_M_end > cdslen) {
      cds = orf_M_end;
      cdslen = len_M_end;
      _start_found = 1;
      _stop_found = 0;
//      cout << "_start_found = 1 et _stop_found = 0 "<<cds.first << " "<<cds.second<<endl;
    }

  }
  
  // find the longest begin->* on a different frame that the existing longest cds, if the begin->* segment is at least 64bp longer
  vector<s32> allStart(start);
  for(s32 i = 0; i < 3; i++) if((i-cds.first+1)%3 !=0) allStart.push_back(i);
  std::sort(allStart.begin(), allStart.end());
  pair<s32,s32> orf_begin_stop;
  s32 len_begin_stop; 

  if(_ref_score) {//XXX Is it possible to have that, I never find it ?
//	  cout <<"score_begin_stop"<<endl;
    s32 score_begin_stop;
    orf_begin_stop = this->bestORF(allStart, allStop, strand, mrna_seq.length(), score_begin_stop,alreadyLookPos, modelPhaseCoord);//XXX I replace stop by allStop
 //   cout << " start orf orf_begin_stop " << orf_begin_stop.first << " stop Orf " << orf_begin_stop.second << " score "<< score_begin_stop<<endl;
    len_begin_stop = orf_begin_stop.second - orf_begin_stop.first + 1;
//    cout <<"score_begin_stop "<<score_begin_stop<<" phaseScore "<<phaseScore<<endl;
    if(score_begin_stop > phaseScore){
      phaseScore = score_begin_stop;
      cds = orf_begin_stop;
      cdslen = len_begin_stop;
      _start_found = 0;
      _stop_found = 1;
 //     cout <<"_start_found = 0 et _stop_found = 1 "<<cds.first<<" "<<cds.second<<endl;
    }
  }
  
  else {
//	  cout<<"else score_begin_stop"<<endl;
    orf_begin_stop = this->longestORF(allStart, stop,alreadyLookPos);
    len_begin_stop = orf_begin_stop.second - orf_begin_stop.first + 1;
 //   cout<<"len_begin_stop "<<len_begin_stop<<" cdslen + 64 "<<cdslen + 64<<endl;
    if(len_begin_stop > cdslen + 64) {
      cds = orf_begin_stop;
      cdslen = len_begin_stop;
      _start_found = 0;
      _stop_found = 1;
 //     cout <<"_start_found = 0 et _stop_found = 1 "<<cds.first<<" "<<cds.second<<endl;
    }
  
  }

  // SECOND STEP IF NO STOP WAS FOUND, TAKE THE FIRST ONE FOUND IN THE SAME FRAME, AT A MAX DISTANCE OF _WINDOW
//  cout << "second step "<<cds.first << " cds.second "<<cds.second << endl;
  s32 frame;
  for (frame = 0; frame < 3; frame++){
    if ((cds.first-1-frame)%3==0) {break;}
  }
  
  string mrna_xlarge, mrna_3p, mrna_5p;
  string* whole_seq = (*_exons.begin())->getWholeSeq();
  s32 start_xlarge, size_xlarge; 
  s32 newStopPos = -1, newStartPos = -1;
  
  if (_stop_found == 0 && cds!=make_pair(0,0)) {
//	  cout << "_stop_found==0"<<endl;
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
    
//    cout << " just before if " << cds.first << " cds.second "<< cds.second << endl;
    if(newStopPos != -1) {
  //  	cout <<"new stop_found=1"<<endl;
      cds = make_pair(cds.first, newStopPos + 3);
      cdslen = cds.second - cds.first + 1;
      _stop_found = 1;
      _3pXL = newStopPos + 3 - (s32)mrna_seq.length();
  //    cout << " if newStopPos != -1 cds.first "<<cds.first << " cds.second "<< cds.second<< "_3pXL " << _3pXL << endl;
    }
    
  }
//  cout << "_3pXL after if " << _3pXL << endl;
  // THIRD STEP IF NO START WAS FOUND, TAKE THE CLOSEST ONE TO THE BEGINNING OF THE CDS IN THE LIMIT OF (-_WINDOW,+_WINDOW)
//  cout << "third step "<<cds.first << " cds.second " << cds.second<< endl;
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
    if (stop.empty()) {lastStop=0;} else {lastStop = stop.back()+3;}

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
  //  	  cout<<" if newStartPos != -1 "<<cds.first <<" "<<cds.second<<" _5pXL "<< _5pXL<<endl;
      }
      else {
    	  _5pXL = 0;
    	  cds = make_pair(newStartPos-(s32)mrna_5p.length()+1,cds.second);
 //   	  cout<<" else newStartPos == -1 "<< cds.first << " " << cds.second<<endl;
      }
      cdslen = cds.second - cds.first + 1;
//      cout <<"new start found=1 "<<cds.first<<endl;
    }
  }
//  cout << " cds select " << cds.first << " " << cds.second<<" 3pXL "<<_3pXL << endl;
  if(cdslen != 0) _cds = cds;
  else _cds = make_pair(0, 0);
  _cds_size = cdslen;

//  cout <<"phase score " << phaseScore<<endl;
  return phaseScore;
}

// to calculate the score of an ORF based on the proteic mapping phase information
s32 GeneModel::getPhaseScore(pair<s32, s32> orf, s32 strand, s32 mrna_len, const map<string, s32>& modelPhaseCoord){
  s32 score=0, rev_j=0;
  for(s32 j=orf.first;j<=orf.second;j++){
    ostringstream oss;
    if(strand>0) { oss << this->getSeqname() << "@" << j ; }
    if(strand<0) {
      rev_j = mrna_len - j + 1;
      oss << this->getSeqname() << "@" << rev_j;
    }
    string key = oss.str();
    map<string, s32>::const_iterator it = modelPhaseCoord.find(key);
    if(it != modelPhaseCoord.end() && it->second){
      s32 phase_diff;
      if(strand>0) phase_diff = (j - orf.first)%3 +1 - it->second;
      else {phase_diff = (j - orf.first)%3 +1 + it->second; }
      if(!phase_diff) score++;
    }
  }
  return score;
}

// to find the longest ORF between all combinaisons of start and stop
pair<s32,s32> GeneModel::longestORF(vector<s32> lStart, vector<s32> lStop,map<string,bool>&alreadyLookPos) {
	 std::clock_t c_start = std::clock();
  vector<s32>::iterator itStart, itStop;
  pair<s32,s32> orf;
  s32 max_len = 0;
  for ( itStart = lStart.begin() ; itStart < lStart.end(); itStart++ ) {
    for ( itStop = lStop.begin() ; itStop < lStop.end(); itStop++ ) {
      s32 pstart = *itStart, pstop = *itStop;
      s32 cds_len = pstop - pstart + 1 + 2;
//     ostringstream oss;
 //          	 oss <<pstart << "@" << pstop;
  //         	 string key = oss.str();
  //        	 map<string, bool>::iterator itAlreadylook = alreadyLookPos.find(key);
 //     cout<<"pstart " <<pstart<<" pstop "<<pstop <<"cds_len "<<cds_len<<" max_len "<<max_len<<endl;
      if(cds_len <=0 || cds_len%3 != 0) { continue; }
 //     if(itAlreadylook != alreadyLookPos.end()){
 //   	  itStop = lStop.end();
 //   	   continue;
 //     }

      else if(cds_len > max_len) {
		orf = make_pair(pstart + 1, pstop + 3);
		max_len = cds_len;
      }
  //    alreadyLookPos.insert( make_pair(key,true) );
      itStop = lStop.end();
 //     cout<<"itStop = lStop.end() "<<*itStop<<endl;
    }
  }
//  cout<<"longest orf " << max_len<<endl;
  std::clock_t c_end = std::clock();
//   cout << " time longestORF " <<  1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << endl;
  return orf;
}

// to find the ORF between all combinaisons of start and stop that has the best phase score
pair<s32,s32> GeneModel::bestORF(vector<s32> lStart, vector<s32> lStop, s32 strand, s32 mrna_len, s32& score,map<string,bool>& alreadyLookPos, const map<string, s32>& modelPhaseCoord ){
	 std::clock_t c_start = std::clock();

	 vector<s32>::iterator itStart, itStop;
	 pair<s32,s32> orf, tmp_orf;
	 s32 max_score=0;
	 for ( itStart = lStart.begin() ; itStart < lStart.end(); itStart++ ) {
		 for ( itStop = lStop.begin() ; itStop < lStop.end(); itStop++ ) {
			 s32 pstart = *itStart, pstop = *itStop;
			 s32 cds_len = pstop - pstart + 1 + 2;
			 if(cds_len <=0 || cds_len%3 != 0) { continue; }
			 tmp_orf = make_pair(pstart + 1, pstop + 3);

			 std::clock_t c_start_getPhaseScore = std::clock();
			 score = this->getPhaseScore(tmp_orf, strand, mrna_len, modelPhaseCoord);
			 std::clock_t c_end_getPhaseScore = std::clock();
			 if(PRINTTIME)	cout << " time getPhaseScore " <<  1000.0 * (c_end_getPhaseScore-c_start_getPhaseScore) / CLOCKS_PER_SEC << endl;
			 if(score > max_score) {
				 orf = tmp_orf;
				 max_score = score;
			 }
			 itStop = lStop.end();
		 }
	 }
	 score = max_score;
	 std::clock_t c_end = std::clock();
	 if(PRINTTIME) cout << " time bestORF " <<  1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << endl;
	 return orf;
}

// to recalculate the absolute coordinate of the ORF on the genomic sequence
void GeneModel::mapORF() {//XXX Why we need to recalculate ?
//	cout <<"mapOrf  _cds.first " << _cds.first << " _cds.second "<< _cds.second << endl;
//	cout << " nb exons "<< _exons.size()<<endl;
  s32 debOrf = _cds.first, finOrf = _cds.second;
  s32 posEp=0; // exonic sequence pos
  s32 posG=0; // genome position
  s32 orfgin=0, orfgout=0;
  if (_strand >= 0) { // if forward or unknown
    posEp=1;
    for(TSSRList::iterator it = _exons.begin() ; it != _exons.end(); it++) {
      s32 debExon = (*it)->start();
      s32 finExon = (*it)->end();
      if(it == --_exons.end()) {finExon = min((finExon + _window), (s32)(*it)->getWholeSeq()->length());}
      if(it == _exons.begin()) {debExon-=_5pXL;}
      for(posG = debExon ; posG <= finExon ; posG++){
    	  if (posEp == debOrf) orfgin=posG;
    	  if (posEp == finOrf) orfgout=posG;
 //   	  cout << " strand >=0 orfgin "<< orfgin << " " << orfgout << endl;
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
      if(rit == --_exons.rend()) {debExon -= _3pXL;} 
      if(rit == _exons.rbegin()) {finExon += _5pXL;}
 //     cout << " _3pXL " << _3pXL << " _5pXL "<< _5pXL << endl;
      for(posG = finExon ; posG >= debExon ; posG--){
    	  if (posEp == debOrf) orfgin=posG;
    	  if (posEp == finOrf) orfgout=posG;
  //  	  cout << " strand < 0 orfgin "<< orfgin << " " << orfgout << " posEp "<< posEp << endl;
    	  if (posEp > finOrf) break;
    	  posEp++;
      }
    }

      _cds = make_pair(orfgout,orfgin);
    if(!this->containCDS())
    {cerr << "error in mapORF _cds.first " << _cds.first << " _cds.second "<< _cds.second <<endl;//FIXME THAT APPEND SOMETIMES !!
 //   exit(1);
    }
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


void selectModel(map<s32,list<GeneModel> >& mapOfGene){
	 for(map<s32,list<GeneModel> >::iterator itMap = mapOfGene.begin(); itMap != mapOfGene.end() ; ++itMap){
		 if(itMap->second.size() == 1 ) {
			 continue;
		 }
		  cout << " we look at gene end " << itMap->first << endl;
		  selectModel(itMap->second);
	 }
}
void selectModel(list<GeneModel>& listOfGene){
	cout<<" selectModel "<<endl;
	cout << "listOfGene size before " << listOfGene.size() << endl;
	bool findCDS1 = false, findCDS2=false;
	TSSRList::iterator it,it2;
	list<GeneModel>::iterator itNext ;
	list<GeneModel> tmpList = listOfGene;
	for(list<GeneModel>::iterator itOfGene = listOfGene.begin(); itOfGene != --listOfGene.end();++itOfGene){
	//	list<GeneModel>::iterator itNext = itOfGene;
		cout << " itOfGene " << itOfGene->getCDS().first << " " << itOfGene->getCDS().second << " " << itOfGene->getStartFound() << " " << itOfGene->getStopFound() << endl;

		for(list<GeneModel>::iterator itNext = listOfGene.begin(); itNext != listOfGene.end();++itNext){

			if(itOfGene->getCDS().first != itNext->getCDS().first){
				cout <<" if(itOfGene->getCDS().first != itNext->getCDS().first) "<<endl;
				continue;
			}
			if(itOfGene->getCDS().second == itNext->getCDS().second){
				cout << "if(itOfGene->getCDS().second == itNext->getCDS().second) "<<endl;
				continue;
			}
/*
				if( itNext->getStopFound() == 1 &&itOfGene->getStopFound()==1){
					cout << " if( itNext->getStartFound() == 1 &&itOfGene->getStartFound()==1) "<<endl;
					continue;
				}
*/			cout << " itNext " << itNext->getCDS().first << " " << itNext->getCDS().second << " " << itNext->getStartFound() << " " << itNext->getStopFound() << endl;
	//		if(itOfGene->getStartFound() && itOfGene->getStartFound() && itNext->getStopFound() && itNext->getStopFound()){

		//		++itNext;
//				continue;
	//		}
	//		if((itOfGene->getStartFound() && ! itOfGene->getStartFound()) || (itNext->getStopFound() && !itNext->getStopFound())){
	//			cout << " we find a stop but no start ! " << endl;
	//			exit(1);
	//		}
			//	if(itOfGene->getCDS().first !=  itNext->getCDS().first){
			//	++itNext;
			//	continue;
		//	}
			list<pair<s32,s32> > exonCDS1, exonCDS2;
			s32 debCDS1 = itOfGene->getCDS().second, finCDS1 = itOfGene->getCDS().first;
			if (itOfGene->getCDS().first < itOfGene->getCDS().second) {
			   debCDS1 = itOfGene->getCDS().first;
			   finCDS1 = itOfGene->getCDS().second;
			 }
			cout << " debCds " << debCDS1 << " finCDS " << finCDS1 << endl;
			for ( it = itOfGene->getExons().begin(); it != itOfGene->getExons().end(); it++){
				if( ( (*it)->start() >= itOfGene->getCDS().first && (*it)->start() <=  itOfGene->getCDS().second ) || ( (*it)->end()>= itOfGene->getCDS().first &&  (*it)->end() <= itOfGene->getCDS().second ) ){
					s32 debExon1 = (*it)->start();
					s32 finExon1 = (*it)->end();
					if( finCDS1 < debExon1 || debCDS1 > finExon1 ) {
						continue;
					}
					else {
						s32 end = (finCDS1 < finExon1) ? finCDS1 : finExon1;
						s32 begin = (debCDS1 < debExon1) ? debExon1 : debCDS1;
						cout << " cds1 " << begin << " " << end << endl;
						exonCDS1.push_back(make_pair(begin,end));
					}
				}
			}
			s32 debCDS2 = itOfGene->getCDS().second, finCDS2 = itOfGene->getCDS().first;
			if (itOfGene->getCDS().first < itOfGene->getCDS().second) {
			   debCDS2 = itOfGene->getCDS().first;
			   finCDS2 = itOfGene->getCDS().second;
			 }
			for( it2= itNext->getExons().begin(); it2 != itNext->getExons().end(); it2++){
				if( ( (*it2)->start() >= itOfGene->getCDS().first && (*it2)->start() <=  itOfGene->getCDS().second ) || ( (*it2)->end()>= itOfGene->getCDS().first &&  (*it2)->end() <= itOfGene->getCDS().second ) ){
					s32 finExonNext = (*it2)->end();
					s32 debExonNext = (*it2)->start();
					if( finCDS2 < debExonNext || debCDS2 > finExonNext ) {
						continue;
					}
					else {
						s32 end2 = (finCDS2 < finExonNext) ? finCDS2 : finExonNext;
						s32 begin2 = (debCDS2 < debExonNext) ? debExonNext : debCDS1;
						cout << " cds2 " << begin2 << " " << end2 << endl;
						exonCDS2.push_back(make_pair(begin2,end2));
					}
			   }
			}
	/*		if(exonCDS1.size() != exonCDS2.size()){
				 ++itNext;
				continue;
			}
	*/		bool sameCDS = true;
			list<pair<s32,s32> >::iterator itExon1,itExon2;
			for(itExon1 = exonCDS1.begin(),itExon2 = exonCDS2.begin(); itExon1!= exonCDS1.end(), itExon2!= exonCDS2.end();++itExon1, ++itExon2){
				if(itExon1->first == itExon2->first && itExon1->second == itExon2->second){
					cout << " itExon1->first == itExon2->first && itExon1->second == itExon2->second ";
					cout << itExon1->first << " " << itExon2->first << " "<< itExon1->second  << " "<< itExon2->second<<endl;
				//	++itNext;
					continue;
				}
				else{
					sameCDS = false;
					cout << " same cds = false " << endl;
					cout << itExon1->first << " " << itExon2->first << " "<< itExon1->second  << " "<< itExon2->second<<endl;
					if(itExon1->second == itOfGene->getCDS().second ||itExon2->second == itNext->getCDS().second ){
						cout << " we were at the end of the cds !! " <<endl;
						exit(1);
					}
					else
						break;
					//exit(1);
				}
			}
		/*	if(sameCDS == true){
				if(itOfGene->getModelSize() > itNext->getModelSize()){
					itNext = listOfGene.erase(itNext);
				}
				else{
					itOfGene = listOfGene.erase(itOfGene);
					++itNext;
				}
			}
			else{
				 ++itNext;
			}
	*/
			//++itNext;
		}
		//++itOfGene;
	}
	cout << "listOfGene size after " << listOfGene.size() << endl;
}



void selectModel( map<string,GeneModel> mapGeneModel){
	for(map<string,GeneModel>::iterator itMap = mapGeneModel.begin(); itMap != mapGeneModel.end();++itMap){
		for(map<string,GeneModel>::iterator itMap2 = mapGeneModel.begin(); itMap2 != mapGeneModel.end();++itMap2){
			if(itMap == itMap2 )continue;
			if(itMap->second.getStartFound() && itMap2->second.getStartFound() && itMap->second.getStopFound() && itMap2->second.getStopFound())continue;
			if((itMap->second.getCDS().first <= itMap2->second.getCDS().first && itMap->second.getCDS().second >= itMap2->second.getCDS().second) || (itMap->second.getCDS().first >= itMap2->second.getCDS().first && itMap->second.getCDS().second <= itMap2->second.getCDS().second)){
				cerr << " overlap CDS "<<itMap->second.getCDS().first << " "<<itMap->second.getCDS().second <<"start found " << itMap->second.getStartFound() << " end " << itMap->second.getStopFound()<< "\n" << itMap2->second.getCDS().first <<" "<< itMap2->second.getCDS().second<<"start found " << itMap2->second.getStartFound() << " end " << itMap2->second.getStopFound() <<endl;
				//exit(1);
			}
		}
	}
	cout << " end select Orf " << endl;
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
	    if (it == _exons.begin()) {debExon -= _5pXL;}
	    if (it == --_exons.end()) {finExon += _3pXL;}
	    }
	    else if(_strand < 0) {
	      if (it == --_exons.end()) {finExon += _5pXL;}
	      if (it == _exons.begin()) {debExon -= _3pXL;}
	    }
	    if( finCDS < debExon || debCDS > finExon ) {
	    	continue;
	    }
	    else {
	    	s32 end = (finCDS < finExon) ? finCDS : finExon;
	    	s32 begin = (debCDS < debExon) ? debExon : debCDS;
			oss <<"@"<< begin << "@" << end ;
	    }
	  }
	string key = oss.str();
	map<string,GeneModel>::iterator itMap = mapGeneModel.find(key);
	if( itMap== mapGeneModel.end()){
		pair<string,GeneModel> tmpPair = make_pair(key,*this);
		mapGeneModel.insert(tmpPair);
	}
	else{
		if(this->_model_size >= itMap->second._model_size){
			itMap->second= *this;
		}
		//else do nothing
	}

}

bool GeneModel::overlapOrf(pair<s32,s32> cds2){

	return ((this->_cds.first <= cds2.first && this->_cds.second >= cds2.first) || (this->_cds.first >= cds2.first && this->_cds.first <= cds2.second));
}
