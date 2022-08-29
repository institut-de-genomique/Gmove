/*******************************************************************************
+
+  GffRecordList.cpp
+
+  Copyright (c) 2017 Genoscope, CEA, CNS, Evry, France
+  Author : Marion Dubarry, jmaury@genoscope.cns.fr
+
 *******************************************************************************/

#include "GffRecordList.h"

using namespace std;

GffRecordList::GffRecordList(char * fileName,map<string, string>fastaDB, s8 datatype){
  fstream file;
  string line;
  map<string,bool> mapAttribute;
  string previousAttribute = "";
  
  checkFormatGff(fileName);
  file.open(fileName, ios::in);
  
  while(getline( file, line )){
    if (line[0]=='#') continue;
    
    GffRecord* record = new GffRecord(line, datatype);
    string scaff = record->getSeqName();
    map<string,string>::iterator itFasta = fastaDB.find(scaff);
    if(itFasta == fastaDB.end()) { delete(record); continue; }
    
    if(datatype == RNA_DATA || datatype == PROTEIN_DATA){
      if(record->getType() != "exon" && record->getType() != "HSP") { 
	delete(record); 
	continue; 
      }
      if(record->getAttribute() != previousAttribute && previousAttribute != "")
	checkUniqId(mapAttribute, previousAttribute);
      _records.push_back(record);
    }
    else if(datatype == ANNOT_DATA || datatype == ABINITIO_DATA){      
      if(record->getType() != "CDS" && record->getType() != "UTR" ) { 
	delete(record); 
	continue; 
      }
      if(record->getAttribute() != previousAttribute && previousAttribute != "")
	checkUniqId(mapAttribute, previousAttribute);
      _records.push_back(record);
    }
    previousAttribute = record->getAttribute();
  }
  
  _records.sort(sortRecordGff);
  file.close();
}

void checkUniqId(map<string,bool>& mapAttribute, string currentAttribute){
  map<string,bool>::iterator itMapAttribute = mapAttribute.find(currentAttribute);
  if( itMapAttribute != mapAttribute.end()){
    cerr << "Error : id is not uniq : " << currentAttribute <<endl;
    exit(1);
  }
  else{
    pair<string,bool> tmpPair = make_pair(currentAttribute,true);
    mapAttribute.insert(tmpPair);
  }
}

void GffRecordList::copyGffRecordList(GffRecordList gffRecord){
  GffRecordL* records = gffRecord.getRecords();
  for(GffRecordL::iterator itG = records->begin() ; itG != records->end();++itG){
    _records.push_back(*itG);
  }
}

void GffRecordList::printGffRecord(){
  for(GffRecordL::iterator itG = _records.begin() ; itG != _records.end();++itG){
    cout << *(*itG) <<   endl;
  }
}

void GffRecordList::checkFormatGff(char * fileName ){
  fstream file;
  file.open(fileName, ios::in);
  string line;
  while(getline( file, line )){
    if (line[0]=='#')//comment line
      continue;
    std::string delimiter = "\t";
    s32 cmpt = 0;
    size_t pos = 0;
    string tokens;
    
    while ((pos = line.find(delimiter)) != string::npos) {
      tokens = line.substr(0, pos);
      ++cmpt;
      line.erase(0, pos + delimiter.length());
    }
    if(!line.empty())
      ++cmpt;
    
    if(cmpt != 9){
      cerr << "error Gff file "<< fileName <<" hasn't 9 columns " << endl;
      exit(1);
    }
  }
  file.close();
}

/*
  We need to keep in mind where come from the exon, we cannot 
  build all combinaison from this data
*/
void GffRecordList::loadGff(map<string, map<s32,s32> >& cds){
  s32 phase;
  map<string,list<pair<s32,s32> > > posCds;
  map<string,s32> selectedCds;
  map<string,bool> idTagRecord;
  string previous_exon_id,tag,previousAttribute, scaff, currentAttribute;
  _records.sort(sortRecordGff);
  for(GffRecordL::iterator itR = _records.begin(); itR != _records.end();++itR){
    GffRecord* record = *itR;
    currentAttribute = record->getAttribute() ;
    if(currentAttribute != previousAttribute && previousAttribute != "") {
       cdsPhase(phase, posCds[previousAttribute], scaff, cds);
    }
    scaff = record->getSeqName();
    pair < int, int > pairPos = make_pair(record->getStart(),record->getEnd());
    posCds[currentAttribute].push_back(pairPos);
    if(record->getStrand()=="+") phase = 1;
    else phase = -1;
    // @@@@ MODIF-JM-JULY
    //if(currentAttribute != previousAttribute && previousAttribute != "")
    //	cdsPhase(phase,posCds[previousAttribute],scaff,cds);
    // @@@@ MODIF-JM-JULY
    previousAttribute = currentAttribute;
  }
}

// use it for reannotation
// We need to keep in mind where come from the exon, we cannot 
// build all combinaison from this data
void GffRecordList::loadAnnotation( map<string, map<s32,s32> >& cds, bool getCds){
  map< string,  list < pair < s32,s32> > > posCds;
  s32 phase; //phase = 1, 2 or 3
  string previousAttribute;
  s32 previousEnd = -1 , previousStart = -1 ;
  for(GffRecordL::iterator itR = _records.begin(); itR != _records.end();++itR){
    GffRecord* record = *itR;
    string scaff = record->getSeqName();
    string currentAttribute = record->getAttribute();
    string type = record->getType();
    
    if(currentAttribute != previousAttribute && previousAttribute != "" && getCds){
      if(record->getStrand()=="+")phase = 1;
      else phase = -1;
      cdsPhase(phase, posCds[previousAttribute], scaff, cds);
    }
    if(type == "CDS" && getCds){
      pair < int, int > pairPos = make_pair(record->getStart(),record->getEnd());
      posCds[currentAttribute].push_back(pairPos);
    }
    if(previousEnd +1 == record->getStart()){
      record->setStart(previousStart);
      
      itR = _records.erase(--itR);
    }
    record->setType("exon");
    //update the last one
    previousAttribute = currentAttribute;
    previousEnd = record->getEnd();
    previousStart = record->getStart();
  }
}

void GffRecordList::intron() {
  string previousId;
  s32 previousEnd = -1, currentStart = -1;
  map< string, s32 > introns;
  s32 color = 0;
  s32 currentColor = 0;
  ostringstream tmpKey;
  GffRecordL recordsOneTranscrit;

  for(GffRecordL::iterator itRecord = _records.begin() ; itRecord != _records.end();++itRecord) {
    if((*itRecord)->getAttribute() == previousId) {
      currentStart = (*itRecord)->getStart();
      tmpKey << previousEnd << "@" << currentStart<<"@";
      recordsOneTranscrit.push_back(*itRecord);
      (*itRecord)->setColor(currentColor);
    }
    else {
      string key = tmpKey.str();
      map<string,s32>::iterator itIntrons = introns.find(key);
      if(itIntrons == introns.end()) {
	color++;
	currentColor = color;
	for(GffRecordL::iterator itPrevious = recordsOneTranscrit.begin() ; itPrevious != recordsOneTranscrit.end();++itPrevious)
	  (*itPrevious)->setColor(currentColor);
	introns.insert(make_pair(key,currentColor));
      }
      else {
	for(GffRecordL::iterator itPrevious = recordsOneTranscrit.begin() ; itPrevious != recordsOneTranscrit.end();++itPrevious)
	  (*itPrevious)->setColor(itIntrons->second);
	currentColor = itIntrons->second;
      }
      tmpKey.str("");
      recordsOneTranscrit.clear();
      tmpKey<<(*itRecord)->getSeqName() <<"@";
      recordsOneTranscrit.push_back(*itRecord);
    }
    previousId = (*itRecord)->getAttribute();
    previousEnd = (*itRecord)->getEnd();
  }
  currentColor++;
  for(GffRecordL::iterator itPrevious = recordsOneTranscrit.begin() ; itPrevious != recordsOneTranscrit.end();++itPrevious)
    (*itPrevious)->setColor(currentColor);
}

void cdsPhase(s32& ph, list< pair< s32,s32> > pos, string name, map<string,map<s32,s32> >& cds) {
  s32 phase = ph;
  s32 len = 0;
  for(list<pair <s32,s32 > >::reverse_iterator itPos = pos.rbegin(); itPos != pos.rend(); ++itPos) {
      s32 start = (*itPos).first;
      s32 end = (*itPos).second;
      len += end - start + 1;
  }
  if(len%3 == 0) {
    if(ph == 1) {
      for(list<pair <s32,s32 > >::iterator itPos = pos.begin(); itPos != pos.end(); ++itPos){
	s32 start = itPos->first;
	s32 end = itPos->second;
	for(int i=start ; i <= end ; ++i){
	  ostringstream oss;
	  oss << name << "@" << i;
	  string key = oss.str();
	  if( cds.find(key) == cds.end()) cds[key].insert(make_pair(phase,1));
	  else {
	    map<s32,s32>::iterator itCds = cds[key].find(phase);
	    if( itCds != cds[key].end()) itCds->second =  itCds->second + 1;
	    else cds[key].insert(make_pair(phase,1));
	  }
	  ++phase;
	  if(phase > 3) phase = 1;
	}
      }
    } 
    else {
      for(list<pair <s32,s32 > >::reverse_iterator itPos = pos.rbegin(); itPos != pos.rend(); ++itPos) {
	s32 start = (*itPos).first;
	s32 end = (*itPos).second;
	for(int i = end ; i >= start ; --i) {
	  ostringstream oss;
	  oss << name << "@" << i;
	  string key = oss.str();
	  if( cds.find(key) == cds.end()) cds[key].insert(make_pair(phase,1));
	  else {
	    map<s32,s32>::iterator itCds = cds[key].find(phase);
	    if( itCds != cds[key].end()) itCds->second =  itCds->second + 1;
	    else cds[key].insert(make_pair(phase,1));
	  }
	  --phase;
	  if(phase < -3) phase = -1;
	}
      }
    }
  }
}

void GffRecordList::cleanMono(){
  map<string, GffRecordL> mapMono = extractMono();
  for(map<string,GffRecordL >::iterator itMap = mapMono.begin(); itMap != mapMono.end();++itMap){    
    GffRecordL cleanRecord = itMap->second;
    fusionMonoExons(cleanRecord);
    for(GffRecordL::iterator itRecords = cleanRecord.begin(); itRecords != cleanRecord.end();++itRecords)
      _records.push_back(*itRecords);
  }
}

map<string, GffRecordL> GffRecordList::extractMono(){
  GffRecordL listMono;
  s32 nbexon = 0;
  string previousId;
  for (GffRecordL::iterator itRecord =_records.begin() ;itRecord != _records.end(); ) {
    if(previousId.empty()) {
      previousId = (*itRecord)->getAttribute();
      ++itRecord;
      continue;
    }
    else if(previousId == (*itRecord)->getAttribute()) {
      ++nbexon;
    }
    else {
      if(nbexon == 0) {
	GffRecord* tmpRecord = *--itRecord;
	itRecord = _records.erase(itRecord);
	listMono.push_back(tmpRecord);
      }
      nbexon = 0;
    }
    previousId = (*itRecord)->getAttribute();
    ++itRecord;
  }
  if(nbexon==0){
    GffRecord* tmpRecord = _records.back();
    GffRecordL::iterator itRecord = _records.end();
    _records.erase(--itRecord);
    listMono.push_back(tmpRecord);
  }
  map<string,GffRecordL> mapRecord;
  cerr << "There are " << listMono.size() << " single-exon genes "<< endl;
  insertMap(mapRecord,listMono);
  return mapRecord;
}

void GffRecordList::insertMap(map<string,GffRecordL>& mapRecord, GffRecordL listMono){
  for( GffRecordL::iterator itList = listMono.begin() ; itList != listMono.end(); ++itList){
    map<string,GffRecordL>::iterator itMapRecord;
    string scaff = (*itList)->getSeqName();
    itMapRecord = mapRecord.find(scaff);
    if(itMapRecord == mapRecord.end()){
      GffRecordL tmpList;
      tmpList.push_back(*itList);
      pair< string,GffRecordL > tmpPair;
      tmpPair = make_pair(scaff,tmpList);
      mapRecord.insert(tmpPair);
    }
    else
      itMapRecord->second.push_back(*itList);
  }
}

void GffRecordList::fusionMonoExons(GffRecordL& listMono){
  //Same idea as /env/ig/soft/rdbioseq/annotation-snapshot/linux-noarch/bin/loci
  listMono.sort(sortMonoGff);
  //--listMono.end() : we stop one element before the end, like size()-1
  for(GffRecordL::iterator itListMono = listMono.begin(); itListMono != --listMono.end() ; ) {
    GffRecordL::iterator itNext = itListMono;
    ++itNext;
    if((*itListMono)->getSeqName() == (*itNext)->getSeqName() && (*itListMono)->getStrand() == (*itNext)->getStrand()){ 
      if( ((*itNext)->getEnd() >= (*itListMono)->getStart() && (*itNext)->getEnd() <= (*itListMono)->getEnd())
	  || ((*itListMono)->getEnd() >= (*itNext)->getStart() && (*itListMono)->getEnd() <= (*itNext)->getEnd())){
	(*itNext)->setStart(min((*itListMono)->getStart(),(*itNext)->getStart()));
	
	(*itNext)->setEnd(max((*itListMono)->getEnd(),(*itNext)->getEnd()));
	itListMono = listMono.erase(itListMono);
      }
      else
	++itListMono;
    }
    else
      ++itListMono;
  }
}


bool sortMonoGff (GffRecord* record1 , GffRecord*  record2){
  return (record1->getStart() <= record2->getStart());
}

bool sortRecordGff (GffRecord*  record1 , GffRecord* record2){
  return ( (record1->getSeqName() < record2->getSeqName()) ||
	   ( record1->getSeqName() == record2->getSeqName() && record1->getAttribute() < record2->getAttribute()) ||
	   ( record1->getSeqName() == record2->getSeqName() && record1->getAttribute() == record2->getAttribute()
	     && record1->getStart() < record2->getStart() )
	   
	   );
}

