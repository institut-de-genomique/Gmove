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


GffRecordList::GffRecordList(char * fileName,map<string,string>fastaDB, string typeData){
	  fstream file;
	  string line;
	  map<string,bool> mapAttribute;
	  string previousAttribute = "";

	  checkFormatGff(fileName);
	  file.open(fileName, ios::in);

	  while(getline( file, line )){
		  if (line[0]=='#')//comment line
		  			continue;
		  GffRecord* record = new GffRecord(line);
		  string scaff = record->getSeqName();
		  map<string,string>::iterator itFasta = fastaDB.find(scaff);
		  if(itFasta==fastaDB.end())continue;

		  if(typeData == "rna" || typeData =="prot"){
			  if(record->getType() != "exon" && record->getType() != "HSP") continue;
			  if(record->getAttribute() != previousAttribute && previousAttribute != "")
						checkUniqId(mapAttribute, previousAttribute);
		//	  cout << "record " << *record << endl;
			  _records.push_back(record);
		  }
		  else if(typeData == "annot" || typeData == "abinitio"){

			  if(record->getType() != "CDS" && record->getType() != "UTR" ) continue;
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


void GffRecordList::loadGff(map<string, map<s32,s32> >& cds){
	/*
	We need to keep in mind where come from the exon, we cannot build all combinaison from this data
	*/
	 s32 phase;
	 map<string,list<pair<s32,s32> > > posCds;
	 map<string,s32> selectedCds;
	 map<string,bool> idTagRecord;
	 string previous_exon_id,tag,previousAttribute;
	 _records.sort(sortRecordGff);
	 for(GffRecordL::iterator itR = _records.begin(); itR != _records.end();++itR){
		GffRecord* record = *itR;
		string scaff = record->getSeqName();
		string currentAttribute = record->getAttribute() ;
		pair < int, int > pairPos = make_pair(record->getStart(),record->getEnd());
		posCds[currentAttribute].push_back(pairPos);
		if(record->getStrand()=="+")phase = 1;
		else phase = -1;
		if(currentAttribute != previousAttribute && previousAttribute != "")
			cdsPhase(phase,posCds[previousAttribute],scaff,cds);
		previousAttribute = currentAttribute;
	}
}


void GffRecordList::loadAnnotation( map<string, map<s32,s32> >& cds, bool getCds){
	// use it for reannotation
//	We need to keep in mind where come from the exon, we cannot build all combinaison from this data
	map< string,  list < pair < s32,s32> > > posCds;
	s32 phase;//phase = 1,2 or 3
	string previousAttribute;
	s32 previousEnd = 0 , previousStart = 0 ;
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
		//	GffRecord* record = *itR;
		}
		record->setType("exon");
	 //update the last one
	 previousAttribute = currentAttribute;
	 previousEnd = record->getEnd();
	 previousStart = record->getStart();
	}
}


void GffRecordList::intron(){
	string previousId;
	s32 previousStart, previousEnd,currentStart;
	map< string, s32 > introns;
	s32 color =0;
	s32 currentColor = 0;
	ostringstream tmpKey ;
	GffRecordL recordsOneTranscrit;
	for(GffRecordL::iterator itRecord = _records.begin() ; itRecord != _records.end();++itRecord){
		if((*itRecord)->getAttribute() == previousId){
			currentStart = (*itRecord)->getStart();
			tmpKey << previousEnd << "@" << currentStart<<"@";
			recordsOneTranscrit.push_back(*itRecord);
			(*itRecord)->setColor(currentColor);
		}
		else{
			string key = tmpKey.str();
	//		cout << "not the same id " << key << endl;
			map<string,s32>::iterator itIntrons = introns.find(key);
			if(itIntrons == introns.end()){
				color++;
								currentColor = color;
		//		cout <<"new transcrit never found "<<endl;
				for(GffRecordL::iterator itPrevious = recordsOneTranscrit.begin() ; itPrevious != recordsOneTranscrit.end();++itPrevious){
					(*itPrevious)->setColor(currentColor);
		//			cout << "color " << (*itPrevious)->getColor()<< " " << **itPrevious<< endl;
				}
		//		cout << "increment color "<< key << endl;
				introns.insert(make_pair(key,currentColor));

			//	(*itRecord)->setColor(color);


			}
			else{
		//		cout << "else transrti already known " << endl;
				for(GffRecordL::iterator itPrevious = recordsOneTranscrit.begin() ; itPrevious != recordsOneTranscrit.end();++itPrevious){
									(*itPrevious)->setColor(itIntrons->second);
		//							cout << "color " << (*itPrevious)->getColor()<< " " << **itPrevious<< endl;
								}
			//	(*itRecord)->setColor(itIntrons->second);
				currentColor = itIntrons->second;
			}
			tmpKey.str("");
			recordsOneTranscrit.clear();
			tmpKey<<(*itRecord)->getSeqName() <<"@";
			recordsOneTranscrit.push_back(*itRecord);
		}
		previousId = (*itRecord)->getAttribute();
		previousStart = (*itRecord)->getStart();
		previousEnd = (*itRecord)->getEnd();
		}
	currentColor++;
	for(GffRecordL::iterator itPrevious = recordsOneTranscrit.begin() ; itPrevious != recordsOneTranscrit.end();++itPrevious){
					(*itPrevious)->setColor(currentColor);
		//			cout << "color " << (*itPrevious)->getColor()<< " " << **itPrevious<< endl;
				}
//	for(GffRecordL::iterator itRecord = _records.begin() ; itRecord != _records.end();++itRecord){
//		cout << "color " << (*itRecord)->getColor()<< " " << **itRecord << endl;
//	}
}
void cdsPhase(s32& phase, list< pair< s32,s32> > pos, string name, map<string,map<s32,s32> >& cds){
	if(phase ==1){
		for(list<pair <s32,s32 > >::iterator itPos = pos.begin(); itPos != pos.end(); ++itPos){
			s32 start = itPos->first;
			s32 end = itPos->second;
			for(int i =start ; i <= end ; ++i){
				 ostringstream oss;
				 oss << name << "@" << i;
				 string key = oss.str();
				 if( cds.find(key) == cds.end()){
					 cds[key].insert(make_pair(phase,1));
				 }
				 else{
					 map<s32,s32>::iterator itCds = cds[key].find(phase);
					 if( itCds != cds[key].end())
						 itCds->second =  itCds->second + 1;
					 else
						 cds[key].insert(make_pair(phase,1));
				 }
				 ++phase;
				 if(phase > 3) phase = 1 ;
			 }
		}
	}
	else{
		for(list<pair <s32,s32 > >::reverse_iterator itPos = pos.rbegin(); itPos != pos.rend(); ++itPos){
				s32 start = (*itPos).first;
				s32 end = (*itPos).second;
			for(int i = end ; i > start ; --i){
				 ostringstream oss;
				 oss << name << "@" << i;
				 string key = oss.str();
				 if( cds.find(key) == cds.end())
					 cds[key].insert(make_pair(phase,1));
				 else{
					 map<s32,s32>::iterator itCds = cds[key].find(phase);
					 if( itCds != cds[key].end())
						 itCds->second =  itCds->second + 1;
					 else
						 cds[key].insert(make_pair(phase,1));
				 }
				--phase;
				if(phase < -3) phase = -1 ;
			 }
		}
	}
}

void GffRecordList::cleanMono(){
	cout << "clean mono " << endl;
	map<string, GffRecordL> mapMono = extractMono();
	for(map<string,GffRecordL >::iterator itMap = mapMono.begin(); itMap != mapMono.end();++itMap){

		GffRecordL cleanRecord = itMap->second;
		fusionMonoExons(cleanRecord);
		for(GffRecordL::iterator itRecords = cleanRecord.begin(); itRecords != cleanRecord.end();++itRecords){
			_records.push_back(*itRecords);
		}
	}
}
map<string, GffRecordL> GffRecordList::extractMono(){
	GffRecordL listMono;
	s32 nbexon = 0;
	string previousId;
	for (GffRecordL::iterator itRecord =_records.begin() ;itRecord != _records.end();){
	//	cout << "all nodes " << **itRecord << endl;
		if(previousId.empty()){
			previousId = (*itRecord)->getAttribute();
			++itRecord;
		//	++nbexon;
			continue;
		}
		else if(previousId == (*itRecord)->getAttribute()){
			++nbexon;
		}
		else{
			if(nbexon == 0){
				GffRecord* tmpRecord = *--itRecord;
				itRecord = _records.erase(itRecord);
				listMono.push_back(tmpRecord);
	//			cout << "mono " << *tmpRecord << endl;
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
//		cout << "mono " << *tmpRecord << endl;
	}
	map<string,GffRecordL> mapRecord;
	cout << "There are " << listMono.size() << " mono exonique genes "<< endl;
	insertMap(mapRecord,listMono);
	return mapRecord;
}
// map<string,list<GffRecord> > mapMono;
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
	for(GffRecordL::iterator itListMono = listMono.begin(); itListMono != --listMono.end() ;){//--listMono.end() : we stop one element before the end, like size()-1
		GffRecordL::iterator itNext = itListMono;
		++itNext;
		if((*itListMono)->getSeqName() == (*itNext)->getSeqName() && (*itListMono)->getStrand() == (*itNext)->getStrand()){ //(itListMono->strand == itNext->strand || itListMono->strand=='.' || itNext->strand == '.')){//&& itListMono->strand == itNext->strand){
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

