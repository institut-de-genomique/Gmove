/*
 * GeneModelList.cpp
 *
 *  Created on: 28 oct. 2016
 *      Author: mdubarry
 */

#include "GeneModelList.h"


void GeneModelList::deleteIncludedModel(){
  TSSRList::iterator it,it2;
  list<GeneModel>::iterator itNext ;
  list<GeneModel> tmpList = _models;
  s32 cdsSize1,cdsSize2;
  for(list<GeneModel>::iterator itOfGene = _models.begin(); itOfGene != --_models.end();++itOfGene){
    if(_models.size() == 1 )
      break;
    if(itOfGene->getToDelete() == true)
      continue;
    bool LOOP = true;
    for(list<GeneModel>::iterator itNext = itOfGene; itNext != _models.end() && LOOP ;++itNext){
      if(itNext->getToDelete() == true)
	continue;
      if(itOfGene == itNext)
	continue;
      
      if(itOfGene->getCDS().first == itNext->getCDS().first && itOfGene->getCDS().second == itNext->getCDS().second)
	//there is a first selection on same selection, if model are still in the list that means internal exons are not the same
	continue;
      
      
      list<pair<s32,s32> > exonCDS1, exonCDS2;
      exonCDS1 = itOfGene->getExonsCds();
      s32 finCDS1 = itOfGene->getCDS().first;
      if (itOfGene->getCDS().first < itOfGene->getCDS().second)
	finCDS1 = itOfGene->getCDS().second;
      cdsSize1 = itOfGene->getCdsSize();
      
      if(exonCDS1.size()==0)
	continue;
      exonCDS2 = itNext->getExonsCds();
      s32 finCDS2 = itNext->getCDS().first;
      if (itNext->getCDS().first < itNext->getCDS().second)
	finCDS2 = itNext->getCDS().second;
      cdsSize2 = itNext->getCdsSize();
      
      if(exonCDS2.size()==0)
	continue;
      list<pair<s32,s32> >::iterator itExon1,itExon2;
      for(itExon1 = exonCDS1.begin(),itExon2 = exonCDS2.begin(); itExon1!= exonCDS1.end(), itExon2!= exonCDS2.end() && LOOP;++itExon1, ++itExon2){
	if(itExon1== exonCDS1.end() || itExon2== exonCDS2.end())
	  break;
	if(itExon1->first < itExon2->first && itExon1->second < itExon2->first ){
	  --itExon2;
	  continue;
	}
	if(itExon1->first > itExon2->second && itExon1->second > itExon2->second){
	  --itExon1;
	  continue;
	}
	if((itExon1->first == itExon2->first && itExon1->second == itExon2->second )
	   || (itExon1 == exonCDS1.begin() && itExon1->first >= itExon2->first && itExon1->second == itExon2->second)
	   || (itExon2 == exonCDS2.begin() && itExon2->first >= itExon1->first && itExon2->second == itExon1->second)){
	  if(itExon1->second+ itOfGene->get3pXL() == finCDS1 && cdsSize1 <= cdsSize2){
	    itOfGene->setToDelete(true);
	    LOOP = false;
	    break;
	  }
	  else if (itExon2->second + itNext->get3pXL() == finCDS2  && cdsSize2 <= cdsSize1){
	    itNext->setToDelete(true);
	    break;
	  }
	  continue;
	}
	else{
	  if(itExon1->second>= finCDS1 && cdsSize1 <= cdsSize2 && itExon1->second < itExon2->second){ 
	    //On n'a plus la meme cds mais on verifie si on n'a pas atteint la fin de l'une d'elles
	    itOfGene->setToDelete(true);
	    LOOP = false;
	    break;
	  }
	  else if (itExon2->second>= finCDS2 && cdsSize2 <= cdsSize1 && itExon2->second < itExon1->second){
	    itNext->setToDelete(true);
	    break;
	  }
	  else
	    break;
	}
      }
    }
  }
}

void GeneModelList::longestCDS(){
  GeneModel* longestModel = NULL;
  for(list<GeneModel>::iterator itModel = _models.begin(); itModel != _models.end(); ++itModel){
    if(itModel== _models.begin()){
      longestModel = &*itModel;
      longestModel->setToDelete(false);
      continue;
    }
    if(itModel->getCdsSize() >= longestModel->getCdsSize()
       && itModel->getModelSize() +itModel->get3pXL() +itModel->get5pXL() >= longestModel->getModelSize() +longestModel->get3pXL() +longestModel->get5pXL()
       && !itModel->getToDelete()){
      if(itModel->getCdsSize() == longestModel->getCdsSize() && itModel->getModelSizeonGeno() > longestModel->getModelSizeonGeno()) {
	itModel->setToDelete(true); // Delete gene with too big intron (case of tandem duplicated gene)
      }
      else {
	if(itModel != _models.begin()) {
	  longestModel->setToDelete(true); //FIXME probleme in this if ?
	}
	longestModel = &*itModel;
	longestModel->setToDelete(false);
      }
    }
    else {
      itModel->setToDelete(true);
    }
  }
  this->catchModelOnUTR(longestModel);
}

void GeneModelList::bestScore(){
  GeneModel* longestModel = NULL;
  for(list<GeneModel>::iterator itModel = _models.begin(); itModel != _models.end(); ++itModel){
    if(itModel== _models.begin()){
      longestModel = &*itModel;
      longestModel->setToDelete(false);
      continue;
    }
    if(itModel->score() >= longestModel->score() && !itModel->getToDelete()){
      longestModel->setToDelete(true);
      longestModel = &*itModel;
      longestModel->setToDelete(false);
    }
    else itModel->setToDelete(true);
  }
  this->catchModelOnUTR(longestModel);
}


void GeneModelList::catchModelOnUTR(GeneModel * longestModel){// catch model taht are on utr from longest orf in cc
  for(list<GeneModel>::iterator itModel = _models.begin(); itModel != _models.end(); ++itModel){
    if(! longestModel->overlapOrf(itModel->getCDS())){
      itModel->setToDelete(false);
    }
  }
}

void GeneModelList::insertModels(GeneModelList smallList){
  std::copy(smallList._models.begin(), smallList._models.end(), std::back_inserter(_models));
}

void GeneModelList::includedMono(){ // a mono CDs overlap a pluri and the mono is < 300
  for(GeneModelL::iterator itModel = _models.begin(); itModel != _models.end(); ++itModel){
    if(itModel->getToDelete())
      continue;
    for(list<GeneModel>::iterator itNext = _models.begin(); itNext != _models.end(); ++itNext){
      if(itNext->getToDelete())
	continue;
      if(itModel == itNext)
	continue;
      if(itNext->cdsIsMono()){
	if(itModel->cdsIsMono() && itModel->overlapOrf(itNext->getCDS())){
	  if(itModel->getCdsSize() >= itNext->getCdsSize())
	    itNext->setToDelete(true);
	  else
	    itModel->setToDelete(true);
	}
	else if(itModel->overlapModel(*itNext) && itNext->getCdsSize() < 300)
	  itNext->setToDelete(true);
      }
    }
  }
}


void GeneModelList::clusterLocation(){ //Cluster location sur les CDS !!
  s32 numCluster = 1;
  map<s32,list<GeneModel*> > mapClusterModel;
  map<s32,list<GeneModel*> >::iterator itMapClusterModel;
  for(list<GeneModel>::iterator itModel1 = _models.begin(); itModel1 != _models.end() ; ++itModel1){
    if(itModel1->getToDelete())
      continue;
    for(list<GeneModel>::iterator itModel2 =  _models.begin(); itModel2 != _models.end() ;++itModel2){
      if(itModel2->getToDelete())
	continue;
      if(itModel1 == itModel2)
	continue;
      if(itModel1->overlapOrf(*itModel2) ){//&& itModel1->getStrand() == itModel2->getStrand()){
	if(itModel1->getCluster()==0 && itModel2->getCluster()==0){
	  itModel1->setCluster(numCluster);
	  itModel2->setCluster(numCluster);
	  list<GeneModel*> tmpList;
	  tmpList.push_back(&*itModel1);
	  tmpList.push_back(&*itModel2);
	  itMapClusterModel = mapClusterModel.find(numCluster);
	  if(itMapClusterModel == mapClusterModel.end())
	    mapClusterModel.insert(make_pair(numCluster,tmpList));
	  else{
	    (itMapClusterModel->second).push_back(&(*itModel1));
	    (itMapClusterModel->second).push_back(&(*itModel2));
	  }
	  ++numCluster;
	}
	else if(itModel1->getCluster()==itModel2->getCluster())
	  continue;
	else if(itModel1->getCluster() == 0 ){
	  itModel1->setCluster(itModel2->getCluster());
	  mapClusterModel[itModel2->getCluster()].push_back(&*itModel1);
				}
	else if(itModel2->getCluster() == 0){
	  itModel2->setCluster(itModel1->getCluster());
	  mapClusterModel[itModel1->getCluster()].push_back(&*itModel2);
	}
	else
	  this->fusionCluster(itModel1->getCluster(),itModel2->getCluster(),mapClusterModel);
      }
    }
  }
  for(map<s32,list<GeneModel*> >::iterator itMap = mapClusterModel.begin() ; itMap != mapClusterModel.end() ; ++itMap){
    itMap->second.sort(sortCdsB);
  }
  for(map<s32,list<GeneModel*> >::iterator itMap = mapClusterModel.begin() ; itMap != mapClusterModel.end() ; ++itMap){
    if(itMap->first == 0 )//at map[0] model with no clusters
      continue;
    itMap->second.sort(sortCdsB);
    list<GeneModel*> keepModel;
    for(list<GeneModel*>::iterator itList = itMap->second.begin() ; itList != itMap->second.end() ; ++itList){
      if(itList == itMap->second.begin()  ){
	(*itList)->setToDelete(false);
	keepModel.push_back(*itList);
      }
      else{
	for(list<GeneModel*>::iterator itKeepModel = keepModel.begin() ; itKeepModel != keepModel.end(); ++itKeepModel){
	  if((*itList)->overlapOrf((*itKeepModel)->getCDS())){
	    (*itList)->setToDelete(true);
	    break;
	  }
	}
	if((*itList)->getToDelete() == false)
	  keepModel.push_back(*itList);
      }
    }
  }
}

void GeneModelList::fusionCluster(s32 formerCluster, s32 newCluster, map<s32,list<GeneModel*> >& mapClusterModel){
	for(list<GeneModel>::iterator itModel1 = _models.begin(); itModel1 != _models.end() ; ++itModel1){
		if(itModel1->getCluster() == formerCluster){
			itModel1->setCluster(newCluster);
		}
	}
	list<GeneModel*> tmpList = mapClusterModel[formerCluster];
	mapClusterModel.erase(formerCluster);
	for(list<GeneModel*>::iterator itTmp = tmpList.begin() ; itTmp != tmpList.end(); ++itTmp){
		mapClusterModel[newCluster].push_back(*itTmp);
	}
}

void GeneModelList::deleteSmallCDS(s32 min_size_cds){
	for(GeneModelL::iterator itModels = _models.begin(); itModels != _models.end();++itModels){
		if(itModels->getToDelete())
			continue;

		if(itModels->getCdsSize() < min_size_cds)
				itModels->setToDelete(true);
	}
}


void GeneModelList::ratioCdsUtr(){
  for(GeneModelL::iterator itModels = _models.begin(); itModels != _models.end();++itModels){
    if(itModels->getToDelete())
      continue;
    f4 ratio =itModels->getCdsSize()*1.0 / itModels->getModelSize();
    if(itModels->getExons().size() == 1 )
      if(ratio < 0.8)
	itModels->setToDelete(true);
  }
}

s32 GeneModelList::printOut(ofstream& ofs_modelsFilename, bool formatDef){
  s32 modelsPrint = 0;
  for(GeneModelL::iterator itModels = _models.begin(); itModels != _models.end() ; ++itModels){
    if(!itModels->getToDelete()){
      ++modelsPrint;
      //	cout << "CDS in printOUT " << itModels->getCDS().first << " " << itModels->getCDS().second << endl;
      //	cout << *itModels << endl;
      itModels->printAnnot(ofs_modelsFilename, formatDef);
    }
  } 
  //	ofs_modelsFilename.close();
  return modelsPrint;
}

void GeneModelList::filter(bool ratio,s32 min_size_cds,char* longReadsFilename){

	    this->includedMono();
	    this->deleteSmallCDS(min_size_cds);

	    if(longReadsFilename == NULL)
	    	this->clusterLocation();
	    if(ratio)
	    	this->ratioCdsUtr();
}
bool sortCds(const GeneModel & model1, const GeneModel & model2) {
	if(model1.getCdsSize() > model2.getCdsSize() )
		return true;
	else if(model1.getCdsSize() == model2.getCdsSize() && model1.getModelSize() >= model2.getModelSize() )
		return true;
	else
		return false;
}
bool sortCdsB(const GeneModel* model1, const GeneModel* model2) {
	if(model1->getCdsSize() > model2->getCdsSize() )
			return true;
		else if(model1->getCdsSize() == model2->getCdsSize() && model1->getModelSize() >= model2->getModelSize() )
			return true;
		else
			return false;
}

