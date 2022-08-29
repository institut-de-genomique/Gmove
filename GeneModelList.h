/*
 * GeneModelList.h
 *
 *  Created on: 28 oct. 2016
 *      Author: mdubarry
 */

#ifndef GENEMODELLIST_H_
#define GENEMODELLIST_H_

#include "GeneModel.h"

#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <list>

using namespace std;

typedef list<GeneModel> GeneModelL;
typedef map<string, GeneModel> MapModel;

struct get_second : public std::unary_function<MapModel::value_type, string>
{
  GeneModel operator()(const MapModel::value_type& value) const
  {
    return value.second;
  }
};

bool sortCds(const GeneModel &model1, const GeneModel& model2) ;
bool sortCdsB(const GeneModel* model1, const GeneModel* model2) ;


class GeneModelList{
 protected:
  GeneModelL _models;
  
 public:
  GeneModelList(){}
  GeneModelList(map<string, GeneModel> mapGeneModel){
    transform(mapGeneModel.begin(), mapGeneModel.end(), back_inserter(_models),  get_second() );//XXX This is slow !!
    _models.sort(sortCds);
  }
  ~GeneModelList(){}
  
  /* Methods */
  void deleteIncludedModel();
  void longestCDS();
  void bestScore();
  void catchModelOnUTR(GeneModel * longestModel);
  void insertModels(GeneModelList );
  void includedMono(); // a mono CDs overlap a pluri and the mono is < 300
  void clusterLocation();//Cluster location sur les CDS !!
  void fusionCluster(s32 formerCluster, s32 newCluster, map<s32,list<GeneModel*> >& mapClusterModel);
  void deleteSmallCDS(s32 min_size_cds);
  void ratioCdsUtr();
  s32 printOut(ofstream& ,bool);
  void filter(bool,s32,char* longReadsFilename);
  
  
  /* Accessors */
  s32 getSize()const {return _models.size();}
  
  GeneModelL getModels() const {return _models;}
  
  
};


#endif /* GENEMODELLIST_H_ */
