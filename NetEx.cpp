/*******************************************************************************
+
+  NetEx.cpp
+
+  
+  
+ 
 *******************************************************************************/

#include "NetEx.h"

// void NetEx::getStats (){
//   ofstream filestat;
//   filestat.open("models.stats");
//   cerr<<"sequence: "<<getFilename()<<endl;
  
//   Graph g_write=getGraph();
//   IndexMap index = get(vertex_index, g_write);
//   Components components=getComponents();
//   Components::value_type::iterator j,jend;
//   s32 cb=0,cp=0,splices=0;
//   map<s32,s32> taille,splice;
//   for (Components::size_type c = 0; c < components.size(); ++c) {//for each graph component (=gene)
    
//     for (j=components[c].begin(); j != components[c].end(); ++j){//for each component element (=exon)
//       //cout<<"noeud "<<*j<<endl;
//       s32 numout=out_degree(index[*j],g_write);
//       if (numout>1){
// 	//splice[numout]++;
// 	splices++;
//       }
//       cp++;
//     }
//     if (cp>1){taille[cp]++;cb++;}
//     cp=0;
//   }
//   cerr<<"gene-models count: "<<cb<<endl<<"splice events: "<<splices<<endl<<endl<<"models size (number of exons : value)"<<endl;
//   map<s32,s32>::iterator itm;
//   for(itm=taille.begin();itm!=taille.end();itm++){
//     cerr<<(*itm).first<<" : "<<(*itm).second<<endl;
//   }
//   cerr<<endl<<endl;
// }


//** graph outputting (dot format)  **//
void NetEx::graphOut (string filename){
  stringstream ss;
  ss << "Model.dot";
  filename += ss.str();
  filename = "dot_files/"+filename;
  ofstream ofs( filename.c_str() );
  write_graphviz(ofs, _graph,make_label_writer(get(&MyVertex::name, _graph) ) , make_label_writer(get(edge_weight, _graph)));
}

// to count the number of vertices (exons) and edges (introns) in the graph
void NetEx::countVerticesAndEgdes(list<s32>& comp, s32& nbVertices, s32& nbEdges) {
	nbVertices = comp.size();
	for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++)
		nbEdges += out_degree(*j,_graph);
}

list<s32> NetEx::sourcesNodes(list<s32>comp){
	 list<s32> sources;
	 for (list<s32>::iterator j = comp.begin(); j != comp.end() ; j++) {
		 if(!in_degree(*j, _graph))
			sources.push_back(*j);
	 }
	 return sources;
}
list<s32> NetEx::endNodes(list<s32> comp){
	list<s32> end;
	 for (list<s32>::iterator j = comp.begin(); j != comp.end() ; j++) {
	    if (!out_degree(*j, _graph)){
			end.push_back(*j);
	    }
	 }
	 return end;
}
// to count the number of paths in a connected component
s32 NetEx::count_allPaths(list<s32>& comp) {
  list<s32> sources;
  map<s32, s32> subs;
  
  for (list<s32>::iterator j = comp.begin(); j != comp.end() ; j++) {
    if (!out_degree(*j, _graph))
      subs[*j] = 1;

    if(!in_degree(*j, _graph))
    	sources.push_back(*j);
  }
  s32 total = 0;
  for (list<s32>::iterator itm = sources.begin() ; itm != sources.end(); itm++) { 
    s32 tmp = count_paths(*itm, subs);
    if(tmp!=-1) total += tmp;
    else return -1;
  }
  return total;
}

// to count the number of paths 
s32 NetEx::count_paths(s32 node, map<s32, s32>& subs) {
  if (subs.find(node) != subs.end())
    return subs[node];

  pair<out_edge_it, out_edge_it> p;
  s32 total = 0;
  for (p = out_edges(node, _graph); p.first != p.second; ++p.first) {
    s32 tmp = count_paths(target(*p.first, _graph), subs);
    total += tmp;
    if (tmp == -1 || tmp > 100000000 || total > 100000000) return -1;
  }

  subs[node] = total;
  return total;
}

//NetEx::
// all paths finder 
list<list<s32> > NetEx::allPathsFinder(list<s32>& comp) {
  list<list<s32> > selectedPaths;
  list<s32> sources;
  // AJOUT MAX
  map<s32, list<list<s32> > > solved;
  list<s32> tlist(1, 0);
  for (list<s32>::iterator j = comp.begin(); j != comp.end(); ++j) {
    if(!out_degree(*j, _graph)) {
		// AJOUT MAX
		// Ajout des noeuds terminaux comme solutions
		s32 idVertex = nodeToVertex(*j);
		tlist.front() = idVertex; //FIXME Why is it a list ?
		solved[idVertex].push_back(tlist);
    }
  }
  sources = sourcesNodes(comp);
  // Liste vide indiquant le chemin parcouru
  list<s32> start;
  for (list<s32>::iterator itm = sources.begin(); itm != sources.end(); ++itm) {
    start.clear();
    compute_paths(*itm, start, solved);
  }

  // AJOUT MAX :
  // Fusion des chemins pour chaque source => résultat
  for (list<s32>::iterator it=sources.begin(); it != sources.end(); ++it) {
	// On utilise splice car on a plus besoin de solved => on ne créé donc pas de copies
	s32 idNode = nodeToVertex(*it);
	selectedPaths.splice(selectedPaths.end(), solved[idNode]);
  }
  return selectedPaths;
}

void NetEx::printMap (map<s32,list < list<s32> > > mapP){
	cout <<"print map map.size() " << mapP.size()<< " " ;
	for(map<s32,list < list<s32> > >::iterator it =  mapP.begin(); it != mapP.end(); ++it){
		cout << " ["<<it->first << "] " ;
		for(list< list <s32> >:: iterator itT = it->second.begin(); itT != it->second.end(); ++itT){
			for(list <s32 >::iterator itTT = itT->begin() ; itTT != itT->end(); ++itTT){
				cout << *itTT <<" ";
			}
			cout <<  "\t\\\t";
		}
		cout << endl;
	}
	cout << endl;
}

s32 NetEx::nodeToVertex(s32 id){
	TSSRList::iterator it1 = _vertices->begin();
	advance(it1,id);
//	cout <<"print in nodeToVertex " << *(*it1) << endl;
	s32 childId = (*it1)->getID();
	return childId;
}
void NetEx::bfs(s32 startId,map<s32,list < list<s32> > >& paths , map<s32,s32> vertexToBGL,s32 currentColor,list<s32>& colorNotAllowed,std::queue<s32> & Qbfs  ){
	cout << "enter in BFS" << endl;
	out_edge_it ei, ei_end;
	bool boolRecurciveCall = false;
	time_t before,after,beforeFor;
	while (! Qbfs.empty()) {
		cout <<"enter in while" <<endl;
		cout << "QBFS " << Qbfs.size() << endl;
		s32 parentId = Qbfs.front(); Qbfs.pop();

		map<s32,s32>::iterator itVB = vertexToBGL.find(parentId);
		if(itVB == vertexToBGL.end()){
			cerr << "error do not find id in vertex "<< parentId<<endl;
			exit(1);
		}
		s32 idNodeParent = itVB->second;
			map<s32,list<list <s32> > >::iterator itPathParent = paths.find(parentId);
			/*			cout <<"number of out edges " << out_degree(idNodeParent,_graph)<< " itPathParent " << parentId << " path end "<< itPathParent->first<< " paths " << paths.size() << endl;

			if(out_degree(idNodeParent,_graph) == 0 && itPathParent == paths.end()){//for mono exoniques
			cout <<" mono " << endl;
			list<s32> newPath;
			list<list<s32> > newList;
			newPath.push_back(parentId);
			newList.push_back(newPath);
			paths.insert(make_pair(parentId,newList));
		}
*/		before = time(NULL);

		for (tie(ei, ei_end) = out_edges(idNodeParent, _graph); ei != ei_end; ++ei) {
			cout <<"loop on out edges "<<endl;
			s32 childId = nodeToVertex(target(*ei, _graph));
			//check the color before doing anything
			//Which color is my child node ?
			TSSRList::iterator it1 = _vertices->begin();
			advance(it1,target(*ei, _graph));
			map<s32,bool> mapIdTranscrit = (*it1)->getIdTranscrit();
			map<s32,bool>::iterator itColorChild = mapIdTranscrit.find(currentColor);

			if(itColorChild == mapIdTranscrit.end()){
				bool boolColorNotAllowed = false;
				beforeFor = time(NULL);
				for(list<s32>::iterator itList = colorNotAllowed.begin(); itList != colorNotAllowed.end();++itList){
					map<s32,bool>::iterator itNotAllowed = mapIdTranscrit.find(*itList);
					if(itNotAllowed != mapIdTranscrit.end()){
						boolColorNotAllowed = true;
						break;
					}
				}
				after = time(NULL);
			if(PRINTTIME)	cout << " time allPathFinder for colorNotAllowed " << difftime(after,beforeFor)<< endl;
				if(boolColorNotAllowed)
					continue;
				else{
					colorNotAllowed.push_back(currentColor);
					boolRecurciveCall = true;
				}

			}
			map<s32,list<list <s32> > >::iterator itPathChild = paths.find(childId);
			if(itPathChild != paths.end()){ // update Qbfs
				beforeFor = time(NULL);
				for(list<list <s32> >::iterator itListPathChild = itPathChild->second.begin(); itListPathChild != itPathChild->second.end();++itListPathChild){
					list<s32>::iterator itLL;
					itLL = std::find(itListPathChild->begin(), itListPathChild->end(), parentId);
					if(itLL != itListPathChild->end())
						break;
					else
						Qbfs.push(childId);
				}
				after = time(NULL);
			if(PRINTTIME)	cout << " time allPathFinder for itPathChild" << difftime(after,beforeFor)<< endl;
			}
			else
				Qbfs.push(childId);

			//Start BFS
			if(itPathChild == paths.end()){//if PathChild doesnt exist
				cout << "no path child " << endl;
				list<list<s32> > newList;
				if(itPathParent != paths.end()){//if Parent contains a path
					beforeFor = time(NULL);
					for(list<list <s32> >::iterator itListPathParent = itPathParent->second.begin(); itListPathParent != itPathParent->second.end();++itListPathParent){
						list<s32> newPath;
						for(list<s32>::iterator itL = itListPathParent->begin() ; itL != itListPathParent->end() ; ++itL){
							newPath.push_back(*itL);
						}
						newPath.push_back(childId);
						newList.push_back(newPath);
					}
					after = time(NULL);
				if(PRINTTIME)	cout << " time allPathFinder for itPathParent" << difftime(after,beforeFor)<< endl;
					paths.insert(make_pair(childId,newList));
				}
				else{
					list<s32> newPath;
					newPath.push_back(parentId);
					newPath.push_back(childId);
					newList.push_back(newPath);
					if (paths.find(childId) != paths.end()){
						cout << "Error child id already in predListList " <<endl;
						exit(1);
					}
					paths.insert(make_pair(childId,newList));
				}
			}
			else{
				cout <<"path child " << endl;
				beforeFor = time(NULL);
				list<list <s32> > copyList(paths[parentId]);
				for(list<list <s32> >::iterator itCopyList = copyList.begin(); itCopyList != copyList.end();++itCopyList){
					itCopyList->push_back(childId);
				}
				itPathChild->second.merge(copyList);
				itPathChild->second.unique();
				after = time(NULL);
				if(PRINTTIME) cout << " time allPathFinder for itCopyList" << difftime(after,beforeFor)<< endl;
			}
			if(boolRecurciveCall){
				cout <<"recurcive call of searchColor " << endl;
				searchAllColor(target(*ei, _graph), paths,colorNotAllowed);
			}

		}
		after = time(NULL);
		if(PRINTTIME)	cout << " time allPathFinder BFS " << difftime(after,before)<< endl;
	}//While

}



void NetEx::dfs(s32 startId,map<s32,list < list<s32> > >& paths , map<s32,s32> vertexToBGL,s32 currentColor,list<s32>& colorNotAllowed,std::stack<s32>& Sdfs ){
	out_edge_it ei, ei_end;
//	std::stack<s32> Sdfs ;
//	Sdfs.push(startId);

	map<s32,s32> mapCurrentColor;
	mapCurrentColor.insert(make_pair(startId,currentColor));
	while (! Sdfs.empty()) {
		s32 parentId = Sdfs.top(); Sdfs.pop();
		map<s32,s32>::iterator itVB = vertexToBGL.find(parentId);
		if(itVB == vertexToBGL.end()){
			cerr << "error do not find id in vertex "<< parentId<<endl;
			exit(1);
		}
		s32 idNodeParent = itVB->second;
		map<s32,list<list <s32> > >::iterator itPathParent = paths.find(parentId);


		if(out_degree(idNodeParent,_graph) == 0 && itPathParent == paths.end()){//for mono exoniques
			list<s32> newPath;
			list<list<s32> > newList;
			newPath.push_back(parentId);
			newList.push_back(newPath);
			paths.insert(make_pair(parentId,newList));
		}
	//	currentColor = mapCurrentColor[parentId];
		for (tie(ei, ei_end) = out_edges(idNodeParent, _graph); ei != ei_end; ++ei) {
			s32 childId = nodeToVertex(target(*ei, _graph));
			//check the color before doing enything
			//Which color is my child node ?
			TSSRList::iterator it1 = _vertices->begin();
			advance(it1,target(*ei, _graph));
			map<s32,bool> mapIdTranscrit = (*it1)->getIdTranscrit();
			cout << "\n\n\t\t * parent " << parentId << " child " << childId <<" current color " << currentColor << endl;
			cout << " sources map id transcrit " <<childId << " : map size " << mapIdTranscrit.size() << endl;
			for(map<s32,bool>::iterator itMap = mapIdTranscrit.begin() ; itMap != mapIdTranscrit.end(); ++itMap){
					cout <<"color child " <<  itMap->first << endl;
			}
			map<s32,bool>::iterator itColorChild = mapIdTranscrit.find(currentColor);

			if(itColorChild == mapIdTranscrit.end()){
				bool boolColorNotAllowed = false;
				for(list<s32>::iterator itList = colorNotAllowed.begin(); itList != colorNotAllowed.end();++itList){
					map<s32,bool>::iterator itNotAllowed = mapIdTranscrit.find(*itList);
					if(itNotAllowed != mapIdTranscrit.end()){
						cout << " no correspondance between colorChild and currentColor, stop "<<endl;
						cout <<itNotAllowed->first << " " << *itList << " break " << endl;
						boolColorNotAllowed = true;
						break;
					}
				}
				if(boolColorNotAllowed)
					continue;
				else{
					cout << "add color to notAllow " << currentColor << " new currentColor is "<< mapIdTranscrit.begin()->first<<endl;
					colorNotAllowed.push_back(currentColor);
					//call recurcive
					cout << "call recurcive searchColor childId " << childId << endl;
					searchAllColor(target(*ei, _graph), paths,colorNotAllowed);
				//	currentColor = mapIdTranscrit.begin()->first; //FIXME how to assign currentColor if multiple color in map ?
				//	cout << "new current color " << currentColor << endl;
				}

			}
			else{
				cout << " child Color correspond to currentColor, continue script " << endl;
			}
			map<s32,list<list <s32> > >::iterator itPathChild = paths.find(childId);
			if(itPathChild != paths.end()){ // update Qbfs
				for(list<list <s32> >::iterator itListPathChild = itPathChild->second.begin(); itListPathChild != itPathChild->second.end();++itListPathChild){
					list<s32>::iterator itLL;
					itLL = std::find(itListPathChild->begin(), itListPathChild->end(), parentId);
					if(itLL != itListPathChild->end())
						break;
					else{
						Sdfs.push(childId);
						cout << "add to stack " << childId << endl;
					}

				}
			}
			else{
				Sdfs.push(childId);
				cout << "add to stack "<< childId << endl;
			}


			//Start BFS
			if(itPathChild == paths.end()){//if PathChild doesnt exist
				list<list<s32> > newList;
				if(itPathParent != paths.end()){//if Parent contains a path
					for(list<list <s32> >::iterator itListPathParent = itPathParent->second.begin(); itListPathParent != itPathParent->second.end();++itListPathParent){
						list<s32> newPath;
						for(list<s32>::iterator itL = itListPathParent->begin() ; itL != itListPathParent->end() ; ++itL){
							newPath.push_back(*itL);
						}
						newPath.push_back(childId);
						newList.push_back(newPath);
					}
					paths.insert(make_pair(childId,newList));
				}
				else{
					list<s32> newPath;
					newPath.push_back(parentId);
					newPath.push_back(childId);
					newList.push_back(newPath);
					if (paths.find(childId) != paths.end()){
						cout << "Error child id already in predListList " <<endl;
						exit(1);
					}
					paths.insert(make_pair(childId,newList));
				}
			}
			else{
				list<list <s32> > copyList(paths[parentId]);
				for(list<list <s32> >::iterator itCopyList = copyList.begin(); itCopyList != copyList.end();++itCopyList){
					itCopyList->push_back(childId);
				}
				itPathChild->second.merge(copyList);
				itPathChild->second.unique();
			}
			mapCurrentColor.insert(make_pair(parentId,currentColor));
		}
	}//While
}
map<s32,s32> NetEx::mapIdsVertextoBGL(){
	map<s32,s32> vertexToBGL;
	vertex_it currentNode,endNode;
	for (tie(currentNode,endNode) = boost::vertices(_graph); currentNode!=endNode; ++currentNode){
		s32 idNode = nodeToVertex(*currentNode);
		if(vertexToBGL.find(idNode) != vertexToBGL.end()){
			cerr << "errors ids already known " << endl;
			exit(1);
		}
		vertexToBGL.insert(make_pair(idNode,*currentNode));
	}
return vertexToBGL;
}


list<list<s32> > NetEx::PathsFinderWithCondition(list<s32> comp){
	list<s32> sources = sourcesNodes(comp);
	map<s32,list < list<s32> > >predM;
	list<s32> endList = endNodes(comp);
	list<list<s32> > paths;
	time_t before, after;
	for(list<s32>::iterator itSources = sources.begin() ; itSources != sources.end() ; ++itSources){
		cout <<"PathsFinderWithCondition Sources " << *itSources  << endl;
		predM.clear();
		list<s32> colorNotAllowed;
		before = time(NULL);
		searchAllColor(*itSources, predM,colorNotAllowed);
		after = time(NULL);
	if(PRINTTIME)	cout << " time allPathFinder searchAllColor " << difftime(after,before)<< endl;
//		printMap(predM);
		before = time(NULL);
		for(list<s32>::iterator itEnd = endList.begin() ; itEnd != endList.end() ; ++itEnd){
			s32 childId = nodeToVertex(*itEnd);
			map<s32,list < list<s32> > >::iterator itPred = predM.find(childId);
			paths.splice(paths.end(), predM[childId]);
		}
		after = time(NULL);
	if(PRINTTIME)	cout << " time allPathFinder for loop " << difftime(after,before)<< endl;
	}
	before = time(NULL);
	paths.sort();
	after = time(NULL);
	if(PRINTTIME) cout << " time allPathFinder sort " << difftime(after,before)<< endl;
	before = time(NULL);
	paths.unique();
	after = time(NULL);
	if(PRINTTIME) cout << " time allPathFinder unique " << difftime(after,before)<< endl;
	 return paths;
}

void NetEx::searchAllColor(s32 idSource,map<s32,list < list<s32> > >& predM,list<s32>& colorNotAllowed){
	s32 idVertex = nodeToVertex(idSource);
	TSSRList::iterator it1 = _vertices->begin();
	advance(it1,idSource);
	SSRContig* ctg = *it1;
	map<s32,s32> vertexToBGL = mapIdsVertextoBGL();
	map<s32,bool> mapIdTranscrit = ctg->getIdTranscrit();
//	std::stack<s32> Sdfs ;
	std::queue<s32> Qbfs ;

//	Sdfs.push(idVertex);
	time_t before = time(NULL);


			map<s32,list<list <s32> > >::iterator itPathParent = predM.find(idVertex);
			cout <<"number of out edges " << out_degree(idSource,_graph)<< " itPathParent " << idVertex << " path end "<< itPathParent->first<< " paths " << predM.size() << endl;
			if(out_degree(idSource,_graph) == 0 && itPathParent == predM.end()){//for mono exoniques
				cout <<" mono " << endl;
				list<s32> newPath;
				list<list<s32> > newList;
				newPath.push_back(idVertex);
				newList.push_back(newPath);
				predM.insert(make_pair(idVertex,newList));
			}
			else if(out_degree(idSource,_graph) != 0 ){
				for(map<s32,bool>::iterator itMap = mapIdTranscrit.begin() ; itMap != mapIdTranscrit.end(); ++itMap){
				cout << "current color "<< itMap->first<<" idVertex "<< idVertex<< " predM " << predM.end()->first << endl;
					s32 currentColor = itMap->first;
						Qbfs.push(idVertex);
						bfs(idVertex,predM , vertexToBGL, currentColor,colorNotAllowed,Qbfs);
				}

	//	dfs( idVertex,predM, vertexToBGL,currentColor,colorNotAllowed,Sdfs);
	}
	time_t after = time(NULL);
if(PRINTTIME)	cout << " time allPathFinder loop on bfs " << difftime(after,before)<< endl;
}


// NOUVEAU allpath - to compute paths a lot faster than before
void NetEx::compute_paths(s32 node, list<s32>& path, map<s32, list<list<s32> > >& solved) {
  list<s32> tmp_path = path; 
  // On est arrivé sur un noeud déjà connu ?
  s32 idVertex = nodeToVertex(node);
  if (solved.find(idVertex) != solved.end()) {
    list<s32> final_path;
    list<s32> *path_ref;
    s32 start_node;
    list<list<s32> >::iterator solved_path;
    // On ajoute tous les sous chemins à la liste résolue
    while(tmp_path.begin() != tmp_path.end()) {
     start_node = tmp_path.front();
      // On concatène le chemin courant avec la solution trouvée
      for (solved_path=solved[idVertex].begin(); solved_path != solved[idVertex].end(); ++solved_path) {
		solved[start_node].push_back(tmp_path);
		path_ref = &solved[start_node].back();
		path_ref->insert(path_ref->end(), solved_path->begin(), solved_path->end());
      }
      tmp_path.pop_front();
    }
  }
  else {
    // Ce qui se passe quand on est pas sur un noeud déjà connu
    tmp_path.push_back(idVertex);
    for ( pair<out_edge_it, out_edge_it> p = out_edges(node, _graph) ; p.first != p.second ; ++p.first ) {
      compute_paths(target(*p.first, _graph), tmp_path, solved);
    }
  }
}

// test allpath
// to retrieve all connected components of the graph
vector<list<s32> > NetEx::getComponents() {
  GraphU gU(_vertices->size());
  pair<edge_it, edge_it> p = edges(_graph);
  for(edge_it it = p.first ; it != p.second ; it++) {
	 add_edge(source(*it, _graph), target(*it, _graph), gU);
  }
  std::vector<int> component(num_vertices(gU));
  _nb_connected_components = connected_components(gU, &component[0]);
  vector<list<s32> > cc(_nb_connected_components);
  for (s32 i=0 ; i < (s32)component.size() ; i++){
	  cc[component[i]].push_back(i);
  }
  return cc;
}


//to find out whether the graph has cycles..
s32 NetEx::nbCycle() {
  bool has_cycle = false;
  int count_cycle=0;
  cycle_detector vis(has_cycle, count_cycle);
  depth_first_search(_graph, visitor(vis));
  return count_cycle;
}

void NetEx::tagExons(map<string,s32>& inExon, map<s32,TSSRList>& mapStartExon,map<s32,TSSRList>&mapEndExon,TSSRList& monoExon){
	for(TSSRList::iterator itComponent = _vertices->begin(); itComponent != _vertices->end(); ++itComponent){
		SSRContig* ctg = *itComponent;
		if(!in_degree(ctg->getID(), _graph) && !out_degree(ctg->getID(), _graph))
			monoExon.push_back(ctg);
		else if(!in_degree(ctg->getID(), _graph))
			mapStartExon[(*itComponent)->end()].push_back(*itComponent);

		else if(!out_degree(ctg->getID(), _graph))
			mapEndExon[(*itComponent)->start()].push_back(*itComponent);
		else{
			 SSRContig* exon = *itComponent;
			 ostringstream oss;
			 oss << exon->seqName() << "@" << exon->start() << "@" << exon->end() << "@" << exon->strand();
			 string key = oss.str();
			 inExon.insert(make_pair(key,exon->getID()));
		}
	}
}


void NetEx::cleanGraph(){
	 cerr << "\tBefore cleaning, number vertices " << num_vertices(_graph) << " number edges "<< num_edges(_graph)<<endl;
	//it will clean the component regarding the start and end exons.
	s32 left,right,left2,right2;
	//Tag exons
	int deleteMono = 0;
	TSSRList listMono; // OutExon = exons at begin or end of the graph, inExon = exons intern of the graph
	map<s32,TSSRList> mapStartExon, mapEndExon;
	map<string,s32> inExon;
	tagExons(inExon,mapStartExon,mapEndExon,listMono); //TODO change tag Exon
	//TODO sort exon by pos and length !

	//fusionne les exons start
	for(map<s32,TSSRList>::iterator itMap = mapStartExon.begin(); itMap != mapStartExon.end(); ++itMap){
		for(TSSRList::iterator itCurrent = itMap->second.begin() ; itCurrent != itMap->second.end();++itCurrent){
			SSRContig* ctgStart = *itCurrent;
			if ( _graph[ctgStart->getID()].guid.deleteVertex == true) continue;
			for(TSSRList::iterator itNext = ++itCurrent ; itNext  != itMap->second.end();++itNext){
				SSRContig* ctgNext = *itNext;
				if(ctgStart->strand() != ctgNext->strand()) continue;
				if(ctgStart->end() == ctgNext->end()){
					if(ctgStart->size() > ctgNext->size() ){ // we keep current and delete next
						for (pair<out_edge_it, out_edge_it> pEdges = out_edges(ctgNext->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
							s32 outEdge = target(*pEdges.first, _graph);
							edge_t e; bool b;
							tie(e,b) = add_edge(ctgStart->getID(),outEdge, _graph); //For the moment id in vertices and graph are equal
							s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
							formerWeight += get(edge_weight_t(), _graph,e);
							put(edge_weight_t(), _graph, e, formerWeight);
						}
						clear_vertex(ctgNext->getID(),_graph); //clear all edges before removing the vertex
						_graph[ctgNext->getID()].guid.deleteVertex = true;
					}
					else {//we delete the current
						for (pair<out_edge_it, out_edge_it> pEdges = out_edges(ctgStart->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
							s32 outEdge = target(*pEdges.first, _graph);
							edge_t e; bool b;
							tie(e,b) = add_edge(ctgNext->getID(),outEdge, _graph);
							s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
							formerWeight += get(edge_weight_t(), _graph,e);
							put(edge_weight_t(), _graph, e, formerWeight);

						}
						clear_vertex(ctgStart->getID(),_graph); //clear all edges before removing the vertex
						_graph[ctgStart->getID()].guid.deleteVertex = true;
					}
				}
			}
			--itCurrent;
			if(_graph[ctgStart->getID()].guid.deleteVertex == true) continue;
			//Fusionne mono FIXME Good idea to do it here ?
			for(TSSRList::iterator itMono = listMono.begin(); itMono != listMono.end(); ++itMono){
				if((*itMono)->strand() != ctgStart->strand()) continue;
				if((*itMono)->start() <= ctgStart->start() && (*itMono)->end()<= ctgStart->end() && overlap(ctgStart,*itMono) ){
					ctgStart->setStart((*itMono)->start());
					 _graph[(*itMono)->getID()].guid.deleteVertex = true;
					 ++deleteMono;
				}
			}
		}
	}
//	ostringstream convert;
//					convert << itMap->first;
//					string tmp = convert.str();
//				this->graphOut("tmp");

	// simplification pour endExon
	for(map<s32,TSSRList>::iterator itMap = mapEndExon.begin(); itMap != mapEndExon.end(); ++itMap){
		for(TSSRList::iterator itCurrent = itMap->second.begin() ; itCurrent != itMap->second.end();++itCurrent){
			SSRContig* ctgStart = *itCurrent;
			if ( _graph[ctgStart->getID()].guid.deleteVertex == true) continue;
			for(TSSRList::iterator itNext = ++itCurrent ; itNext  != itMap->second.end();++itNext){
				SSRContig* ctgNext = *itNext;
			//	cout << "Change Next " << endl;
				if(ctgStart->strand() != ctgNext->strand()) continue;
				leftRight(ctgStart,left,right);
				leftRight(ctgNext,left2,right2);
				if(overlap(ctgStart, ctgNext) && ctgStart->start() == ctgNext->start()){
					if(ctgStart->size() > ctgNext->size()){ // we keep current and delete next
				//		cout << "endExon start "<< ctgStart->start() << " " << ctgStart->end() << "\t";
				//		cout << " end " << ctgNext->start() << " " << ctgNext->end() << endl;
						for (pair<in_edge_it, in_edge_it> pEdges = in_edges(ctgNext->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
							s32 inEdge = source(*pEdges.first, _graph);
							edge_t e; bool b;
				//			cout << "inEdge " << inEdge << " ctgStart->getID() " << ctgStart->getID() << " start "<< ctgStart->start() << " end " << ctgStart->end() << endl;
							tie(e,b) = add_edge(inEdge,ctgStart->getID(), _graph);
							s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
				//			cout << "formerWeight "<< formerWeight << "\t";
							formerWeight += get(edge_weight_t(), _graph,e);
				//			cout << " newWeight " << formerWeight << endl;
							put(edge_weight_t(), _graph, e, formerWeight);
				//			ostringstream convert;
				//			convert << inEdge;
				//			string tmp = convert.str();
				//		this->graphOut(tmp);
						}
					 clear_vertex(ctgNext->getID(),_graph); //clear all edges before removing the vertex
					 _graph[ctgNext->getID()].guid.deleteVertex = true;

					}
					else{ //if (deleteList.find(itCurrent) == deleteList.end()){//we delete the current
				//		cout << "endExon end "<< ctgNext->start() << " " << ctgNext->end() << endl;
						for (pair<in_edge_it, in_edge_it> pEdges = in_edges(ctgStart->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
							s32 inEdge = source(*pEdges.first, _graph);
							edge_t e; bool b;
							tie(e,b) = add_edge(inEdge,ctgNext->getID(), _graph);
							s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
							formerWeight += get(edge_weight_t(), _graph,e);
							put(edge_weight_t(), _graph, e, formerWeight);
						}
						clear_vertex(ctgStart->getID(),_graph); //clear all edges before removing the vertex
						_graph[ctgStart->getID()].guid.deleteVertex = true;
					}
			 }
		 }
			--itCurrent;
			if(_graph[ctgStart->getID()].guid.deleteVertex == true)
				continue;
			//Fusionne mono FIXME Good idea to do it here ?
			for(TSSRList::iterator itMono = listMono.begin(); itMono != listMono.end(); ++itMono){
				if((*itMono)->strand() != ctgStart->strand())
					continue;
				if((*itMono)->start() >= ctgStart->start() && (*itMono)->end()>= ctgStart->end() && overlap(ctgStart,*itMono)){
					ctgStart->setEnd((*itMono)->end());
					_graph[(*itMono)->getID()].guid.deleteVertex = true;
					++deleteMono;
					}
			}
		 if(_graph[ctgStart->getID()].guid.deleteVertex == true)
			 continue;
		 //fusionne stat/end exons
		 for(map<s32,TSSRList>::iterator itMapStart = mapStartExon.begin(); itMapStart != mapStartExon.end(); ++itMapStart){
			 for(TSSRList::iterator itCurrentStart = itMapStart->second.begin() ; itCurrentStart != itMapStart->second.end();++itCurrentStart){
				 if(_graph[ (*itCurrentStart)->getID()].guid.deleteVertex == true)
					 continue;
		 		 if(!overlap(*itCurrent,*itCurrentStart) || (*itCurrent)->strand() != (*itCurrentStart)->strand())
		 			 continue;
	 			 if(ctgStart->start() <= (*itCurrentStart)->start() && ctgStart->end() <= (*itCurrentStart)->end() &&ctgStart->end() !=  (*itCurrentStart)->start() ){
		 			 SSRContig* newCtg = new SSRContig(*ctgStart);//FIXME memory leak
		 			 newCtg->setEnd((*itCurrentStart)->end());
		 			 newCtg->setID(_vertices->size());
		 			 ostringstream oss;
		 			 oss << newCtg->seqName() << "@" << newCtg->start() << "@" << newCtg->end() << "@" << newCtg->strand();
		 			 string key = oss.str();
		 			 map<string,s32>::iterator itKey = inExon.find(key);
		 			 if( itKey == inExon.end()){
		 				 inExon.insert(make_pair(key,true));
		 				 _vertices->push_back(newCtg);
		 				 Vertex u = boost::add_vertex(_graph);
		 				 _graph[u].vertex = newCtg;
		 				 for (pair<out_edge_it, out_edge_it> pEdges = out_edges((*itCurrentStart)->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
		 					 s32 outEdge = target(*pEdges.first, _graph);
		 					 edge_t e; bool b;
		 					 tie(e,b) = add_edge(newCtg->getID(),outEdge, _graph); //For the moment id in vertices and graph are equal
		 					s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
		 					formerWeight += get(edge_weight_t(), _graph,e);
		 					put(edge_weight_t(), _graph, e, formerWeight);
		 				 }
		 				 for (pair<in_edge_it, in_edge_it> pEdges = in_edges(ctgStart->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
		 					 s32 inEdge = source(*pEdges.first, _graph);
		 					 edge_t e; bool b;
		 					 tie(e,b) = add_edge(inEdge,newCtg->getID(), _graph);
		 					s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
		 					formerWeight += get(edge_weight_t(), _graph,e);
		 					put(edge_weight_t(), _graph, e, formerWeight);
		 				}
		 			}
		 			else{
		 				for (pair<out_edge_it, out_edge_it> pEdges = out_edges((*itCurrentStart)->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
		 					s32 outEdge = target(*pEdges.first, _graph);
		 					edge_t e; bool b;
		 					tie(e,b) = add_edge(itKey->second,outEdge, _graph); //For the moment id in vertices and graph are equal
		 					s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
		 					formerWeight += get(edge_weight_t(), _graph,e);
		 					put(edge_weight_t(), _graph, e, formerWeight);
		 				}
		 				for (pair<in_edge_it, in_edge_it> pEdges = in_edges(ctgStart->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
		 					s32 inEdge = source(*pEdges.first, _graph);
		 					edge_t e; bool b;
		 					tie(e,b) = add_edge(inEdge,itKey->second, _graph);
		 					s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
		 					formerWeight += get(edge_weight_t(), _graph,e);
		 					put(edge_weight_t(), _graph, e, formerWeight);
		 				}
		 			}
		 		}
			 }
		 }
		}

	}
	 deleteNode();
	 cerr << "\tAfter cleaning, number vertices " << num_vertices(_graph) << " number edges "<< num_edges(_graph)<<endl;
}

void NetEx::deleteNode(){
	for(u32 vi = 0 ; vi < num_vertices(_graph); ){
		if(_graph[vi].guid.deleteVertex == true){
			TSSRList::iterator it1;
			it1 = _vertices->begin();
			int nb =vi;
			advance(it1,nb);
			_vertices->erase (it1);
			clear_vertex(vi,_graph);
			remove_vertex(vi, _graph);
		}
		else
			++vi;
	}

//	vertex_it b,e;
 //   			    					    	for (boost::tie(b,e) = boost::vertices(_graph); b!=e; ++b){
  //  			    					    		TSSRList::iterator it1 = _vertices->begin();
   // 			    					    				advance(it1,*b);
   // 			    					    		std::cout << "The vertex ID is: " <<  (*it1)->getID() << " id in graph " << *b<<endl;     			    					    	}


    			    					//    	exit(1);
}



void NetEx::simplifyBigGraph(list<s32>comp,s32 threshold){
	string seqname;
	map<Edge,s32> weight;
	map<int,int> degreeNode;
	s32 nbEdges=0, nbVertices=0;
	nbVertices = comp.size();
	s32 getWeight;
	for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++){
		int degreeBefore = degree(*j,_graph);
		degreeNode.insert(make_pair(*j,degreeBefore));
	}

		for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++){
			nbEdges += out_degree(*j,_graph);
			for (pair<out_edge_it, out_edge_it> pEdges = out_edges(*j, _graph); pEdges.first != pEdges.second; ++pEdges.first){
				getWeight =  get(edge_weight_t(), _graph,*pEdges.first);
				if(getWeight < threshold )
					remove_edge(*pEdges.first, _graph);
			}
		}
		for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++){
	//		TSSRList::iterator it1;
	//		it1 = _vertices->begin();
	//		int nb =*j;
	//		advance(it1,nb);
			int newDegree =  degree(*j,_graph);
			int degreeBefore = degreeNode.find(*j)->second;

			if(degreeBefore !=0 && newDegree ==0)
				_graph[*j].guid.deleteVertex = true;
		}
		deleteNode();


}
