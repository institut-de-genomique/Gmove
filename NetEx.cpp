/*******************************************************************************
+
+  NetEx.cpp
+
+  
+  
+ 
 *******************************************************************************/

#include "NetEx.h"
//bool VERBOSE = 0;


//** graph outputting (dot format)  **//
void NetEx::graphOut (string directory, string filename){
  if (directory == ""){
    directory = "out";
  }
  filename = directory +"/"+ filename +"Model.dot";
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
    if(!in_degree(*j, _graph)){
      sources.push_back(*j);
    }
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
      //	s32 idVertex = nodeToVertex(*j);
      //	tlist.front() = idVertex; //FIXME Why is it a list ?
      //	solved[idVertex].push_back(tlist);
      tlist.front() = *j;
      solved[*j].push_back(tlist);
    }
  }
  sources = sourcesNodes(comp);
  // Liste vide indiquant le chemin parcouru
  list<s32> start;
  for (list<s32>::iterator itm = sources.begin(); itm != sources.end(); ++itm) {
    start.clear();
    compute_paths(*itm, start, solved);
    
    //time_t  t = time(0);
    //cerr << ctime(&t)<<endl;
  }
  
  // AJOUT MAX :
  // Fusion des chemins pour chaque source => résultat
  for (list<s32>::iterator it=sources.begin(); it != sources.end(); ++it) {
    // On utilise splice car on a plus besoin de solved => on ne créé donc pas de copies
    //s32 idNode = nodeToVertex(*it);
    //selectedPaths.splice(selectedPaths.end(), solved[idNode]);
    selectedPaths.splice(selectedPaths.end(), solved[*it]);
  }
  return selectedPaths;
}

list<list<s32> > NetEx::pathWithHigherWeight(list<s32>& comp){
  list<list<s32> > allPath;
  list <s32> sources;
  sources = sourcesNodes(comp);
  for (list<s32>::iterator startNode = sources.begin(); startNode != sources.end(); ++startNode) { //Find Path for each start Node
    map<s32,s32> previousNode;
    map<s32,s32> distance;
    initDijkstra(comp,*startNode,distance,previousNode);
    list<s32> Q = comp;
    
    
    while (!Q.empty()){
      s32 distS = Q.front();
      
      Q.pop_front();
      for ( pair<out_edge_it, out_edge_it> p = out_edges(distS, _graph) ; p.first != p.second ; ++p.first ) {
	updateDist(distS,target(*p.first, _graph),distance,previousNode);
      }
    }
    
    list<s32> stop = endNodes(comp);
    list<s32> path;
    for(list<s32>::iterator itStop = stop.begin(); itStop != stop.end();++itStop){
      path = shortestPath(*itStop,*startNode,previousNode);
      allPath.push_back(path);
    }
  }
  return allPath;
}

list<s32> NetEx::shortestPath(s32 stop,s32 start,map<s32,s32> previousNode){
  s32 currentNode = stop;
  list<s32> path;
  while (currentNode != start){
    if (currentNode == -1){
      path.clear();
      return path;
    }
    path.push_back(currentNode);
    currentNode = previousNode[currentNode];
  }
  path.push_back(start);
  return path;
}

void NetEx::initDijkstra(list<s32>comp,s32 startNode,map<s32,s32>& distance, map<s32,s32>& previousNode){
  for (list<s32>::iterator j = comp.begin(); j != comp.end(); ++j) {
    distance.insert(pair<s32,s32>(*j,0));
    previousNode.insert(pair<s32,s32>(*j,-1));
  }
  distance[startNode] = 0;
}

void NetEx::updateDist(s32 s1, s32 s2,map<s32,s32>& distance, map<s32,s32>& previousNode){
  edge_t e =  edge(s1,s2,_graph).first;
  s32 weight=	get(edge_weight_t(), _graph,e);
  if( distance[s2] < distance[s1] + weight){
    distance[s2] = distance[s1] + weight;
    previousNode[s2] = s1;
  }
}

s32 NetEx::findMinDist(list<s32> Q, map<s32,s32>distance){
  s32 min = -1;
  s32 saveNode = -1;
  
  for(list<s32>::iterator itQ = Q.begin(); itQ != Q.end(); ++itQ){
    if(min ==-1 || distance[*itQ] < min ){
      min = distance[*itQ];
      saveNode = *itQ;
    }
  }
  return saveNode;
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

void NetEx::bfs(s32 startId, map<s32,list < list<s32> > >& paths, s32 currentColor,
		set<s32>& colorNotAllowed, std::queue<s32> & Qbfs, set<s32> visited){
  
  out_edge_it ei, ei_end;
  time_t before, after, beforeFor;
  
  while (! Qbfs.empty()) {
    bool boolRecurciveCall = false;
    s32 parentId = Qbfs.front(); Qbfs.pop();
    if(visited.find(parentId) != visited.end()){
      continue;
    }
    
    if(SSRContig::VERBOSE) cout << "loop on Qbfs size "<< Qbfs.size() << "startId= " << startId << endl;
    map<s32,list<list <s32> > >::iterator itPathParent = paths.find(parentId);
    
    before = time(NULL);
    if(SSRContig::VERBOSE) cout << "in BFS out_degree " << out_degree(parentId, _graph) << " node " << parentId << endl;
    for (tie(ei, ei_end) = out_edges(parentId, _graph); ei != ei_end; ++ei) {
      TSSRList::iterator it1 = _vertices->begin();
      std::advance(it1,target(*ei, _graph));
      s32 childId = (*it1)->getID();
      if(SSRContig::VERBOSE) cout << "child " << *(*it1) << " out degree " << out_degree(childId, _graph)<< endl;
      map<s32,bool> mapIdTranscrit = (*it1)->getIdTranscrit();
      set<s32> tmpMapIdTranscrit;
      if(SSRContig::VERBOSE) cout <<"currentColor " << currentColor << " child Color it as "
				  << mapIdTranscrit.size() <<" color ";
      for(map<s32,bool>::iterator itMapIdTranscrit = mapIdTranscrit.begin(); 
	  itMapIdTranscrit != mapIdTranscrit.end(); ++itMapIdTranscrit){
	tmpMapIdTranscrit.insert(itMapIdTranscrit->first);
	if(SSRContig::VERBOSE) cout << itMapIdTranscrit->first << " ";
      }
      if(SSRContig::VERBOSE) cout << endl;
      
      map<s32,bool>::iterator itColorChild = mapIdTranscrit.find(currentColor);
      
      if(itColorChild == mapIdTranscrit.end()){//the current color is not in the mapIdTranscrit
	//bool boolColorNotAllowed = false;
	beforeFor = time(NULL);
	std::vector<int> v_intersection;
	std::set_intersection(colorNotAllowed.begin(), colorNotAllowed.end(),
			      tmpMapIdTranscrit.begin(), tmpMapIdTranscrit.end(),
			      std::back_inserter(v_intersection));
	if(v_intersection.size() != 0 ){
	  visited.insert(parentId);
	  continue;
	}
	else{
	  cout << " else colorNotAllowed.push_back(currentColor); "<<endl;
	  //Put all colorFrom my current Node in colorNotAllowed
	  TSSRList::iterator it1 = _vertices->begin();
	  std::advance(it1,parentId);
	  SSRContig* ctg = *it1;
	  if(SSRContig::VERBOSE) cout <<"parent " <<  *ctg << endl;
	  map<s32,bool> mapIdTranscritParent = ctg->getIdTranscrit();
	  //mapA.insert(mapB.begin(), mapB.end())
	  //	colorNotAllowed.insert(mapIdTranscrit.begin(),mapIdTranscrit.end());
	  //put colr not Allowed from Parent !
	  //		for(map<s32,bool>::iterator itMapIdTranscritParent = mapIdTranscritParent.begin(); itMapIdTranscritParent != mapIdTranscritParent.end(); ++itMapIdTranscritParent){
	  //			colorNotAllowed.insert(itMapIdTranscritParent->first);
	  //		}
	  colorNotAllowed.insert(currentColor);
	  if(mapIdTranscrit.size() > 1){
	    boolRecurciveCall = true;
	    if(SSRContig::VERBOSE) cout << "size colored not allow " << colorNotAllowed.size()<< endl;
	  }
	}
	after = time(NULL);
	if(SSRContig::VERBOSE) cout << "time allPathFinder for colorNotAllowed " << difftime(after,beforeFor)<< endl;
      }
      map<s32,list<list <s32> > >::iterator itPathChild = paths.find(childId);
      //Start BFS
      if(visited.find(childId) == visited.end()){
	if(itPathChild == paths.end()){//if PathChild doesnt exist
	  //		cout << "itPathChild == paths.end()"<<endl;
	  list<list<s32> > newList;
	  
	  if(itPathParent != paths.end()){//if Parent contains a path
	    beforeFor = time(NULL);
	    if(SSRContig::VERBOSE) cout << "itPathParent->second " << itPathParent->second.size() << endl;
	    for(list<list <s32> >::iterator itListPathParent = itPathParent->second.begin(); itListPathParent != itPathParent->second.end();++itListPathParent){
	      list<s32> newPath;
	      for(list<s32>::iterator itL = itListPathParent->begin() ; itL != itListPathParent->end() ; ++itL){
		newPath.push_back(*itL);
	      }
	      newPath.push_back(childId);
	      newList.push_back(newPath);
	    }
	    after = time(NULL);
	    if(SSRContig::VERBOSE) cout << "time allPathFinder for itPathParent" << difftime(after,beforeFor)<< endl;
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
	  
	  
	  //Qbfs.push(childId);
	  //			cout << "childId " << childId << " push in QBFS " << Qbfs.size()<< endl;
	  beforeFor = time(NULL);
	  list<list <s32> > copyList(paths[parentId]);
	  if(SSRContig::VERBOSE) cout << "copyList " << copyList.size() << endl;
	  for(list<list <s32> >::iterator itCopyList = copyList.begin(); itCopyList != copyList.end();++itCopyList){
	    itCopyList->push_back(childId);
	  }
	  itPathChild->second.merge(copyList);
	  itPathChild->second.unique();
	  after = time(NULL);
	  if(SSRContig::VERBOSE) cout << "time allPathFinder for itCopyList" << difftime(after,beforeFor)<< endl;
	  
	  
	}
	Qbfs.push(childId);
      }
      
      visited.insert(parentId);
      if(boolRecurciveCall && out_degree(target(*ei, _graph),_graph) > 0 && mapIdTranscrit.size()> 1){
	if(SSRContig::VERBOSE) cout << "recurcive call searchColor " << endl;
	searchAllColor(target(*ei, _graph), paths,colorNotAllowed,visited);
      }
    }
    
    if(out_degree(parentId,_graph)== 0){
      
      TSSRList::iterator it1 = _vertices->begin();
      std::advance(it1,parentId);
      //	cout << "in bfs node as no out_edges : " << **it1 << endl;
      //	for (list < list<s32> >::iterator it = paths[idNodeParent].begin(); it != paths[idNodeParent].end(); ++it){
      //		cout <<"path of no edges ";
      //		for(list<s32>::iterator it2 = it->begin(); it2 != it->end(); ++it2){
      //			cout<<  *it2 << " ";
      //		}
      cout << endl;
      //
      //	}
    }
    
    after = time(NULL);
    
    if(SSRContig::VERBOSE) printf ("time allPathFinder BFS %.f \n",difftime(after,before) );
  }//While
  //	cout << "exit BFS "<< endl;
}


/*
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
*/

void NetEx::printAllNodes(list<s32> comp){
  for (list<s32>::iterator j = comp.begin(); j != comp.end() ; j++) {
    cout << "j " << *j << endl;
    TSSRList::iterator it1 = _vertices->begin();
    std::advance(it1,*j);
    SSRContig* ctg = *it1;
    cout << *ctg << " color ";
    
    map<s32,bool> mapIdTranscrit = ctg->getIdTranscrit();
    for(map<s32,bool>::iterator itMapIdTranscrit = mapIdTranscrit.begin(); itMapIdTranscrit != mapIdTranscrit.end(); ++itMapIdTranscrit){
      cout << itMapIdTranscrit->first << " ";
    }
    cout << endl;
  }
  //		exit(1); 
}


list<list<s32> > NetEx::PathsFinderWithCondition(list<s32> comp){
  list<s32> sources = sourcesNodes(comp);
  map<s32,list < list<s32> > >predM;
  list<s32> endList = endNodes(comp);
  //	list<s32> endListVertex;
  //	for(list<s32>::iterator itEnd = endList.begin() ; itEnd != endList.end() ; ++itEnd){
  //		s32 childId = nodeToVertex(*itEnd);
  //		endListVertex.push_back(childId);
  //	}
  
  list<list<s32> > paths;
  time_t before, after;
  //	int tmpCmp = 0;
  //	map<s32,s32> vertexToBGL = mapIdsVertextoBGL();
  for(list<s32>::iterator itSources = sources.begin() ; itSources != sources.end() ; ++itSources){
    if(SSRContig::VERBOSE) cout << "PathsFinderWithCondition number of sources nodes " << sources.size()<< " size cc " << comp.size()<< " source "<< *itSources << endl;
    predM.clear();
    set<s32> colorNotAllowed;
    set<s32> visited;
    before = time(NULL);
    searchAllColor(*itSources, predM,colorNotAllowed,visited);
    after = time(NULL);
    if(SSRContig::VERBOSE) cout << " time allPathFinder searchAllColor " << difftime(after,before)<< endl;
    //		printMap(predM);
    before = time(NULL);
    for(list<s32>::iterator itEnd = endList.begin() ; itEnd != endList.end() ; ++itEnd){
      //s32 childId = nodeToVertex(*itEnd);
      map<s32,list < list<s32> > >::iterator itPred = predM.find(*itEnd);
      paths.splice(paths.end(), predM[*itEnd]);
    }
    after = time(NULL);
    if(SSRContig::VERBOSE) 	cout << " time allPathFinder for loop " << difftime(after,before)<< endl;
    //time_t  t = time(0);
    //cerr << ctime(&t)<<endl;
    //	if (tmpCmp > 4)
    //		exit(1);
    //	++tmpCmp;
  }
  paths.sort();
  paths.unique();
  return paths;
}

void NetEx::searchAllColor(s32 idSource,map<s32,list < list<s32> > >& predM,set<s32>& colorNotAllowed, set<s32>& visited){
  //time_t before1 = time(NULL);
  TSSRList::iterator it1 = _vertices->begin();
  std::advance(it1,idSource);
  SSRContig* ctg = *it1;
  
  if(SSRContig::VERBOSE) cout <<"parent " <<  *ctg << endl;
  map<s32,bool> mapIdTranscrit = ctg->getIdTranscrit();
  //time_t before = time(NULL);
  map<s32,list<list <s32> > >::iterator itPathParent = predM.find(idSource);
  if(out_degree(idSource,_graph) == 0 && itPathParent == predM.end()){//for mono exoniques
    //	cout << "enter in if mono exonique " << endl;
    list<s32> newPath;
    list<list<s32> > newList;
    newPath.push_back(idSource);
    newList.push_back(newPath);
    predM.insert(make_pair(idSource,newList));
  }
  else if(out_degree(idSource,_graph) != 0 ){
    for(map<s32,bool>::iterator itMap = mapIdTranscrit.begin() ; itMap != mapIdTranscrit.end(); ++itMap){
      //colorNotAllowed.clear();
      //	cout << "\nloop in searchColor " << endl;
      std::queue<s32> Qbfs ;
      s32 currentColor = itMap->first;
      set<s32>::iterator it;
      it = find (colorNotAllowed.begin(), colorNotAllowed.end(), currentColor);
      //		  if (it == colorNotAllowed.end()){
      if(visited.find(idSource) == visited.end())
	Qbfs.push(idSource);
      bfs(idSource,predM , currentColor,colorNotAllowed,Qbfs,visited);
      //		  }
    }
  }
  //	else{
  //	 cout << "out degree == 0 end path " << *ctg << endl;
  //	}
  //time_t after = time(NULL);
}


// NOUVEAU allpath - to compute paths a lot faster than before
void NetEx::compute_paths(s32 node, list<s32>& path, map<s32, list<list<s32> > >& solved) {
  list<s32> tmp_path = path; 
  // On est arrivé sur un noeud déjà connu ?
  // s32 idVertex = nodeToVertex(node);
  if (solved.find(node) != solved.end()) {
    //  if (solved.find(idVertex) != solved.end()) {
    list<s32> final_path;
    list<s32> *path_ref;
    s32 start_node;
    list<list<s32> >::iterator solved_path;
    // On ajoute tous les sous chemins à la liste résolue
    while(tmp_path.begin() != tmp_path.end()) {
      start_node = tmp_path.front();
      // On concatène le chemin courant avec la solution trouvée
      //    for (solved_path=solved[idVertex].begin(); solved_path != solved[idVertex].end(); ++solved_path) {
      for (solved_path=solved[node].begin(); solved_path != solved[node].end(); ++solved_path) {
	solved[start_node].push_back(tmp_path);
	path_ref = &solved[start_node].back();
	path_ref->insert(path_ref->end(), solved_path->begin(), solved_path->end());
      }
      tmp_path.pop_front();
    }
  }
  else {
    // Ce qui se passe quand on est pas sur un noeud déjà connu
    tmp_path.push_back(node);
    // tmp_path.push_back(idVertex);
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
    //  cout << "cc " << component[i] << " i "<<i <<  endl;
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
    //	cout << "tagExon continue " << **itComponent << " in degree " <<in_degree((*itComponent)->getID(), _graph) << " out degree "<< out_degree((*itComponent)->getID(), _graph) <<endl;

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

void NetEx::updateWeightStartNode(s32 nodeId,s32 formerNode ){
  for (pair<out_edge_it, out_edge_it> pEdges = out_edges(nodeId, _graph); pEdges.first != pEdges.second; ++pEdges.first){
    s32 outEdge = target(*pEdges.first, _graph);
    edge_t e; bool b;
    tie(e,b) = add_edge(formerNode,outEdge, _graph); //For the moment id in vertices and graph are equal
    s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
    formerWeight += get(edge_weight_t(), _graph,e);
    put(edge_weight_t(), _graph, e, formerWeight);
  }
}

void NetEx::updateWeightEndNode(s32 nodeId, s32 newNode){
  for (pair<in_edge_it, in_edge_it> pEdges = in_edges(nodeId, _graph); pEdges.first != pEdges.second; ++pEdges.first){
    s32 inEdge = source(*pEdges.first, _graph);
    edge_t e; bool b;
    tie(e,b) = add_edge(inEdge,newNode, _graph);
    s32 formerWeight =  get(edge_weight_t(), _graph,*pEdges.first);
    formerWeight += get(edge_weight_t(), _graph,e);
    put(edge_weight_t(), _graph, e, formerWeight);
  }
}

void NetEx::cleanGraph(){
  cerr << "Before cleaning, number vertices " << num_vertices(_graph) << " number edges "<< num_edges(_graph)<<endl;
  //it will clean the component regarding the start and end exons.
  s32 left,right,left2,right2;
  //Tag exons
  int deleteMono = 0;
  TSSRList listMono; // OutExon = exons at begin or end of the graph, inExon = exons intern of the graph
  map<s32,TSSRList> mapStartExon, mapEndExon;
  map<string,s32> inExon;
  tagExons(inExon,mapStartExon,mapEndExon,listMono);
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
	    updateWeightStartNode(ctgNext->getID(),ctgStart->getID());
	    clear_vertex(ctgNext->getID(),_graph); //clear all edges before removing the vertex
	    _graph[ctgNext->getID()].guid.deleteVertex = true;
	    //			cout << "clear vertex " << *ctgNext << endl;
	  }
	  else {//we delete the current
	    updateWeightStartNode(ctgStart->getID(),ctgNext->getID());
	    clear_vertex(ctgStart->getID(),_graph); //clear all edges before removing the vertex
	    _graph[ctgStart->getID()].guid.deleteVertex = true;
	    //			cout << "clear node " << *ctgStart << endl;
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
      //		cout << "start continue " << *ctgStart << " in degree " <<in_degree(ctgStart->getID(), _graph) << " out degree "<< out_degree(ctgStart->getID(), _graph) <<endl;
      
    }
  }
  // simplification pour endExon
  for(map<s32,TSSRList>::iterator itMap = mapEndExon.begin(); itMap != mapEndExon.end(); ++itMap){
    for(TSSRList::iterator itCurrent = itMap->second.begin() ; itCurrent != itMap->second.end();++itCurrent){
      SSRContig* ctgEnd = *itCurrent;
      if ( _graph[ctgEnd->getID()].guid.deleteVertex == true) continue;
      for(TSSRList::iterator itNext = ++itCurrent ; itNext  != itMap->second.end();++itNext){
	SSRContig* ctgNext = *itNext;
	//	cout << "Change Next " << endl;
	if(ctgEnd->strand() != ctgNext->strand()) continue;
	leftRight(ctgEnd,left,right);
	leftRight(ctgNext,left2,right2);
	if(overlap(ctgEnd, ctgNext) && ctgEnd->start() == ctgNext->start()){
	  if(ctgEnd->size() > ctgNext->size()){ // we keep current and delete next
	    updateWeightEndNode(ctgNext->getID(), ctgEnd->getID());
	    clear_vertex(ctgNext->getID(),_graph); //clear all edges before removing the vertex
	    _graph[ctgNext->getID()].guid.deleteVertex = true;
	    //		cout << "clear Node end " << *ctgNext << endl;
	    
	  }
	  else{ //if (deleteList.find(itCurrent) == deleteList.end()){//we delete the current
	    updateWeightEndNode(ctgEnd->getID(),ctgNext->getID());
	    clear_vertex(ctgEnd->getID(),_graph); //clear all edges before removing the vertex
	    _graph[ctgEnd->getID()].guid.deleteVertex = true;
	    //			cout << "clear node end else " << *ctgEnd << endl;
	  }
	}
      }
      
      --itCurrent;
      if(_graph[ctgEnd->getID()].guid.deleteVertex == true)
	continue;
      //Fusionne mono FIXME Good idea to do it here ? FIXME Pourquoi on fusionne 2 fois les monos ?
      for(TSSRList::iterator itMono = listMono.begin(); itMono != listMono.end(); ++itMono){
	if((*itMono)->strand() != ctgEnd->strand())
	  continue;
	if((*itMono)->start() >= ctgEnd->start() && (*itMono)->end()>= ctgEnd->end() && overlap(ctgEnd,*itMono)){
	  ctgEnd->setEnd((*itMono)->end());
	  _graph[(*itMono)->getID()].guid.deleteVertex = true;
	  ++deleteMono;
	  //		cout << "clear mono 1 " << **itMono << endl;
	}
      }
      if(_graph[ctgEnd->getID()].guid.deleteVertex == true)
	continue;
      //fusionne start/end exons
      for(map<s32,TSSRList>::iterator itMapStart = mapStartExon.begin(); itMapStart != mapStartExon.end(); ++itMapStart){
	for(TSSRList::iterator itCurrentStart = itMapStart->second.begin() ; itCurrentStart != itMapStart->second.end();++itCurrentStart){
	  //		cout << "start/end continue " << **itCurrentStart << " in degree " <<in_degree((*itCurrentStart)->getID(), _graph) << " out degree "<< out_degree((*itCurrentStart)->getID(), _graph) <<endl;
	  
	  if(_graph[ (*itCurrentStart)->getID()].guid.deleteVertex == true)
	    continue;
	  if(!overlap(*itCurrent,*itCurrentStart) || (*itCurrent)->strand() != (*itCurrentStart)->strand())
	    continue;
	  if(ctgEnd->start() <= (*itCurrentStart)->start() && ctgEnd->end() <= (*itCurrentStart)->end() &&ctgEnd->end() !=  (*itCurrentStart)->start() ){
	    SSRContig* newCtg = new SSRContig(*ctgEnd);//FIXME memory leak
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
	      updateWeightStartNode((*itCurrentStart)->getID(),newCtg->getID());
	      updateWeightEndNode(ctgEnd->getID(),newCtg->getID());
	      //FIXME test
	      
	      clear_vertex(ctgEnd->getID(),_graph); //clear all edges before removing the vertex
	      _graph[ctgEnd->getID()].guid.deleteVertex = true;
	      //			 cout << "clear vertex " << *ctgEnd <<  " in degree " <<in_degree(ctgEnd->getID(), _graph) << " out degree "<< out_degree(ctgEnd->getID(), _graph) <<endl;
	      
	    }
	    else{
	      delete(newCtg);
	      updateWeightStartNode((*itCurrentStart)->getID(),itKey->second);
	      updateWeightEndNode(ctgEnd->getID(),itKey->second);
	    }
	  }
	  
	}
	
      }
      //	 cout << "end continue " << *ctgEnd << " in degree " <<in_degree(ctgEnd->getID(), _graph) << " out degree "<< out_degree(ctgEnd->getID(), _graph) <<endl;
      
    }
    
  }
  
  //Just for long reads
  //TODO arranger le code, pas ultra efficace car on reparcourir les start et end déjà nettoyé plus haut
  for(TSSRList::iterator itComponent = _vertices->begin(); itComponent != _vertices->end(); ++itComponent){
    for(map<s32,TSSRList>::iterator itMapStart = mapStartExon.begin(); itMapStart != mapStartExon.end(); ++itMapStart){
      for(TSSRList::iterator itCurrentStart = itMapStart->second.begin() ; itCurrentStart != itMapStart->second.end();++itCurrentStart){
	
	if(_graph[ (*itCurrentStart)->getID()].guid.deleteVertex == true)
	  continue;
	if(!in_degree((*itComponent)->getID(), _graph)){
	  //			cout << "start Exon, continue " << **itComponent << "in degree " <<in_degree((*itComponent)->getID(), _graph) << " out degree "<< out_degree((*itComponent)->getID(), _graph) <<endl;
	  continue;
	}
	
	if((*itCurrentStart)->end() == (*itComponent)->end() && (*itCurrentStart)->start() >= (*itComponent)->start() && (*itComponent)->strand() == (*itCurrentStart)->strand()){
	  //delete start exon
	  updateWeightStartNode((*itCurrentStart)->getID(),(*itComponent)->getID());
	  clear_vertex((*itCurrentStart)->getID(),_graph); //clear all edges before removing the vertex
	  _graph[(*itCurrentStart)->getID()].guid.deleteVertex = true;
	  //			cout << "new clean start" << **itCurrentStart << "in degree " <<in_degree((*itCurrentStart)->getID(), _graph) << " out degree "<< out_degree((*itCurrentStart)->getID(), _graph) <<endl;
	  
	}
      }
    }
    for(map<s32,TSSRList>::iterator itMap = mapEndExon.begin(); itMap != mapEndExon.end(); ++itMap){
      for(TSSRList::iterator itEnd = itMap->second.begin() ; itEnd != itMap->second.end();++itEnd){
	if(_graph[ (*itEnd)->getID()].guid.deleteVertex == true)
	  continue;
	if( !out_degree((*itComponent)->getID(), _graph))
	  continue;
	
	if((*itEnd)->start() == (*itComponent)->start() && (*itEnd)->end() <= (*itComponent)->end() && (*itEnd)->strand() == (*itComponent)->strand()){
	  updateWeightEndNode((*itEnd)->getID(),(*itComponent)->getID());
	  clear_vertex((*itEnd)->getID(),_graph); //clear all edges before removing the vertex
	  _graph[(*itEnd)->getID()].guid.deleteVertex = true;
	  //		cout << "new clean end " << **itEnd << "in degree " <<in_degree((*itEnd)->getID(), _graph) << " out degree "<< out_degree((*itEnd)->getID(), _graph) <<endl;
	  
	}
      }
    }
    
  }
  deleteNode();
  cerr << "After cleaning, number vertices " << num_vertices(_graph) << " number edges "<< num_edges(_graph)<<endl;
}

void NetEx::deleteNode(){
  for(u32 vi = 0 ; vi < num_vertices(_graph); ){
    if(_graph[vi].guid.deleteVertex == true){
      TSSRList::iterator it1;
      it1 = _vertices->begin();
      int nb =vi;
      std::advance(it1,nb);
      //		cout << " delete vertex " << **it1 << "in degree " <<in_degree((*it1)->getID(), _graph) << " out degree "<< out_degree((*it1)->getID(), _graph) <<endl;
      
      _vertices->erase (it1);
      clear_vertex(vi,_graph);
      remove_vertex(vi, _graph);
    }
    else
      ++vi;
  }
}

void NetEx::synchronisedId(){ //boost library update ids after deleting a node. We should update our vertices to keep track of them
  for(u32 vi = 0 ; vi < num_vertices(_graph); ++vi){
    TSSRList::iterator it1;
    it1 = _vertices->begin();
    std::advance(it1,vi);
    (*it1)->setID(vi);
  }
}

void NetEx::simplifyBigGraph(list<s32>comp, s32 threshold){
  string seqname;
  map<Edge,s32> weight;
  map<int,int> degreeNode;
  s32 nbEdges = 0; 
  s32 getWeight;

  for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++){
    int degreeBefore = degree(*j, _graph);
    degreeNode.insert(make_pair(*j, degreeBefore));
  }
  
  for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++){
    nbEdges += out_degree(*j, _graph);
    for (pair<out_edge_it, out_edge_it> pEdges = out_edges(*j, _graph); pEdges.first != pEdges.second; ++pEdges.first){
      getWeight =  get(edge_weight_t(), _graph,*pEdges.first);
      if(getWeight < threshold )
	remove_edge(*pEdges.first, _graph);
    }
  }
  for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++){
    int newDegree =  degree(*j, _graph);
    int degreeBefore = degreeNode.find(*j)->second;
    
    if(degreeBefore !=0 && newDegree ==0)
      _graph[*j].guid.deleteVertex = true;
  }
  deleteNode();
}

// get weight of a given path
s32 NetEx::pathWeight(list<s32>& path) {
  s32 weight = 0;
  list<s32>::iterator itA = path.begin();
  for(list<s32>::iterator itB = path.begin() ; itB != path.end() ; itB++){
    if(itA == itB) continue;
    // edge could have been removed in large graphs
    if(edge(*itA, *itB, _graph).second)
      weight += get(edge_weight_t(), _graph, edge(*itA, *itB, _graph).first);
    itA = itB;
  }
  return weight;
}

// get maximum weigth in the CC
s32 NetEx::maxWeigth(list<s32>& comp) {
  s32 max = 0;
  list<s32>::iterator itA = comp.begin();
  for(list<s32>::iterator itB = comp.begin() ; itB != comp.end() ; itB++){
    if(itA == itB) continue;
    if(edge(*itA, *itB, _graph).second) {
      s32 weight = get(edge_weight_t(), _graph, edge(*itA, *itB, _graph).first);
      if(weight > max) { max = weight; }
    }
    itA = itB;
  }
  return max;
}
