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
void NetEx::graphOut (){
	/*
	 * // Graph structure with customized property output
		template < typename VertexAndEdgeListGraph, typename VertexPropertyWriter >
	void write_graphviz(std::ostream& out, const VertexAndEdgeListGraph& g, VertexPropertyWriter vpw);
	 *
	 * typedef boost::adjacency_list<boost::listS, boost::vecS, boost::directedS,
    boost::property<boost::vertex_name_t, unsigned int,
    boost::property<boost::vertex_index_t, unsigned int> >,
    boost::property<boost::edge_name_t, unsigned int> > Graph;

  *  typedef boost::property_map<Graph, boost::vertex_name_t>::type VertexNameMap;
  *  VertexNameMap vname_map1 = get(boost::vertex_name, graph1);
  * write_graphviz(file_graph1, graph1,
                     make_label_writer(vname_map1),
                     make_label_writer(ename_map1));
	 */
  string filename = _seqname;
  stringstream ss;
  ss << "Model.dot";
  filename += ss.str();
  filename = "dot_files/"+filename;
  ofstream ofs( filename.c_str() );
  int i = 0;
  list<string> listLabel;
  for(TSSRList::iterator it = _vertices->begin(); it != _vertices->end(); ++it) {

	SSRContig* ctg = *it;
	stringstream tmpLabel;
	tmpLabel<< ctg->getID()<<" "<< ctg->start() <<"-" << ctg->end();
	string tmpLabel2 = tmpLabel.str();//convert( ctg->start(),ctg->end());//tmpLabel.str();
	listLabel.push_back(tmpLabel2);
	++i;
   }
  char* name;
  for(std::list<string>::const_iterator it = listLabel.begin(); it != listLabel.end(); ++it){
	  strcat(name, it->c_str());
  }
 write_graphviz(ofs, _graph,make_label_writer((const char *)name));
}

// to count the number of vertices (exons) and edges (introns) in the graph
void NetEx::countVerticesAndEgdes(list<s32>& comp, s32& nbVertices, s32& nbEdges) {
  for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++) {
    nbVertices++;
    for (pair<out_edge_it, out_edge_it> p = out_edges(*j, _graph) ; p.first != p.second ; ++p.first ) nbEdges++;
  }
}

// to count the number of paths in a connected component
s32 NetEx::count_allPaths(list<s32>& comp) {
  list<s32> sources;
  map<s32, s32> subs;
  
  for (list<s32>::iterator j = comp.begin(); j != comp.end(); ++j) {
    if (!out_degree(*j, _graph))
      subs[*j] = 1;
    
    if(!in_degree(*j, _graph)) {
    	sources.push_back(*j);
    }
  }

  s32 total = 0;
  for (list<s32>::iterator itm = sources.begin() ; itm != sources.end(); itm++) { 
    s32 tmp = count_paths(*itm, subs);
    if(tmp!=-1) total += tmp;
    else {return -1;}
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

/* Ancienne version
s32 NetEx::count_allPaths(list<s32>& comp) {
  list<s32> sources;
  map<s32,s32> targets;
  s32 count = 0;

  for (list<s32>::iterator j = comp.begin(); j != comp.end(); j++) {
    if(!out_degree(*j, _graph)) targets[*j]++;
    if(!in_degree(*j, _graph)) sources.push_back(*j);
  }
  
  for (list<s32>::iterator itm = sources.begin() ; itm != sources.end(); itm++) { 
    if(count == 0) count++;
    count_paths(*itm, count, targets);
  }
  return count;
}

void NetEx::count_paths(s32 root, s32& count, map<s32,s32>& ends) {
  for (pair<out_edge_it, out_edge_it> p = out_edges(root, _graph) ; p.first != p.second ; ++p.first ){
    s32 vtarget = target(*p.first, _graph);
    if( ends.find(vtarget) != ends.end() ) count++;
    count_paths(vtarget, count, ends);
  }
}
*/

// all paths finder 
list<list<s32> > NetEx::allPathsFinder(list<s32>& comp) {
//	cout << " allPathsFinder "<< endl;
  list<list<s32> > selectedPaths;
  list<s32> sources;
  // AJOUT MAX
  map<s32, list<list<s32> > > solved;
  list<s32> tlist(1, 0);
  for (list<s32>::iterator j = comp.begin(); j != comp.end(); ++j) {

    if(!out_degree(*j, _graph)) {
      // AJOUT MAX
      // Ajout des noeuds terminaux comme solutions
    	TSSRList::iterator it1 = _vertices->begin();
    	       advance(it1,*j);
    	        int idVertex = (*it1)->getID();
    	   //     cout << "j " << *j << *(*it1) << endl;
      tlist.front() = idVertex;
       solved[idVertex].push_back(tlist);
    }
    if(!in_degree(*j, _graph)) {
    	 sources.push_back(*j);
    //	 cout << " sources j " << *j << endl;
    }

  }

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
		TSSRList::iterator it1 = _vertices->begin();
	    advance(it1,*it);
	//    cout << " path sources " << *it << " " << *(*it1) << endl;
    selectedPaths.splice(selectedPaths.end(), solved[(*it1)->getID()]);
  }
  return selectedPaths;
}

// NOUVEAU allpath - to compute paths a lot faster than before
void NetEx::compute_paths(s32 node, list<s32>& path, map<s32, list<list<s32> > >& solved) {
  list<s32> tmp_path = path; 
  // On est arrivé sur un noeud déjà connu ?
  TSSRList::iterator it1 = _vertices->begin();
       advance(it1,node);
        int idVertex = (*it1)->getID();
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
/*
void NetEx::all_paths(s32 root, list<s32>& path, list<list<s32> >& path_list, map<s32,s32>& ends) {
  list <s32> loc_list;
  for (pair<out_edge_it, out_edge_it> p = out_edges(root, _graph) ; p.first != p.second ; ++p.first ){
    loc_list = path;
    s32 vtarget = target(*p.first, _graph);
    loc_list.push_back(vtarget);
    if( ends.find(vtarget) != ends.end() ) path_list.push_back(loc_list);
    all_paths(vtarget, loc_list, path_list, ends);
  }
}
*/

// to retrieve all connected components of the graph
vector<list<s32> > NetEx::getComponents() {

  GraphU gU(_vertices->size());
  pair<edge_it, edge_it> p = edges(_graph);
  for(edge_it it = p.first ; it != p.second ; it++) {
	  add_edge(source(*it, _graph), target(*it, _graph), gU);
	//  cout << " add edge *it " <<source(*it, _graph) << " " <<  target(*it, _graph)<< endl;
  }

  std::vector<int> component(num_vertices(gU));
  _nb_connected_components = connected_components(gU, &component[0]);
  vector<list<s32> > cc(_nb_connected_components);
  for (s32 i=0 ; i < (s32)component.size() ; i++){
	//  TSSRList::iterator it1 = _vertices->begin();
	 // advance(it1,i);
	  //cc[component[i]].push_back((*it1)->getID());
//	  cout << " i " << i << endl;
 	  cc[component[i]].push_back(i);
  }
//  cout << " get component "<< cc.size() << endl;
  return cc;
}


//to find out whether the graph has cycles..
s32 NetEx::nbCycle() {
//	cout << " enter in nbCycle "<<endl;
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
		else if(!in_degree(ctg->getID(), _graph)){
			if((*itComponent)->tag() != "d"){
				//cerr << "\nerror tagExons start is " << (*itComponent)->tag()  << endl;
				//exit(1);
			}

		//	startExon.push_back(*itComponent);
	//		if((*itComponent)->strand()== SSRContig::FORWARD) mapStartExon[(*itComponent)->end()].push_back(*itComponent);
	//		else endExon.push_back(*itComponent);
			mapStartExon[(*itComponent)->end()].push_back(*itComponent);
		}

		else if(!out_degree(ctg->getID(), _graph)){
			if((*itComponent)->tag() != "f"){
						//	cerr << "\nerror tagExons end is " << (*itComponent)->tag()<< " "<< *(*itComponent)  << endl;
						//	exit(1);
						}
		//	if((*itComponent)->strand()== SSRContig::FORWARD) endExon.push_back(*itComponent);
		//	else mapStartExon[(*itComponent)->end()].push_back(*itComponent);
			mapEndExon[(*itComponent)->start()].push_back(*itComponent);
		}
		else{
			if((*itComponent)->tag() != "i"){
						//	cerr << "\n error tagExons in is " << (*itComponent)->tag()  << endl;
						//	exit(1);
						}
			 SSRContig* exon = *itComponent;
			 ostringstream oss;
			 oss << exon->seqName() << "@" << exon->start() << "@" << exon->end() << "@" << exon->strand();
			 string key = oss.str();
			 inExon.insert(make_pair(key,exon->getID()));
		}
	}
//	cout << " ***** number of mono " << monoExon.size() << endl;
}


void NetEx::cleanGraph(){
	 cerr << " Before cleaning, number vertices " << num_vertices(_graph) << " number edges "<< num_edges(_graph)<<endl;
	//it will clean the component regarding the start and end exons.
	s32 left,right,left2,right2;
	//list<s32> deleteList;
	//Tag exons
	int deleteMono = 0;
	TSSRList listMono; // OutExon = exons at begin or end of the graph, inExon = exons intern of the graph
//	TSSRList startExon,endExon;
	map<s32,TSSRList> mapStartExon, mapEndExon;
	map<string,s32> inExon;
	tagExons(inExon,mapStartExon,mapEndExon,listMono); //TODO change tag Exon
//	cout << " number of mono " << listMono.size() << endl;
	//TODO sort exon by pos and length !
//	s32 newLeft ;

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
						}
						clear_vertex(ctgNext->getID(),_graph); //clear all edges before removing the vertex
						_graph[ctgNext->getID()].guid.deleteVertex = true;
					//	deleteList.push_back(ctgNext->getID());
					}
					else {//we delete the current
						for (pair<out_edge_it, out_edge_it> pEdges = out_edges(ctgStart->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
							s32 outEdge = target(*pEdges.first, _graph);
							edge_t e; bool b;
							tie(e,b) = add_edge(ctgNext->getID(),outEdge, _graph);
						}
						clear_vertex(ctgStart->getID(),_graph); //clear all edges before removing the vertex
						_graph[ctgStart->getID()].guid.deleteVertex = true;

						//deleteList.push_back(ctgStart->getID());
					}
				}
			}
			--itCurrent;
			if(_graph[ctgStart->getID()].guid.deleteVertex == true) continue;
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
	// simplification pour endExon
	for(map<s32,TSSRList>::iterator itMap = mapEndExon.begin(); itMap != mapEndExon.end(); ++itMap){
		for(TSSRList::iterator itCurrent = itMap->second.begin() ; itCurrent != itMap->second.end();++itCurrent){
			SSRContig* ctgStart = *itCurrent;
			if ( _graph[ctgStart->getID()].guid.deleteVertex == true) continue;
			for(TSSRList::iterator itNext = ++itCurrent ; itNext  != itMap->second.end();++itNext){
				SSRContig* ctgNext = *itNext;
				if(ctgStart->strand() != ctgNext->strand()) continue;
				leftRight(ctgStart,left,right);
				leftRight(ctgNext,left2,right2);
				if(overlap(ctgStart, ctgNext) && ctgStart->start() == ctgNext->start()){
					if(ctgStart->size() > ctgNext->size()){ // we keep current and delete next
						for (pair<in_edge_it, in_edge_it> pEdges = in_edges(ctgNext->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
							s32 inEdge = source(*pEdges.first, _graph);
							edge_t e; bool b;
							tie(e,b) = add_edge(inEdge,ctgStart->getID(), _graph);
						}
					 clear_vertex(ctgNext->getID(),_graph); //clear all edges before removing the vertex
					 _graph[ctgNext->getID()].guid.deleteVertex = true;
			//		 cout << "delete "<<*ctgNext <<endl;
				//	 deleteList.push_back(ctgNext->getID());
					}
					else{ //if (deleteList.find(itCurrent) == deleteList.end()){//we delete the current
						for (pair<in_edge_it, in_edge_it> pEdges = in_edges(ctgStart->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
							s32 inEdge = source(*pEdges.first, _graph);
							edge_t e; bool b;
							tie(e,b) = add_edge(inEdge,ctgNext->getID(), _graph);
						}
				//		deleteList.push_back(ctgStart->getID());
						clear_vertex(ctgStart->getID(),_graph); //clear all edges before removing the vertex
						_graph[ctgStart->getID()].guid.deleteVertex = true;
					}
			 }
		 }
			--itCurrent;
			if(_graph[ctgStart->getID()].guid.deleteVertex == true) continue;
			for(TSSRList::iterator itMono = listMono.begin(); itMono != listMono.end(); ++itMono){
				if((*itMono)->strand() != ctgStart->strand()) continue;
							if((*itMono)->start() >= ctgStart->start() && (*itMono)->end()>= ctgStart->end() && overlap(ctgStart,*itMono)){
							ctgStart->setEnd((*itMono)->end());
								_graph[(*itMono)->getID()].guid.deleteVertex = true;
								 ++deleteMono;
							}
						}

		 if(_graph[ctgStart->getID()].guid.deleteVertex == true){
	//		 cout << "deleteVertex = True " << *ctgStart << endl;
			 continue;
		 }

		 //fusionne stat/end exons
		 for(map<s32,TSSRList>::iterator itMapStart = mapStartExon.begin(); itMapStart != mapStartExon.end(); ++itMapStart){
			 for(TSSRList::iterator itCurrentStart = itMapStart->second.begin() ; itCurrentStart != itMapStart->second.end();++itCurrentStart){
				 if(_graph[ (*itCurrentStart)->getID()].guid.deleteVertex == true){
		// 				 cout << "deleteVertex = True " << *(*itCurrentStart) << endl;
					 continue;
				 }
		 		 if(!overlap(*itCurrent,*itCurrentStart) || (*itCurrent)->strand() != (*itCurrentStart)->strand()) continue;
	 			 if(ctgStart->start() <= (*itCurrentStart)->start() && ctgStart->end() <= (*itCurrentStart)->end() &&ctgStart->end() !=  (*itCurrentStart)->start() ){
		 		//		 cout << "*********** We have to clean that "<< ctgStart->start() << "  " << ctgStart->end() << " " <<  (*itCurrentStart)->start()<< " " <<(*itCurrentStart)->end() << endl;
		 		//		cout << "vertices size " << _vertices->size() << " num vertices in graph " << num_vertices(_graph);
	 				 Tstrand strandT = ctgStart->strand();
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
		 		//		 cout << " vertices size " << _vertices->size() << " num vertices in graph " << num_vertices(_graph)<< endl;
		 				 for (pair<out_edge_it, out_edge_it> pEdges = out_edges((*itCurrentStart)->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
		 					 s32 outEdge = target(*pEdges.first, _graph);
		 					 edge_t e; bool b;
		 					 tie(e,b) = add_edge(newCtg->getID(),outEdge, _graph); //For the moment id in vertices and graph are equal
		 		//							 cout << " we put the edge " << e << endl;
		 				 }
		 				 for (pair<in_edge_it, in_edge_it> pEdges = in_edges(ctgStart->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
		 					 s32 inEdge = source(*pEdges.first, _graph);
		 					 edge_t e; bool b;
		 					 tie(e,b) = add_edge(inEdge,newCtg->getID(), _graph);
		 		//								 cout << " we put the edge " << e  << endl;
		 				}
		 			}
		 			else{
		 				for (pair<out_edge_it, out_edge_it> pEdges = out_edges((*itCurrentStart)->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
		 					s32 outEdge = target(*pEdges.first, _graph);
		 					edge_t e; bool b;
		 					tie(e,b) = add_edge(itKey->second,outEdge, _graph); //For the moment id in vertices and graph are equal
		 				//			 									 cout << " we put the edge " << e << endl;
		 				}
		 				for (pair<in_edge_it, in_edge_it> pEdges = in_edges(ctgStart->getID(), _graph); pEdges.first != pEdges.second; ++pEdges.first){
		 					s32 inEdge = source(*pEdges.first, _graph);
		 					edge_t e; bool b;
		 					tie(e,b) = add_edge(inEdge,itKey->second, _graph);
		 				//			 										 cout << " we put the edge " << e  << endl;
		 				}
		 			}
		 		}
			 }
		 }
		}
	}
	 /*
	 	 vertex_hash hasher(_graph);
	 	 vertex_set removed(0,hasher);
	 	 FilteredGraph fg(_graph,boost::keep_all(),removed);

	 	    BGL_FORALL_VERTICES(_vertices,_graph,Graph)
	 	    {
	 	        if (_graph[_vertices].guid.deleteVertex) {
	 	            removed.insert(_vertices);
	 	        }
	 	    }
	 	    std::cout << "After: " << filtered_num_vertices(fg) << " vertices and " << filtered_num_edges(fg) << " edges\n"; //num_vertices and num_edges don't work for filtered_graphs
	 */
	 deleteNode();

	 cerr << " After cleaning, number vertices " << num_vertices(_graph) << " number edges "<< num_edges(_graph)<<endl;
	 cerr << "delete mono " << deleteMono<<endl;
//	 for (TSSRList::iterator it = _vertices->begin(); it != _vertices->end() ; ++it){
//		 cout << "vertex after clean " << *(*it) << endl;
//	 }
}

void NetEx::deleteNode(){
	time_t before = time(NULL);
//	cout << " enter in deleteNode " << endl;

/*	for(TSSRList::iterator itVertex = _vertices->begin() ; itVertex != _vertices->end();){
		cout << " for loop on "  << *(*itVertex) << endl;
		if(find(deleteList.begin(), deleteList.end(),(*itVertex)->getID())!= deleteList.end() ){

			itVertex = _vertices->erase(itVertex);
			cout << " delete "<< _vertices->size()<< endl;
		}
		else
			++itVertex;
	}
*/
//	vertex_it vi, vi_end, next;
//	tie(vi, vi_end) = vertices(_graph);

	s32 numEdge = num_edges(_graph);
//	for (next = vi; vi != vi_end; vi = next) {
/*	someVector.erase(std::remove_if(someVector.begin(), someVector.end(),
	                                [](decltype(*someVector.begin()) element){
	                                    return !element.get()->update();
	                                },
	                 someVector.end());
	*/
	for(int vi = 0 ; vi < num_vertices(_graph); ){
//	cout << vi << endl;
	//	cout << " _graph[*vi].guid.deleteVertex "<< _graph[vi].guid.deleteVertex << endl;
	//	++next;
		time_t newBefore = time(NULL);
		if(_graph[vi].guid.deleteVertex == true){
//			cout << "delete "<< vi << " size vertices "<<num_vertices(_graph)<<" " << _vertices->size()<< endl;
		//
			TSSRList::iterator it1;

			it1 = _vertices->begin();
			int nb =vi;
			advance(it1,nb);
	//		cout << "delete " << *(*it1) <<endl;


			_vertices->erase (it1);
		//	clear_vertex(vi,_graph);
			clear_vertex(vi,_graph);
			remove_vertex(vi, _graph);
			time_t newAfter = time(NULL);
	//		cout << "time delete in for " << difftime(newAfter,newBefore)<<endl;
		}
		else
			++vi;
	}
	time_t after = time(NULL);
	//if(PRINTTIME) cout << " time deleteNode CleanGraph " << difftime(after,before)<<endl;
}
