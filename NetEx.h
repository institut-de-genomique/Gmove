/*******************************************************************************
+
+  NetEx.h
+
+  
+  
+ 
 *******************************************************************************/


#ifndef NETEX_H
#define NETEX_H

#include "LocalType.h"
#include "SSRContig.h"
#include "DnaDictionary.h"

#include <boost/config.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/astar_search.hpp>

/*
#include <boost/graph/random.hpp>
#include <boost/random.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/filtered_graph.hpp>
*/


#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <algorithm>


using namespace std;
using namespace boost;

struct Guid{ bool deleteVertex;};
struct MyVertex{
	Guid guid;
	SSRContig* vertex;
};




typedef adjacency_list<setS, vecS, bidirectionalS,MyVertex> Graph; //TODO put setS
typedef adjacency_list< vecS, vecS, undirectedS> GraphU;
typedef property_map<Graph, vertex_index_t>::type IndexMap;
typedef graph_traits<Graph>::vertex_iterator vertex_it;
typedef graph_traits<Graph>::edge_iterator edge_it;
typedef graph_traits<Graph>::out_edge_iterator out_edge_it;
typedef graph_traits<Graph>::in_edge_iterator in_edge_it;
typedef graph_traits<Graph>::vertices_size_type size_type;
//typedef graph_traits <Graph>::vertex_property vertex_property;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::edge_descriptor edge_t; // (source,target)
typedef map<Vertex, SSRContig*> IMap;

/*

struct vertex_hash:std::unary_function<Vertex, std::size_t>
{
    vertex_hash(Graph const& g):g(g){}

    std::size_t operator()(Vertex const& v) const {
    return g[v].guid.deleteVertex;
  }

  Graph const& g;
};

typedef boost::unordered_set<Vertex, vertex_hash> vertex_set;
typedef filtered_graph<Graph,keep_all,is_not_in_subset<vertex_set> > FilteredGraph;
*/

typedef size_type* Rank;
typedef Vertex* Parent;

typedef pair<s32,s32> Edge;
typedef pair<s32,SSRContig*> Vertice;
typedef list<SSRContig*> TSSRList;

/*
template <typename GraphType>
std::size_t filtered_num_vertices(const GraphType& g)
{
    std::size_t total=0;
    BGL_FORALL_VERTICES_T(v,g,GraphType){
        ++total;
    }
    return total;
}

template <typename GraphType>
std::size_t filtered_num_edges(const GraphType& g)
{
    std::size_t total=0;
    BGL_FORALL_EDGES_T(e,g,GraphType)
    {
        ++total;
    }
    return total;
}

*/

class NetEx {
  
  TSSRList* _vertices;
  list<Edge>* _edges;
  s32 _strand;
  string _seqname;
  Graph _graph;
  s32 _nb_connected_components;
  map<s32, SSRContig*> VertexMap;
  
 public:
  
  /* Constructors and Destructors*/  
  NetEx(TSSRList* exons, list<Edge>* edges, string &seqname) {
    _vertices = exons;
    _edges = edges;
    _seqname = seqname;
    _graph = Graph(_edges->begin(), _edges->end(), _vertices->size());

   // map<int,Vertex >Vmap;
  /*  IMap idxMap;
    associative_property_map<IndexMap> iMap(idxMap);
    int i = 0 ;
*/

    vertex_it vi, vi_end, next;
    	tie(vi, vi_end) = vertices(_graph);


    	for (next = vi; vi != vi_end; vi = next) {
    		TSSRList::iterator it1;
			it1 = _vertices->begin();
			int nb = *vi;
			advance(it1,nb);
    		++next;
    //		cout << num_vertices(_graph) << " numVertices " << *vi << " *vi" << endl;
    		_graph[*vi].guid.deleteVertex = false;
    		_graph[*vi].vertex = *it1;
    	}


    for(TSSRList::iterator it = exons->begin();it != exons->end();++it){

    	VertexMap.insert(make_pair((*it)->getID(),*it));
   // 	cout << " VertexMap " << (*it)->getID() << " " << *(*it) << endl;
    }

  }
  
  ~NetEx () {
	  delete _edges;
	  delete _vertices;
  }
  
  
  /* Accessors */
  TSSRList* getVertices() const { return _vertices; }
  list<Edge>* getEdges() const { return _edges; }
  string getSeqName() const { return _seqname; }
  Graph* getGraph() { return &_graph; }
  s32 getNb_cc() { return _nb_connected_components; }
  
  /* Methods */
  vector<list<s32> > getComponents();

  list<list<int> > allPathsFinder(list<s32>&);
  void countVerticesAndEgdes(list<s32>&, s32&, s32&);
  s32 count_allPaths(list<s32>&);
  s32 count_paths(s32, map<s32, s32>&);
  void compute_paths(s32, list<s32>&, map<s32, list<list<s32> > >&);
  void graphOut();
  s32 nbCycle();
  void tagExons(map<string,s32>&,map<s32,TSSRList>&,map<s32,TSSRList>&,TSSRList& );
  void cleanGraph();
  void deleteNode();
};

//Finding out cycles
struct cycle_detector : public dfs_visitor<> {  
  cycle_detector( bool& has_cycle, int& num):_has_cycle(has_cycle),cycle(num) { }
  
  template <class Edge, class Graph> void back_edge(Edge, Graph&) {
    cycle++;      
    _has_cycle = true;
  }
  
  protected:
  bool& _has_cycle;
  int& cycle;
}; 

#endif

