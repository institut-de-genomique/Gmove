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

#include <boost/config.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/astar_search.hpp>
#include <boost/graph/filtered_graph.hpp>

#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_utility.hpp>

//include for BFS
#include<boost/graph/breadth_first_search.hpp>
#include<boost/graph/visitors.hpp>
//#include <boost/pending/indirect_cmp.hpp>
//#include <boost/pending/integer_range.hpp>
#include <boost/graph/named_function_params.hpp>

#include <boost/config.hpp>
#include <vector>
#include <boost/pending/queue.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/graph/named_function_params.hpp>


#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <list>
#include <queue>
#include <stack>
#include <algorithm>
#include <time.h>

extern bool PRINTTIME;

using namespace std;
using namespace boost;

struct Guid{ bool deleteVertex;};
typedef boost::property<boost::vertex_color_t, boost::default_color_type> ColorPropertyType ;
struct MyVertex{
  Guid guid;
  SSRContig* vertex;
  string name;
  ColorPropertyType color;
};


typedef boost::property<boost::edge_weight_t, int> EdgeWeightProperty;


typedef adjacency_list<setS, vecS, bidirectionalS,MyVertex,EdgeWeightProperty > Graph;//setS vecS
typedef property_map<Graph, edge_weight_t>::type EdgeWeightMap; //for filter graph
typedef adjacency_list< vecS, vecS, undirectedS> GraphU; //vecS vecS
typedef graph_traits<Graph>::vertex_iterator vertex_it;
typedef graph_traits<Graph>::edge_iterator edge_it;
typedef graph_traits<Graph>::out_edge_iterator out_edge_it;
typedef graph_traits<Graph>::in_edge_iterator in_edge_it;
typedef graph_traits<Graph>::vertices_size_type size_type;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef map<Vertex, size_t> IndexMappp;
typedef graph_traits<Graph>::edge_descriptor edge_t; // (source,target)
typedef size_type* Rank;
typedef Vertex* Parent;

typedef pair<s32,s32> Edge;
typedef pair<s32,SSRContig*> Vertice;
typedef list<SSRContig*> TSSRList;

template < typename TimeMap > class bfs_time_visitor:public default_bfs_visitor {
  typedef typename property_traits < TimeMap >::value_type T;
 public:
 bfs_time_visitor(TimeMap tmap, T & t):m_timemap(tmap), m_time(t) { }
  template < typename Vertex, typename Graph >
    void discover_vertex(Vertex u, const Graph & g) const
  {
    put(m_timemap, u, m_time++);
  }
  TimeMap m_timemap;
  T & m_time;
};


class NetEx {
  
  TSSRList* _vertices;
  list<Edge>* _edges;
  map<Edge,s32> _weightEdge;
  string _seqname;
  Graph _graph;
  s32 _nb_connected_components;

 public:
  
  /* Constructors and Destructors*/  
  NetEx(TSSRList* exons, list<Edge>* ed, string &seqname,  map<Edge,s32> w) {
    IndexMappp mapIndexxx;
    associative_property_map<IndexMappp> propmapIndex(mapIndexxx);
    _nb_connected_components = 0 ;
    _vertices = exons;
    _edges = ed;
    _seqname = seqname;
    _graph = Graph(_edges->begin(), _edges->end(), _vertices->size());
    
    _weightEdge = w;
    
    
    pair<edge_it, edge_it> p = edges(_graph);
    for(edge_it it = p.first ; it != p.second ; it++) {
      s32 startEdge, endEdge;
      startEdge = source(*it, _graph);
      endEdge = target(*it, _graph);
      
      Edge pairEdge = make_pair(startEdge,endEdge);
      map<Edge,s32>::iterator itWeightEdge = _weightEdge.find(pairEdge);
      if(itWeightEdge != _weightEdge.end())
	put(edge_weight_t(), _graph, *it, itWeightEdge->second);
    }
    vertex_it vi, vi_end, next;
    tie(vi, vi_end) = vertices(_graph);
    int i=0;
    for (next = vi; vi != vi_end; vi = next) {
      TSSRList::iterator it1;
      it1 = _vertices->begin();
      int nb = *vi;
      std::advance(it1,nb);
      ++next;
      _graph[*vi].guid.deleteVertex = false;
      _graph[*vi].vertex = *it1;
      _graph[*vi].name = (*it1)->getName();
      put(propmapIndex, (*vi), i++);
    }
    
    vertex_it b,e;
  }
  
  ~NetEx () {
    delete _edges;
    delete _vertices;
  }
  
  
  /* Accessors */
  TSSRList* getVertices() const { return _vertices; }
  list<Edge>* getEdges() const { return _edges; }
  string getSeqName() const { return _seqname; }
  s32 getNb_cc() { return _nb_connected_components; }
  
  /* Methods */
  vector<list<s32> > getComponents();

  list<list<int> > allPathsFinder(list<s32>&);

 // s32 nodeToVertex(s32 idNode);
  void bfs(s32 startId,map<s32,list < list<s32> > >& paths , s32 currentColor,set<s32>& colorNotAllowed,
	   std::queue<s32> & Qbfs,set<s32> visited  );
  void dfs(s32 startId,map<s32,list < list<s32> > >& paths , s32,list<s32> &,std::stack<s32> & );
  void printMap (map<s32,list < list<s32> > > mapP);
  void printAllNodes(list<s32> comp); //XXX Just for test
  list<list<s32> > PathsFinderWithCondition(list<s32> comp);
  list<list<s32> > pathWithHigherWeight(list<s32>& comp);
  void initDijkstra(list<s32>comp,s32 startNode,map<s32,s32>& distance,map<s32,s32>&);
  s32 findMinDist(list<s32> Q, map<s32,s32>distance);
  list<s32> shortestPath(s32 stop,s32 start,map<s32,s32> previousPath);
  void updateDist(s32 s1, s32 s2,map<s32,s32>& distance, map<s32,s32>& previousNode);
  void searchAllColor(s32 idSource,map<s32,list < list<s32> > >& predM,set<s32>& colorNotAllowed,set<s32>&);
 // map<s32,s32> mapIdsVertextoBGL();
  void countVerticesAndEgdes(list<s32>&, s32&, s32&);
  list<s32> sourcesNodes(list<s32>comp);
  list<s32> endNodes(list<s32> comp);
  s32 count_allPaths(list<s32>&);
  s32 count_paths(s32, map<s32, s32>&);
  void compute_paths(s32, list<s32>&, map<s32, list<list<s32> > >&);
  void graphOut (string directory, string name);
  s32 nbCycle();
  void tagExons(map<string,s32>&,map<s32,TSSRList>&,map<s32,TSSRList>&,TSSRList& );
  void updateWeightStartNode(s32 , s32 );
  void updateWeightEndNode(s32, s32 );
  void cleanGraph();
  void deleteNode();
  void synchronisedId();
  void simplifyBigGraph(list<s32>, s32);
  s32 pathWeight(list<s32>&);
  s32 maxWeigth(list<s32>&);
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
/*  */
#endif

