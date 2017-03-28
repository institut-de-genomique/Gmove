#include <iostream>
#include <fstream>
#include <cstdio>


#include <boost/graph/adjacency_list.hpp>

#include "DnaDictionary.h"
#include "SSRContig.h"
#include "SSRContigLists.h"
#include "ReadFile.h"


using namespace std;
using namespace boost;

void usage(string msg);
map<string, string> loadFasta(char* filename);

int main(int argc, char** argv)
{
  s32 c, verbose=0, word_size=20, min_nb_overlap_junction=10, min_size_inside_exon=8, extend_for_splicesites=200, min_size_exon=30, nb_neighbour=20, min_size_intron=9, max_size_intron=50000, min_cov_word2orientate = 5; 
  char *readsFilename = NULL, *contigsFilename = NULL, *fastaFilename;

  // Invokes member function `int operator ()(void);'
  while ((c = getopt(argc, argv, "w:r:c:f:h:e:i:m:o:n:x:v")) != -1) {
    switch (c) {  
    case 'r': 
      readsFilename = optarg; 
      break; 
    case 'c': 
      contigsFilename = optarg; 
      break;  
    case 'f': 
      fastaFilename = optarg; 
      break; 
    case 'w': 
      word_size = atoi(optarg); 
      break; 
    case 'e': 
      min_size_exon = atoi(optarg); 
      break; 
    case 'i': 
      min_size_intron = atoi(optarg); 
      break; 
    case 'm': 
      max_size_intron = atoi(optarg); 
      break; 
    case 'n': 
      nb_neighbour = atoi(optarg); 
      break;
    case 'o': 
      min_size_inside_exon = atoi(optarg); 
      break;
    case 'x': 
      extend_for_splicesites = atoi(optarg); 
      break;
    case 'v': 
      verbose = 1;
      break;
    case 'h': 
      usage(""); 
      break;
    default : 
      abort();
    }
  }

  min_nb_overlap_junction = (word_size - 2 * min_size_inside_exon);
  SSRContigList::MINOVERLAPJUNCTION = min_size_inside_exon;
  SSRContigList::NBNEIGHBOUR = nb_neighbour;
  SSRContigList::MINSIZEINTRON = min_size_intron;
  SSRContigList::MAXSIZEINTRON = max_size_intron;
  SSRContig::EXTEND = extend_for_splicesites;
  SSRContig::MINSIZEEXON = min_size_exon;
  SSRContig::MINCOVSPLICESITES = min_nb_overlap_junction;
  SSRContig::MINCOVWORD2ORIENTATE = min_cov_word2orientate;

  s32 nbSeq = 0;
  DnaDictionary dict(word_size);
  string sSeq;

  if (verbose) cerr << "Load reads and populate word dictionnary (wordsize= " << word_size << ")...\n";
  ReadFile reads(readsFilename);
  nbSeq = reads.loadAndCount(dict);

    
  if (verbose) cerr << "  word dictionnary completed : nbseq= " << 
		   nbSeq << " nbWords= " << dict.getNbWords() <<
		   " nbDiffWords= " << dict.getNbDiffWords() << 
		 " nbBadChars= " << dict.getNbBadChars() << "...\n";


  map<string, string> toto = loadFasta(fastaFilename);
  SSRContigLists liste(contigsFilename, toto);
  liste.testJunctions(dict);


  return 0;

}

map<string, string> loadFasta(char* filename) {
  char ch[1000];
  char c;
  string name, seqline;
  string sequence;
  map<string, string> fastaDB;
  
  fstream fstrm;
  fstrm.open(filename, ios_base::in | ios_base::binary);
  while(!fstrm.eof()) {
    fstrm>>c;
    if('>' == c) {
      if(!name.empty()) { 
		fastaDB.insert(make_pair(name,sequence));
		sequence.erase();
      }
      fstrm>>name; 
      fstrm.getline(ch,1000);
    } else {
		  fstrm>>seqline;
		  seqline.insert(0, 1, c);
		  transform(seqline.begin(), seqline.end(), seqline.begin(), (s32(*)(s32)) tolower);
		  sequence.append(seqline);
		  seqline.erase();
    }
  }
  if(!name.empty())
    fastaDB.insert(make_pair(name,sequence));
  fstrm.close();
  return fastaDB;
}

void usage(string msg) {
  cerr << msg << endl;
  cerr << "+---------------------------------------------------------------+" << endl;
  cerr << "| Creation de blocs de syntenie a partir d'HSP                  |" << endl;
  cerr << "+---------------------------------------------------------------+" << endl;
  cerr << "| genomodel [ -d file -q file ] {-e} {-E Ecore_file}            |" << endl;
  cerr << "|           {-D int} hsp_file                                   |" << endl;
  cerr << "|                                                               |" << endl;
  cerr << "|       -d file : les modeles sur la sequence Data sont cr%G�%@es   |" << endl;
  cerr << "|                 dans le fichier file.                         |" << endl;
  cerr << "|       -q file : les mod%G�%@les sur la sequence Query sont cr%G�%@es  |" << endl;
  cerr << "|                 dans le fichier file.                         |" << endl;
  cerr << "|                                                               |" << endl;
  cerr << "|       -e : si sortie lspextend ou back_translate              |" << endl;
  cerr << "|       -E file : Donne la composition des ecores dans le       |" << endl;
  cerr << "|                 fichier file                                  |" << endl;
  cerr << "|       -D int : Distance limite de fusion (distance entre      |" << endl;
  cerr << "|                ecores, 1 par d%G�%@faut)                          |" << endl;
  cerr << "|       hsp_file : fichier d'HSP (TBlastX Lassap)               |" << endl;
  cerr << "+---------------------------------------------------------------+" << endl;
  exit(1);
}

