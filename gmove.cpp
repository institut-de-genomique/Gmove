#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include <utility>
#include <time.h>
#include <ctime>
#include <malloc.h>
#include <seqan/basic.h>
#include <seqan/gff_io.h>
#include <seqan/stream.h>

#include "DnaDictionary.h"
#include "SSRContig.h"
#include "SSRContigLists.h"
#include "NetEx.h"
#include "GeneModel.h"
#include "ReadFile.h"
#include "GeneModelList.h"


using namespace std;
using namespace boost;
using namespace seqan;

void usage();
void error(string msg);
map<string, s32> loadProtPhase(char *filename);
map<string, string> loadFasta(char* filename);
list<pair <GffRecord, string > > loadGff(char* filename,bool b,map<string, map<s32,s32> >& cds,map<string,string>);
void  loadAnnotation(char* filename, list< pair< GffRecord,string > >& , map<string, map<s32,s32> >& cds,map<string,string>);
void fusionMonoExons(list<GffRecord> listMono,list< pair < GffRecord,string> >& tagRecord);
list< pair < GffRecord,string> > extractMono(list< pair < GffRecord,string> >& tagRecord);
bool sortRecordGff (const pair <GffRecord,string > & record1 , const pair <GffRecord,string> & record2);
bool sortMonoGff (const GffRecord & record1 , const GffRecord & record2);
bool existence(char* name);
void deleteIncludeMono(list< GffRecord >& listMono, list< pair < GffRecord,string> > tagRecord);
void cdsPhase(s32& phase, list<pair<s32,s32> > pos, s32 strand, string name, map<string,map<s32,s32> >& cds);
map<string,s32>  selectPhase(map<string,map<s32,s32> > cds);



int main(int argc, char** argv) {
  
  // default values of parameters
  s32 c, reduce_exon_list=100000, extend_for_splicesites=0, min_size_exon=9, min_size_intron=9, max_size_intron=50000, max_nb_paths=10000, search_window=30, nb_neighbour=20, verbose=1, min_size_cds=100;
  char *fastaFilename = NULL, *annotationFilename = NULL, *protPhaseFilename = NULL, *dataFilename = NULL, *dataProtFilename = NULL;
  string out = "",outRaw,outFilter;
  bool formatDef=true, keep_single_exon_gene=true, unoriented_junctions=false, rawData = false, ratioMono = false;
  s32 word_size = 25;

  // options given as parameters
  const struct option longopts[] =
    {
      {"rna",   required_argument,        0, 'r'},
	  {"prot",   required_argument,        0, 'd'},
	  {"annot",   required_argument,        0, 'a'},
	  {"out",   required_argument,        0, 'o'},
	  {"raw",no_argument ,0,'z'},
	  {"cds", required_argument, 0 ,'c'},
	  {"ratio",no_argument,0,'w'},
       {0,0,0,0},
    };
  int index;

  while ((c =  getopt_long(argc, argv, "f:c:a:r:d:o:SL:x:e:i:m:p:P:b:l:t:u:v:z:w:h",longopts,&index)) != -1) {//getopt(argc, argv, "f:c:j:a:D:d:G:C:J:SL:x:e:i:m:p:P:b:n:l:t:u:v:h")) != -1) {
    switch (c) {
    case 'f':
      fastaFilename = optarg; 
      break;
    case 'c': 
      min_size_cds = atoi(optarg);
      break;
    case 'a':
      annotationFilename = optarg;
      break;
    case 'r':
    	dataFilename = optarg;
    	break;
    case 'd':
    	dataProtFilename = optarg;
    	break;
    case 'o':
    	out = optarg;
      break;
    case 'P':
      protPhaseFilename = optarg;
      break;
    case 'S': 
      keep_single_exon_gene = false;
      break;
    case 'L':
      reduce_exon_list = atoi(optarg);
      break;
    case 'x': 
      extend_for_splicesites = atoi(optarg); 
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
    case 'p':
      max_nb_paths = atoi(optarg);
      break;
    case 'b':
      search_window = atoi(optarg);
      break;
    case 't':
      formatDef=false;
      break;
    case 'u':
      unoriented_junctions = true;
      break;
    case 'v': 
      verbose = 1;
      break;
    case 'z':
       	rawData = true;
       	break;
    case 'w':
    	ratioMono = true;
    	break;
    case 'h': 
      usage();
    }
  }

  // general constants
  SSRContigList::NBNEIGHBOUR = nb_neighbour;
  SSRContigList::MINSIZEINTRON = min_size_intron;
  SSRContigList::MAXSIZEINTRON = max_size_intron;
  SSRContigList::REDUCEEXONLIST = reduce_exon_list;
  SSRContigList::VERBOSE = verbose;
  SSRContig::EXTEND = extend_for_splicesites;
  SSRContig::MINSIZEEXON = min_size_exon;

  // OPTION FOR BLAT MODELS
  GeneModel::UNORIENTED = unoriented_junctions;

  // Check input and output options
  if(fastaFilename == NULL || !existence(fastaFilename)) error("Input genomic sequence : empty filename or file not found.");
  if(dataFilename != NULL && !existence(dataFilename)) error("Input rna file (gff) : empty file or file not found.");
  if(dataProtFilename != NULL && !existence(dataProtFilename)) error("Input prot file (gff) : empty file or file not found.");
  if(annotationFilename != NULL && !existence(annotationFilename)) error("Input annotation file (gff) : empty file or file not found.");
  if(protPhaseFilename != NULL && !existence(protPhaseFilename)) error("Input proteic phase file : file not found");

  // create network
  list<NetEx*> listReseau;
  if(verbose) cerr << "Load fasta file " << fastaFilename << "..." << endl;
  map<string, string> fastaDB = loadFasta(fastaFilename);

  // create output files
  	  //Load Rna seq Data
  time_t before = time(NULL);
  bool proteinBool = false;
  map<string,map<s32,s32> > cds;
  list<pair < GffRecord, string > > listRecord ;
  if(existence(dataFilename)){
	  if(verbose) cerr << "Load data file (gff) " << dataFilename << "..."<<endl;
	  listRecord = loadGff(dataFilename,proteinBool,cds,fastaDB);
	  //TODO sort data in Gmove
	  cerr << "load rnaSeq " << listRecord.size() << endl;
  }
  else
	  cerr<<"No rna file" << endl;

  if(existence(dataProtFilename)){
	  //Load protein data
	  if(verbose) cerr << "Load protein (gff) "<< dataProtFilename<<"..."<<endl;
	  proteinBool = true;
	  list<pair < GffRecord, string > > listRecord2 = loadGff(dataProtFilename,proteinBool,cds,fastaDB);
	  cerr << " load protein "<< listRecord2.size() << endl;
	  listRecord.insert(listRecord.end(),listRecord2.begin(),listRecord2.end());
	  listRecord2.clear();
  }
  else
	  cerr << "No protein file" << endl;

  time_t after = time(NULL);
 if(PRINTTIME) cout <<"time loadGff " << difftime(after,before) << endl;
  if(existence(annotationFilename)){
	  before = time(NULL);
	  if(verbose) cerr << "Load annotation file (gff) " << annotationFilename << "..."<< endl;
		loadAnnotation(annotationFilename,listRecord,cds,fastaDB);
		after = time(NULL);
		if(PRINTTIME) cout << "time load Annotation " << difftime(after,before) << endl;
  }
  else
	  cerr<<"No annotation file"<<endl;

  map<string, s32 > proteinPhaseCoord = selectPhase(cds);
  before = time(NULL);
  ofstream ofs_rawModelsFilename, ofs_filterModelsFilename;
  if(listRecord.empty()){
	  cerr<<"Warning : no data load "<<endl;
	  exit(0);
  }
  else{//Create ouput
	  //output raw data
	  if(out.empty()){
		cout << "No directory, need to create one "<<endl;
		system("mkdir -p out");
		out = "out/";
		system("mkdir -p out/filter");
	  }
	else {
		string cmd = "mkdir -p ";
		cmd += out;
		system(cmd.c_str());
		out +="/";

		string cmdFilter = "mkdir -p ";
		cmdFilter += out;
		cmdFilter += "filter";
		system(cmdFilter.c_str());
	}
	if(rawData){
		if(out.empty())
			system("mkdir -p out/raw");

		else{
			string cmdRaw;
			cmdRaw += "mkdir -p ";
			cmdRaw += out;
			cmdRaw +="raw";
			system(cmdRaw.c_str());
		}
		outRaw=out;
		outRaw += "raw/gmove_";
		string randName = fastaDB.begin()->first;
		outRaw += randName;
		char *tmpnameRaw = strdup(outRaw.c_str());
		mkstemp(tmpnameRaw);
		ofs_rawModelsFilename.open(tmpnameRaw, std::ofstream::out);
	}

	outFilter =out;
	cout << " out filter "<<outFilter<< endl;
	outFilter += "filter/gmove_";
	out +=  "gmove_";
	std::srand(std::time(0)); // use current time as seed for random generator
	out += "tmpfileXXXXXX";
	string randName = "tmpfileXXXXXX";
	outFilter += randName;
	char *tmpnameFilter = strdup(outFilter.c_str());
	mkstemp(tmpnameFilter);

	ofs_filterModelsFilename.open(tmpnameFilter, std::ofstream::out);
  }
//sortRecordGff
  //TODO loop on map scaff
  list< pair < GffRecord,string > > listMono = extractMono(listRecord);//TODO extract mono and put them in a map by scaff !
  map<string,list<GffRecord> > mapMono;
  for( list<pair < GffRecord, string > >::iterator itList = listMono.begin() ; itList != listMono.end(); ++itList){
	  map<string,list<GffRecord> >::iterator itMapMono;
	  string scaff;
	  assign(scaff,itList->first.ref);
	  itMapMono = mapMono.find(scaff);
	  if(itMapMono == mapMono.end()){

		  list<GffRecord> tmpList;
		  tmpList.push_back(itList->first);
		  pair< string,list<GffRecord> > tmpPair;
		  tmpPair = make_pair(scaff,tmpList);
		  mapMono.insert(tmpPair);
	  }
	  else
		  itMapMono->second.push_back(itList->first);
  }

  for(map<string,list<GffRecord> >::iterator itMap = mapMono.begin(); itMap != mapMono.end();++itMap){
	  fusionMonoExons(itMap->second,listRecord);
  }
  after = time(NULL);
  if(PRINTTIME) cout <<"time fusionMonoExons " << difftime(after,before) << endl;
  SSRContigLists liste(listRecord, fastaDB);
  DnaDictionary dict(word_size); // TODO doit pouvoir etre enleve --> changer fonctions correspondantes
  before = time(NULL);
  listReseau = liste.testJunctions(dict);
  after = time(NULL);
  if(PRINTTIME) cout << " time test junction " << difftime(after,before)<<endl;

  // create a directory for every .dot file
  system("mkdir -p dot_files");
  list<NetEx*>::iterator itNet;
   s32 cb=0;
  // test for cycles in the graph (not too sure if we keep that)
  for (itNet = listReseau.begin() ; itNet != listReseau.end() ; itNet++) {
	cout << " change network " << listReseau.size() << endl;
    NetEx* oneRes = *itNet;
    s32 num = std::distance(listReseau.begin(), itNet);
    ++num;
    cout << " Reseaux name " << oneRes->getSeqName()<< " iterator "<< num << endl;
    before = time(NULL);
    oneRes->cleanGraph();

    oneRes->graphOut("afterClean");
    after = time(NULL);
    if(PRINTTIME) cout << " time cleaning graph " << difftime(after,before) << endl;

    cout << "graph out... "<<endl;
    cout << "finish graph out "<<endl;
    map<string, list<s32> > probExons;
    s32 count_cycle = oneRes->nbCycle();
    if(count_cycle) {
      cerr << "[Error] The graph is cyclic in sequence " << oneRes->getSeqName() << endl;
      cerr << "[Error] Number of cycles in the graph is " << count_cycle << endl;
      cerr << "[Error] Continue with next input sequence..." << endl;
      continue;
    }
    vector<list<s32> > components = oneRes->getComponents();
    list<list<s32> > allpaths;
    s32 nb_models=0, nb_used_cc=0, nb_models_filter =0;
    u32 cutoff_nbExons=1;
    if(keep_single_exon_gene) cutoff_nbExons=0;
    GeneModelList allModelFromOneRes;
    for (u32 c = 0; c < components.size(); c++) {
    	time_t before_cc = time(NULL);
    	cout << " change cc "<< c<< endl;
      list<s32> component =  components[c];
      if (component.size() > cutoff_nbExons){ //1 -> nb membre de la composante /! 1 after cleaned
    	  cb++;
    	  s32 nbedges=0, nbvertices=0;
    	  oneRes->countVerticesAndEgdes(component, nbvertices, nbedges);

    	  if(nbedges > 10 * nbvertices) {
    		  if(verbose) cerr << "[Warning] found one connected component with number of edges greater than 10 times the number of vertices (nbvertices= "
			   << nbvertices << " nbedges= " << nbedges << ")" << endl;
    		  continue;
    	  }
		  s32 nbPaths = oneRes->count_allPaths(component);
		  cerr << " nbPaths " << nbPaths <<" size component "<< component.size()<< endl;
		  // calculate coordinates of too furnished region
		  s32 threshold =2;
		  while(nbPaths > max_nb_paths || nbPaths == -1) {

			  if(verbose) {
				  s32 begReg = -1, endReg = 0;
				  for(list<s32>::const_iterator itComp = component.begin(); itComp != component.end(); itComp++){
					  TSSRList::iterator it1 = oneRes->getVertices()->begin();
					  advance(it1,*itComp);
					  SSRContig* ctgTag = *it1;
					  if (begReg == -1 || begReg > ctgTag->start() ) begReg = ctgTag->start();
					  	  if (ctgTag->end() > endReg) endReg = ctgTag->end();
				  }

				  cerr <<"threshold "<< threshold <<" [Warning] Too many paths ( " << nbPaths << " ) : " << oneRes->getSeqName() << "\t" << begReg << "\t" << endReg << endl;

				  oneRes->simplifyBigGraph(component,threshold); //simplify big graph
				  components = oneRes->getComponents();
				  component = components[c];//replace former component with the new simplify one
				  nbPaths = oneRes->count_allPaths(component);

			  ++threshold;
			  }
		  }


		  TSSRList* tmpVertex;
		  tmpVertex = oneRes->getVertices();
		  cout << "out of while new nbPath " << nbPaths<< endl;
		  before = time(NULL);

		  allpaths= oneRes->PathsFinderWithCondition(component);
	//	  allpaths = oneRes->allPathsFinder(component);
		  after = time(NULL);
		  if(PRINTTIME) cout << " time allPathFinder " << difftime(after,before)<< " number of paths " << allpaths.size() << endl;
		  nb_used_cc++;
		  s32 nbPath=1;
		  // find a model for each path in the connected component
		  before = time(NULL);
		  MapModel mapGeneModel;
		  cout << "\n\n";
		  for(list<list<s32> >::iterator itPa = allpaths.begin(); itPa != allpaths.end(); itPa++) {
	//		  for(list<s32>::iterator itPaa = itPa->begin(); itPaa != itPa->end(); itPaa++) {
	//			  cout << *itPaa << " " ;
	//		  }
	//		  cout << endl;
			  GeneModel gene = GeneModel(*itPa, oneRes->getVertices(), oneRes->getSeqName(), c+1, search_window, probExons);//TODO save all gene from one cc
			  std::clock_t c_start = std::clock();
			  bool contain_orf = gene.findORF(proteinPhaseCoord);
			  std::clock_t c_end = std::clock();
			 if(PRINTTIME)
				  cout << " time findORF " <<  1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << endl;
			  // On imprime une seule CDS par transcrit
			  if(contain_orf != gene.containCDS())
				  cerr << "Error contain_orf " << contain_orf << " gene.containCDS() " << gene.containCDS()<< (*gene.getExons().begin())->start()<< " end exon " << (*gene.getExons().end())->end()<<endl;
			  if(gene.containCDS()) {
				  gene.keepModel(mapGeneModel);
			  }
		  }
		  before = time(NULL);
		  vector<GeneModel> vectorOfGene;
		  nbPath=1;
		  for(MapModel::iterator itMap = mapGeneModel.begin(); itMap != mapGeneModel.end(); ++itMap){
			  itMap->second.setID(nbPath);
			  itMap->second.printAnnot(ofs_rawModelsFilename, formatDef);
			  nbPath++;
			  nb_models++;

		  }
		  std::clock_t beforeClean = std::clock();
		  GeneModelList listOfGene(mapGeneModel);//FIXME another way to convert map to list ? splice ??
		  std::clock_t afterClean = std::clock();

		  beforeClean = time(NULL);
		  listOfGene.deleteIncludedModel();
		  afterClean = time(NULL);
	if(PRINTTIME)	  cout << " time deleteIncludedModel " << 1000.0 * (afterClean-beforeClean) / CLOCKS_PER_SEC  << endl;


	  beforeClean = time(NULL);
	  	  if(listOfGene.getSize() > 1)
			  listOfGene.longestCDS();
		  afterClean = time(NULL);
	if(PRINTTIME)	  cout << " time deleteSmallCDS and print annot " << 1000.0 * (afterClean-beforeClean) / CLOCKS_PER_SEC << endl;
		  allModelFromOneRes.insertModels(listOfGene); //FIXME What is it for ? copy de list, pb ! mapGeneModel -> listOfGene -> allModelFromOneRes Pourquoi listOfGene existe ?
     }
      after = time(NULL);
      if(PRINTTIME) cout << " time to loop on cc " << difftime(after,before_cc) << endl;
    }
    cout << "allModelfromOneRes " << allModelFromOneRes.getSize() << endl;


    before= time(NULL);
    allModelFromOneRes.includedMono();
    after = time(NULL);
    cout << " time includedMono " << 1000.0 * (after-before) / CLOCKS_PER_SEC << endl;

    before= time(NULL);
    allModelFromOneRes.deleteSmallCDS(min_size_cds);
    after = time(NULL);
  if(PRINTTIME)  cout << " time deleteSmallCds " << 1000.0 * (after-before) / CLOCKS_PER_SEC << endl;

    before= time(NULL);
    allModelFromOneRes.clusterLocation();
    after = time(NULL);
  if(PRINTTIME)  cout << " time ClusterLocation " << 1000.0 * (after-before) / CLOCKS_PER_SEC << endl;


    if(ratioMono)
    	allModelFromOneRes.ratioCdsUtr();

    GeneModelL tmpModels = allModelFromOneRes.getModels() ; //FIXME do not create a copy !
    for(GeneModelL::iterator itList = tmpModels.begin(); itList != tmpModels.end() ; ++itList){//FIXME print directly in GeneModelL !
		if(!itList->getToDelete()){
			itList->printAnnot(ofs_filterModelsFilename, formatDef);
			++nb_models_filter;
		}
  	  }

    if(verbose) {
      list<pair<s32,s32> > mypbCC;
      cerr <<"Number of usable connected component(s): " << nb_used_cc << endl;
      cerr <<"Number of gene model(s): " << nb_models << endl;
      cerr << "Number of gene model(s) after filtering : " << nb_models_filter << endl;
      for(map<string, list<s32> >::iterator itPE = probExons.begin();itPE != probExons.end(); itPE++){
    	  itPE->second.sort();
    	  itPE->second.unique();
    	  string tmp_cc_nb;
    	  if(itPE->second.front() != itPE->second.back())
    		  mypbCC.push_back( make_pair(itPE->second.front(),itPE->second.back()) );
      }
      mypbCC.sort();
      mypbCC.unique();
      list<pair<s32,s32> >::iterator itCC;
      for(itCC = mypbCC.begin(); itCC != mypbCC.end(); itCC++){
    	  cerr << "[Warning] Split connected components : mRNA." << oneRes->getSeqName() << "." << (*itCC).first << "\tmRNA." << oneRes->getSeqName() << "." << (*itCC).second << endl;
      }
    }
    delete oneRes;
  }
  return 0;
}

// to load into a hash table the proteic mapping phase information
map<string, s32> loadProtPhase(char *filename){
  fstream fstrm;
  string seqname;
  s32 coord, phase;//phase = 1,2 or 3
  map<string, s32> protPhase; // hash table [key]=phase

  fstrm.open(filename, ios_base::in | ios_base::binary);
  while(!fstrm.eof()) {
    seqname.erase();
    fstrm>>seqname;
    if(seqname.empty()) { break; }
    fstrm>>coord;
    fstrm>>phase;
    ostringstream oss;
    oss << seqname << "@" << coord;
    string key = oss.str();
    map<string, s32>::iterator itSeq = protPhase.find(key);
    if(itSeq == protPhase.end())
      protPhase.insert( make_pair(key,phase) );
  }
  fstrm.close();
  return protPhase;
}

// to load the fasta file of the genome into a hash table
map<string, string> loadFasta(char* filename) {
  char ch[1000];
  string name, seqline;
  string sequence;
  map<string, string> fastaDB;// hash table [name]=sequence ( sequence in lower case atcg)
  fstream fstrm;
  fstrm.open(filename, ios_base::in | ios_base::binary);
  while(!fstrm.eof()) {
    fstrm>>seqline;
    if('>' == seqline[0]) {
      if(!name.empty()) { 
    	  cerr << "look at scaff " << name << endl;
    	  fastaDB.insert(make_pair(name,sequence));
    	  sequence.erase();
      }
      seqline.erase(0,1);
      name= seqline;
      fstrm.getline(ch,1000);
    } else {
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

list< pair< GffRecord,string > > loadGff(char* filename,bool protein,map<string, map<s32,s32> >& cds,map<string,string>fastaDB){
	/*
	We need to keep in mind where come from the exon, we cannot build all combinaison from this data
	*/
	GffFileIn gffIn(filename);
	 // Copy the file record by record.
	 GffRecord record;
	 list< pair <GffRecord,string > > tagRecord;
	 s32 phase,lengthCds = 0;
	 map<string,list<pair<s32,s32> > > posCds;
	 map<string,s32> selectedCds;
	 map<string,bool> idTagRecord;
	 string previous_exon_id,tag,previous_tag;
	 while (!atEnd(gffIn)){
		 readRecord(record, gffIn);

		 string scaff;
		 assign(scaff,record.ref);

		 if(!tagRecord.empty() && tagRecord.back().first.tagValues[0] == ""){
			 cerr << " error in LoadGff : name empty (last column) "<< endl;
		 }
		 map<string,string>::iterator itFasta = fastaDB.find(scaff);
		 if(itFasta==fastaDB.end())continue;

		 if(record.type=="exon" || record.type=="HSP"){
			 if(!tagRecord.empty() && record.tagValues[0] == tagRecord.back().first.tagValues[0]){
				 tag = 'i';
				 if(protein){
					lengthCds += record.endPos - record.beginPos;
					pair < int, int > pairPos = make_pair(record.beginPos,record.endPos);
					string tmp ;
					assign(tmp,record.tagValues[0]);
					posCds[tmp].push_back(pairPos);
				}
			 }
			 else{
				 if(!tagRecord.empty()){
					 string tmp;
					 assign(tmp,tagRecord.back().first.tagValues[0]);
					 map<string,bool>::iterator itIdTagRecord = idTagRecord.find(tmp);
					 if(itIdTagRecord != idTagRecord.end()){
						 cerr << "Error : id is not uniq : " << itIdTagRecord->first<<endl;
						 exit(1);
					 }
					 else{
						 pair<string,bool> tmpPair = make_pair(tmp,true);
						 idTagRecord.insert(tmpPair);
					 }
				 }
				 tag = 'd';
				 if(!tagRecord.empty() && tagRecord.back().second == "d"){
					 tagRecord.back().second = 'm';
					 if(protein ){
						 if(tagRecord.back().first.strand=='+')phase = 1;
						 else phase = -1;
						 string name;
						 assign(name,tagRecord.back().first.ref);
						 string tmp ;
						 assign(tmp,tagRecord.back().first.tagValues[0]);
						 cdsPhase(phase,posCds[tmp],tagRecord.back().first.strand,name,cds);
					 }
				 }

				 else{
					if(!tagRecord.empty()){
						tagRecord.back().second = 'f';
						 if(protein ){
							 if( lengthCds%3 == 0 ){
								 if(tagRecord.back().first.strand=='+')phase = 1;
								 else phase = -1;
								 string name;
								 assign(name,tagRecord.back().first.ref);
								 list<pair < s32,s32 > > pos;
								 string tmp ;
								 assign(tmp,tagRecord.back().first.tagValues[0]);
								 cdsPhase(phase,posCds[tmp],tagRecord.back().first.strand,name,cds);
							 }
						 }
					}
				}
				 if(protein){
					lengthCds = record.endPos - record.beginPos;
					pair < int, int > pairPos = make_pair(record.beginPos,record.endPos);
					string tmp ;
					assign(tmp,record.tagValues[0]);
					posCds[tmp].push_back(pairPos);
				}
		 }
		pair<GffRecord, string> tmpPair = make_pair(record,tag);
		tagRecord.push_back(tmpPair);
		clear(record);
		 }
	 }
	 //update the last one
	 if(tagRecord.empty()){
		 cerr << " Warning : we did not charge any data from "<< filename << endl;
	 }
	 else{
		 if(tagRecord.back().second =="i")
			 tagRecord.back().second ='f';
		 else tagRecord.back().second ='m';
	 }
	 close(gffIn);
	 return tagRecord;
}

//-------------------------------------------------------------------------------------

void loadAnnotation(char* filename, list< pair< GffRecord,string > >& tagRecord, map<string, map<s32,s32> >& cds,map<string,string>fastaDB){
	/* use it for reannotation
	We need to keep in mind where come from the exon, we cannot build all combinaison from this data
	*/
	map< string,  list < pair < s32,s32> > > listExon,posCds;

	s32 phase, lengthCds = 0 ;//phase = 1,2 or 3
	map<string,s32> selectedCds;
	GffFileIn gffIn(filename);
	 // Copy the file record by record.
	 GffRecord record;
	 string previous_exon_id,tag,previous_tag;
	 while (!atEnd(gffIn)){
		readRecord(record, gffIn);
		string scaff;
		assign(scaff,record.ref);
		map<string,string>::iterator itFasta = fastaDB.find(scaff);
		if(itFasta==fastaDB.end())continue;
		if(record.type!="CDS" && record.type!="exon" && record.type!="UTR"){ continue;}
		if(!tagRecord.empty() && record.tagValues[0] == tagRecord.back().first.tagValues[0]){
			tag = 'i';
			string name;
			assign(name,record.tagValues[0]);
			if(record.type == "CDS"){
				lengthCds += record.endPos - record.beginPos;
				pair < int, int > pairPos = make_pair(record.beginPos,record.endPos);
				posCds[name].push_back(pairPos);
			}
			if(tagRecord.back().first.endPos == record.beginPos){
				if(listExon.find(name) == listExon.end()){
					pair < int, int > pos (tagRecord.back().first.beginPos, record.endPos);
					listExon[name].push_back(pos);
				}
				else{
					map< string, list<pair < s32,s32> > >::iterator itListPos = listExon.find(name);
					for( list<pair < s32,s32 > >::iterator itPos = itListPos->second.begin() ; itPos != itListPos->second.end(); ++itPos){
						if(itPos->second == record.beginPos){
							itPos->second = record.endPos;
							break;
						}
					}
					tagRecord.back().first.endPos = record.endPos;
					continue;
				}
			}
		}
		else{
			tag = 'd';
			string name;
			assign(name,record.tagValues[0]);
			pair < int, int > pos (record.beginPos,record.endPos);
			listExon[name].push_back(pos);
			if(!tagRecord.empty() && tagRecord.back().second == "d")
				tagRecord.back().second = 'm';
			else{
				if(!tagRecord.empty()){
					tagRecord.back().second = 'f';
						if( lengthCds%3 == 0 ){
							if(tagRecord.back().first.strand=='+')phase = 1;
							else phase = -1;
							string name;
							assign(name,tagRecord.back().first.ref);
							list<pair < s32,s32 > > pos;
							string tmp ;
							assign(tmp,tagRecord.back().first.tagValues[0]);
							cdsPhase(phase, posCds[tmp],  tagRecord.back().first.strand, name, cds);
						}
				}
			}
			lengthCds = 0;
			if(record.type == "CDS"){
				lengthCds = record.endPos - record.beginPos;
				pair < int, int > pairPos = make_pair(record.beginPos,record.endPos);
				string tmp ;
				assign(tmp,record.tagValues[0]);
				posCds[tmp].push_back(pairPos);
			}
		}
		pair<GffRecord, string> tmpPair = make_pair(record,tag);
		tagRecord.push_back(tmpPair);
		clear(record);
	}
	 //update the last one
	 if(tagRecord.empty()){
			 cerr << " Warning : we did not charge any data from "<< filename << endl;
		 }
		 else{
	 if(tagRecord.back().second =="i") tagRecord.back().second ='f';
	 else tagRecord.back().second ='m';
		 }
	 close(gffIn);

}

void cdsPhase(s32& phase, list< pair< s32,s32> > pos, s32 strand, string name, map<string,map<s32,s32> >& cds){
	if(strand=='+'){
		for(list<pair <s32,s32 > >::iterator itPos = pos.begin(); itPos != pos.end(); ++itPos){
			s32 start = (*itPos).first;
			s32 end = (*itPos).second;
			for(int i =start+1 ; i <= end ; ++i){
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

map<string,s32> selectPhase(map<string,map<s32,s32> > cds){
	map<string,s32> phase;
	for(map<string,map<s32,s32 > >::iterator itCds  = cds.begin() ; itCds != cds.end() ; ++itCds){
		s32 max = 0,secondMax = 0 ,maxPhase = 0 ;
		for(map<s32,s32>::iterator itPhase = itCds->second.begin() ; itPhase != itCds->second.end(); ++itPhase){
			if(itPhase->second >= max){
				secondMax = max;
				max = itPhase->second;
				maxPhase = itPhase->first;
			}
		}
		if(max == secondMax )
			phase[itCds->first] = 0;

		else
			phase.insert(make_pair(itCds->first,maxPhase));
	}
	return phase;
}


list< pair < GffRecord,string> > extractMono(list< pair < GffRecord,string> >& tagRecord){
	list< pair < GffRecord,string> > listMono;
	for (list< pair < GffRecord,string> >::iterator it=tagRecord.begin() ;it != tagRecord.end();){
		if(it->second == "m"){
			listMono.push_back(*it);
			it = tagRecord.erase(it);
		}
		else
			++it;
	}
	return listMono;
}

void fusionMonoExons(list<GffRecord> listMono,list< pair < GffRecord,string> >& tagRecord){//TODO move that to GeneModel !
	//Same idea as /env/ig/soft/rdbioseq/annotation-snapshot/linux-noarch/bin/loci
	listMono.sort(sortMonoGff);
	for(list< GffRecord>::iterator itListMono = listMono.begin(); itListMono != --listMono.end() ;){//--listMono.end() : we stop one element before the end, like size()-1
		list< GffRecord >::iterator itNext = itListMono;
		++itNext;
		if(itListMono->ref == itNext->ref && itListMono->strand == itNext->strand){ //(itListMono->strand == itNext->strand || itListMono->strand=='.' || itNext->strand == '.')){//&& itListMono->strand == itNext->strand){ FIXME fusionne mono with . or + or - ??
			if( (itNext->endPos >= itListMono->beginPos && itNext->endPos <= itListMono->endPos)
			|| (itListMono->endPos >= itNext->beginPos && itListMono->endPos <= itNext->endPos)){
				itNext->beginPos = min(itListMono->beginPos,itNext->beginPos);
				itNext->endPos = max(itListMono->endPos,itNext->endPos);
				itListMono = listMono.erase(itListMono);
			}
			else
				++itListMono;
		}
		else
			++itListMono;
	}
	for(list<GffRecord>::iterator itMono  = listMono.begin(); itMono != listMono.end(); ++itMono){
		tagRecord.push_back(make_pair(*itMono,"m"));
	}
}


void deleteIncludeMono(list<GffRecord>& listMono, list< pair < GffRecord,string> > tagRecord){
	//TODO before add at tagRecord, look if mono is inside a pluriexon : delete mono
	for(list< GffRecord>::iterator itMono = listMono.begin(); itMono != listMono.end();){
		bool boolDelete =false;
		for(list< pair < GffRecord,string > >::iterator itOther = tagRecord.begin(); itOther != tagRecord.end(); ++itOther){
			if(itMono->beginPos >= itOther->first.beginPos && itMono->endPos <= itOther->first.endPos && itMono->ref == itOther->first.ref && itMono->strand == itOther->first.strand){
				boolDelete = true;
				break;
			}
		}
			if(boolDelete == true)
				itMono = listMono.erase(itMono);
			else
				++itMono;
	}
}
bool sortRecordGff (const pair <GffRecord,string >& record1 , const pair <GffRecord,string>& record2){
	return (record1.first.ref <= record2.first.ref &&  record1.first.beginPos <= record2.first.beginPos);
}

bool sortMonoGff (const GffRecord& record1 , const GffRecord & record2){
	return (record1.beginPos <= record2.beginPos);
}
// to send error message
void error(string msg) {
  cerr << "[Error] " << msg << endl;
  cerr << "See gmove -h for more details." << endl;
  exit(1);
}

// to check for existence
bool existence(char* name) {
  ifstream f(name);
  if (f.good()) {
    f.close();
    return true;
  }
  else {
    f.close();
    return false;
  }   
}

// the result of gmove -h
void usage() {
	cerr << "--------------------------------------------------------------------------------------------" << endl;

	cerr << "gmove - Gene modelling using various evidence." << endl << endl;
	cerr << "Usage : gmove -f <reference sequence> -c <covtigs> -j <junctions> -G <output gene models> {Options}" << endl;
	cerr << "!! Note : Arguments with * are required, the other are optionnal !!" << endl << endl;
	cerr << "  INPUT FILES" << endl;
	cerr << "*    -f <file> 		: fasta file which contains genome sequence(s)." << endl << endl;
	cerr << "     --annot <file> 	: annotation file in gff "<<endl;
	cerr << "     --rna <file> 		: rna file in gff "<<endl;
	cerr << "     --prot <file> 	: prot file in gff"<<endl;


	cerr << "  OUTPUT FILES" << endl;
	cerr << "     -o <folder> 	: output folder, by default (./out) " << endl;
	cerr << "     --raw 		: output raw data " << endl;
	cerr << endl;
	cerr << "     -S       		: do not output single exon models." << endl;
	cerr << "     -e <int>  	: minimal size of exons, default is 9 nucleotides." << endl;
	cerr << "     -i <int>  	: minimal size of introns, default is 9 nucleotides." << endl;
	cerr << "     -m <int>  	: maximal size of introns, default is 50.000 nucleotides." << endl;
	cerr << "     -p <int>  	: maximal number of paths inside a connected component, default is 10,000." << endl;
	cerr << "     -b <int>  	: number of nucleotides around exons boundaries where to find start and stop codons, default is 30." << endl;
	cerr << "     -t        	: gtf format annotation file - default is gff3" << endl;
	cerr << "     -u        	: choose model strand according to longest ORF - only works if input junctions are non-oriented" << endl;
	cerr << "     --cds         : min size CDS, by default 100 " << endl;
	cerr << "     -v        	: silent mode" << endl;
	cerr << "     -h        	: this help" << endl;
	cerr << " --ratio 			: ratio CDS/UTR min 80% de CDS" << endl;
	cerr << "--------------------------------------------------------------------------------------------" << endl;
	exit(1);
}
