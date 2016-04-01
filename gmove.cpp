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


using namespace std;
using namespace boost;
using namespace seqan;

void usage();
void error(string msg);
map<string, s32> loadProtPhase(char *filename);
map<string, string> loadFasta(char* filename);
list<pair <GffRecord, string > > loadGff(char* filename,bool b,map<string, map<s32,s32> >& cds,map<string,string>);
void  loadAnnotation(char* filename, list< pair< GffRecord,string > >& , map<string, map<s32,s32> >& cds,map<string,string>);
void fusionMonoExons(list< pair < GffRecord,string> >& tagRecord);
list< pair < GffRecord,string> > extractMono(list< pair < GffRecord,string> >& tagRecord);
bool sortRecordGff (const pair <GffRecord,string > & record1 , const pair <GffRecord,string> & record2);
bool sortMonoGff (const pair <GffRecord,string > & record1 , const pair <GffRecord,string> & record2);
bool existence(char* name);
void deleteIncludeMono(list< pair < GffRecord,string > >& listMono, list< pair < GffRecord,string> > tagRecord);
void cdsPhase(s32& phase, list<pair<s32,s32> > pos, s32 strand, string name, map<string,map<s32,s32> >& cds);
map<string,s32>  selectPhase(map<string,map<s32,s32> > cds);

int main(int argc, char** argv) {
  
  // default values of parameters
  s32 c, reduce_exon_list=100000, extend_for_splicesites=0, min_size_exon=25, min_size_intron=9, max_size_intron=50000, max_nb_paths=10000, search_window=30, nb_neighbour=20, verbose=1;
  char *fastaFilename = NULL, *contigsFilename = NULL, *junctionsFilename = NULL, *annotationFilename = NULL, *out_genemodelsFilename = NULL, *out_contigsFilename = NULL, *out_junctionsFilename = NULL, *protPhaseFilename = NULL, *dataFilename = NULL, *dataProtFilename = NULL;
  string out = "";
  bool formatDef=true, keep_single_exon_gene=true, unoriented_junctions=false;
  s32 word_size = 25;

  // options given as parameters
  const struct option longopts[] =
    {
      {"rna",   required_argument,        0, 'r'},
	  {"prot",   required_argument,        0, 'd'},
	  {"annot",   required_argument,        0, 'a'},
	  {"out",   required_argument,        0, 'o'},
       {0,0,0,0},
    };
  int index;

  while ((c =  getopt_long(argc, argv, "f:c:j:a:r:d:o:C:J:SL:x:e:i:m:p:P:b:n:l:t:u:v:h",longopts,&index)) != -1) {//getopt(argc, argv, "f:c:j:a:D:d:G:C:J:SL:x:e:i:m:p:P:b:n:l:t:u:v:h")) != -1) {
    switch (c) {
    case 'f':
      fastaFilename = optarg; 
      break;
    case 'c': 
      contigsFilename = optarg; 
      break;
    case 'j': 
      junctionsFilename = optarg; 
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
     // out_genemodelsFilename = optarg;
      break;
    case 'C' :
      out_contigsFilename = optarg;
      break;
    case 'J': 
      out_junctionsFilename = optarg; 
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
    case 'n': 
      nb_neighbour = atoi(optarg); 
      break;
    case 't':
      formatDef=false;
      break;
    case 'u':
      unoriented_junctions = true;
      break;
    case 'v': 
      verbose = 0; // TODO verbose =1
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
  if(contigsFilename != NULL && !existence(contigsFilename)) error("Input contigs file : empty filename or file not found.");
  if(junctionsFilename != NULL && !existence(junctionsFilename)) error("Input junctions file : empty filename or file not found.");
  if(fastaFilename == NULL || !existence(fastaFilename)) error("Input genomic sequence : empty filename or file not found.");
  if(dataFilename != NULL && !existence(dataFilename)) error("Input rna file (gff) : empty file or file not found.");
  if(dataProtFilename != NULL && !existence(dataProtFilename)) error("Input prot file (gff) : empty file or file not found.");
  if(annotationFilename != NULL && !existence(annotationFilename)) error("Input annotation file (gff) : empty file or file not found.");
 // if(out_genemodelsFilename == NULL) error("Need to provide an output file for gene models.");
  
  if(out_junctionsFilename != NULL && junctionsFilename == out_junctionsFilename) error("Need to provide two distinct files for input and output of junctions.");
  if(out_contigsFilename != NULL && contigsFilename == out_contigsFilename) error("Need to provide two distinct files for input and output of covtigs.");
  if(protPhaseFilename != NULL && !existence(protPhaseFilename)) error("Input proteic phase file : file not found");

  //prepare output :

  //  if(out_junctionsFilename!=NULL){
  ofstream ofs_junctionsFilename;
  ofs_junctionsFilename.open(out_junctionsFilename);
  //}
  //if(out_contigsFilename !=NULL){
  ofstream ofs_contigsFilename;
  ofs_contigsFilename.open(out_contigsFilename);
  //}

  // create network
  list<NetEx*> listReseau;
  if(verbose) cerr << "Load fasta file " << fastaFilename << "..." << endl;
  map<string, string> fastaDB = loadFasta(fastaFilename);


  // create output files
  	if(out.empty()){
  		cout << "No directory, need to create one "<<endl;
  		 system("mkdir -p out");
  		 out = "out/";

  	}
  	else if (*(out.rbegin())!= '/'){
  		char cmd[1000];
  		strcat(cmd,"mkdir -p ");
  		strcat(cmd,out.c_str());
  		 system(cmd);
  		out +="/";
  	}
  	else{
  		char cmd[1000];
  		strcat(cmd, "mkdir -p ");
  		  		strcat(cmd,out.c_str());
  		  		 system(cmd);
  	}

    ofstream ofs_genemodelsFilename;
    out +=  "gmove_";
    string s;
    std::srand(std::time(0)); // use current time as seed for random generator
    out += "tmpfileXXXXXX";
    char *tmpname = strdup(out.c_str());
    mkstemp(tmpname);
    cout << "File name of the output "<< tmpname << endl;
    ofs_genemodelsFilename.open(tmpname, std::ofstream::out);
    //  ofs_genemodelsFilename.open((char*)(fastaDB.begin()->first).c_str(), std::ofstream::out);
 // if(formatDef==true) ofs_genemodelsFilename << "##gff-version 3" << endl;
  //prepare gff files
 // char cmd[1000];
 // strcat(cmd,"sort -k1,1 -k10,10 -k4,4n ");
 // strcat(cmd,dataFilename);
  char * tmpDataFilename = dataFilename;//"tmpSortData.gff";
 // strcat(cmd," > tmpSortData.gff");
 // strcat(cmd,tmpDataFilename);
 // system(cmd);
 // if(system(cmd) != 0 )
//	  cerr << " error cannot sort and create a temporary file "<< endl;
  //Load Rna seq Data
  time_t before = time(NULL);
  bool proteinBool = false;
  map<string,map<s32,s32> > cds;
  list<pair < GffRecord, string > > listRecord ;
  if(existence(dataFilename)){
	  if(verbose) cerr << "Load data file (gff) " << dataFilename << "..."<<endl;
	  listRecord = loadGff(dataFilename,proteinBool,cds,fastaDB);
	  cerr << "load rnaSeq " << listRecord.size() << endl;
  }
  else{
	  cerr<<"No rna file" << endl;
  }
  if(existence(dataProtFilename)){
	  //Load protein data
	  if(verbose) cerr << "Load protein (gff) "<< dataProtFilename<<"..."<<endl;
	  char * tmpDataFilename = dataProtFilename;//"tmpSortData.gff";
	  proteinBool = true;
	  list<pair < GffRecord, string > > listRecord2 = loadGff(dataProtFilename,proteinBool,cds,fastaDB);
	  cerr << " load protein "<< listRecord2.size() << endl;
	  listRecord.insert(listRecord.end(),listRecord2.begin(),listRecord2.end());
	  listRecord2.clear();
  }
  else{
	  cerr << "No protein file " << endl;
  }

  time_t after = time(NULL);
 if(PRINTTIME) cout <<"time loadGff " << difftime(after,before) << endl;
  if(existence(annotationFilename)){
	  before = time(NULL);
	  if(verbose) cerr << "Load annotation file (gff) " << annotationFilename << "..."<< endl;
		loadAnnotation(annotationFilename,listRecord,cds,fastaDB);
		after = time(NULL);
		if(PRINTTIME) cout << "time load Annotation " << difftime(after,before) << endl;
  }
  else{
	  cerr<<"No annotation file"<<endl;
  }
  map<string, s32 > proteinPhaseCoord = selectPhase(cds);
  before = time(NULL);
  if(listRecord.empty()){
	  cerr<<"Error : no data load"<<endl;
  }
  fusionMonoExons(listRecord);
  after = time(NULL);
  if(PRINTTIME) cout <<"time fusionMonoExons " << difftime(after,before) << endl;
  SSRContigLists liste(listRecord, fastaDB);

  if(protPhaseFilename != NULL) {
    if(verbose) cerr << "Load input protein mapping phase information from " << protPhaseFilename << " ..." << endl;
    map<string, s32 > proteinPhase = loadProtPhase(protPhaseFilename);
    proteinPhaseCoord.insert(proteinPhase.begin(),proteinPhase.end());
  }

  DnaDictionary dict(word_size); // TODO doit pouvoir etre enleve --> changer fonctions correspondantes
   // to test junctions between the covtigs
  before = time(NULL);
  listReseau = liste.testJunctions(dict, ofs_junctionsFilename);
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
    before = time(NULL);
    oneRes->cleanGraph();
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
    s32 nb_models=0, nb_used_cc=0;
    u32 cutoff_nbExons=1;
    if(keep_single_exon_gene) cutoff_nbExons=0;
    for (u32 c = 0; c < components.size(); c++) { //TODO ne pas faire de boucle mais découper les clusters
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
		  cout << " nbPaths " << nbPaths << endl;
		  // calculate coordinates of too furnished region
		  if(nbPaths > max_nb_paths || nbPaths == -1) {
			  if(verbose) {
				  s32 begReg = 0, endReg = 0;
				  for(TSSRList::iterator itTag = oneRes->getVertices()->begin(); itTag != oneRes->getVertices()->end(); itTag++){ //XXX we always need double loop? To look at all vertices and just get the one that correspond to my component ?
					  SSRContig* ctgTag = *itTag;
					  for(list<s32>::const_iterator itComp = component.begin(); itComp != component.end(); itComp++){
						  s32 iD = *itComp;
						  if(ctgTag->getID() == iD) {
							  if (begReg == 0) begReg = ctgTag->start();
							  if (ctgTag->end() > endReg) endReg = ctgTag->end();
							  break;
						  }
					  }
				  }
				  cerr << "[Warning] Too many paths ( " << nbPaths << " ) : " << oneRes->getSeqName() << "\t" << begReg << "\t" << endReg << endl;
				  continue;
			  }
		  }
		  if(nbPaths < 1) continue;
		  before = time(NULL);

		  allpaths = oneRes->allPathsFinder(component);
		  after = time(NULL);
		  if(PRINTTIME) cout << " time allPathFinder " << difftime(after,before)<< " number of paths " << allpaths.size() << endl;
		  nb_used_cc++;
		  s32 nbPath=1;
		  // find a model for each path in the connected component
		  before = time(NULL);
		  map<string,GeneModel> mapGeneModel;
		  for(list<list<s32> >::iterator itPa = allpaths.begin(); itPa != allpaths.end(); itPa++) {

			  GeneModel gene = GeneModel(*itPa, oneRes->getVertices(), oneRes->getSeqName(), cb, search_window, probExons);//TODO save all gene from one cc
			  time_t before2 = time(NULL);
			  std::clock_t c_start = std::clock();
			  bool contain_orf = gene.findORF(proteinPhaseCoord);
			  std::clock_t c_end = std::clock();
			  time_t after2 = time(NULL);
			  if(PRINTTIME)  cout << " time findORF " <<  1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << endl;
			  // On imprime une seule CDS par transcrit
			  if(contain_orf != gene.containCDS()){
				  cerr << "Error contain_orf " << contain_orf << " gene.containCDS() " << gene.containCDS()<< (*gene.getExons().begin())->start()<< " end exon " << (*gene.getExons().end())->end()<<endl;
			  }
			  if(gene.containCDS()) {
				  ostringstream oss;
				  oss << gene.getSeqname()<< "@" << gene.getCDS().first << "@" ;
				  string key = oss.str();
				  gene.keepModel(mapGeneModel);
			  }
		  }
		  before = time(NULL);
		  for(map<string,GeneModel >::iterator itMap = mapGeneModel.begin(); itMap != mapGeneModel.end() ; ++itMap){
			  itMap->second.printAnnot(ofs_genemodelsFilename, nbPath, formatDef);
			  nbPath++;
			  nb_models++;
		  }
		//  selectModel(listOfGene); //TODO do not create include model
		  after = time(NULL);
     }
      after = time(NULL);
      if(PRINTTIME) cout << " time to loop on cc " << difftime(after,before_cc) << endl;
    }
    if(verbose) {
      list<pair<s32,s32> > mypbCC;
      cerr <<"Number of usable connected component(s): " << nb_used_cc << endl;
      cerr <<"Number of gene model(s): " << nb_models << endl;
      for(map<string, list<s32> >::iterator itPE = probExons.begin();itPE != probExons.end(); itPE++){
    	  itPE->second.sort();
    	  itPE->second.unique();
    	  string tmp_cc_nb;
    	  if(itPE->second.front() != itPE->second.back())	mypbCC.push_back( make_pair(itPE->second.front(),itPE->second.back()) );
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
    if(itSeq == protPhase.end()){
      protPhase.insert( make_pair(key,phase) );
    }
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
    	  cout << "look at scaff " << name << endl;
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
  if(!name.empty()) {
	cout << "look at scaff " << name << endl;
    fastaDB.insert(make_pair(name,sequence));
  }
  fstrm.close();
  return fastaDB;
}

list< pair< GffRecord,string > > loadGff(char* filename,bool protein,map<string, map<s32,s32> >& cds,map<string,string>fastaDB){
	/* use it for reannotation
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
		 if(!tagRecord.empty() && tagRecord.back().first.tagValues[0] == ""){
			 cerr << " error in LoadGff : name empty (last column) "<< endl;
		 }

		 string idFasta;
		 assign(idFasta,record.ref);
		 map<string,string>::iterator itFasta = fastaDB.find(idFasta);
		 if(itFasta==fastaDB.end())continue;
		 else
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
				 if(!tagRecord.empty() && tagRecord.back().second == "d")
						 tagRecord.back().second = 'm';
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
		 cerr << " Error : we did not charge any data from "<< filename << endl;
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
	 int nbExon = 0;
	 while (!atEnd(gffIn)){
		readRecord(record, gffIn);
		string idFasta;
		assign(idFasta,record.ref);
		map<string,string>::iterator itFasta = fastaDB.find(idFasta);
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
	 if(tagRecord.back().second =="i") tagRecord.back().second ='f';
	 else tagRecord.back().second ='m';
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

			 if( cds.find(key) == cds.end()){
				 cds[key].insert(make_pair(phase,1));
			 }
			 else{
				 map<s32,s32>::iterator itCds = cds[key].find(phase);
				 if( itCds != cds[key].end()){
					 itCds->second =  itCds->second + 1;
				 }
				 else{
					cds[key].insert(make_pair(phase,1));
				 }
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
			 if( cds.find(key) == cds.end()){
				 cds[key].insert(make_pair(phase,1));
			 }
			 else{
				 map<s32,s32>::iterator itCds = cds[key].find(phase);
				 if( itCds != cds[key].end()){
					 itCds->second =  itCds->second + 1;
				 }
				 else{
					cds[key].insert(make_pair(phase,1));
				 }
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
		if(itCds->second.size() > 1) {
		}
		for(map<s32,s32>::iterator itPhase = itCds->second.begin() ; itPhase != itCds->second.end(); ++itPhase){
			if(itPhase->second >= max){
				secondMax = max;
				max = itPhase->second;
				maxPhase = itPhase->first;
			}
		}
		if(max == secondMax ){
			phase[itCds->first] = 0;
		}

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

void fusionMonoExons(list< pair < GffRecord,string> >& tagRecord){
	//Same idea as /env/ig/soft/rdbioseq/annotation-snapshot/linux-noarch/bin/loci
	time_t before,after;
	before = time(NULL);
	list< pair < GffRecord,string > > listMono = extractMono(tagRecord);

	after = time(NULL);
	before = time(NULL);
	listMono.sort(sortMonoGff);
	after = time(NULL);
	time_t beforeAll = time(NULL);
	std::clock_t c_start = std::clock();
	for(list<pair < GffRecord,string> >::iterator itListMono = listMono.begin(); itListMono != --listMono.end() ;){//--listMono.end() : we stop one element before the end, like size()-1
		list<pair < GffRecord,string> >::iterator itNext = itListMono;
		++itNext;
		if(itListMono->first.ref == itNext->first.ref && itListMono->first.strand == itNext->first.strand){
			if( (itNext->first.endPos >= itListMono->first.beginPos && itNext->first.endPos <= itListMono->first.endPos)
					|| (itListMono->first.endPos >= itNext->first.beginPos && itListMono->first.endPos <= itNext->first.endPos)){

				itNext->first.beginPos = min(itListMono->first.beginPos,itNext->first.beginPos);
				itNext->first.endPos = max(itListMono->first.endPos,itNext->first.endPos);
				itListMono = listMono.erase(itListMono);
			}
			else
				++itListMono;
		}
		else
			++itListMono;
	}
	after = time(NULL);
	std::clock_t c_end = std::clock();
	before = time(NULL);
	deleteIncludeMono(listMono,tagRecord);
	after = time(NULL);
	before = time(NULL);
	tagRecord.insert(tagRecord.end(),listMono.begin(),listMono.end());
	after = time(NULL);
}


void deleteIncludeMono(list< pair < GffRecord,string > >& listMono, list< pair < GffRecord,string> > tagRecord){
	//TODO before add at tagRecord, look if mono is inside a pluriexon : delete mono
	for(list< pair < GffRecord,string > >::iterator itMono = listMono.begin(); itMono != listMono.end();){
		bool boolDelete =false;
		for(list< pair < GffRecord,string > >::iterator itOther = tagRecord.begin(); itOther != tagRecord.end(); ++itOther){
			if(itMono->first.beginPos >= itOther->first.beginPos && itMono->first.endPos <= itOther->first.endPos){
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

bool sortMonoGff (const pair <GffRecord,string >& record1 , const pair <GffRecord,string>& record2){
	return (record1.first.beginPos <= record2.first.beginPos);
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
	cerr << "*    -f <file> : fasta file which contains genome sequence(s)." << endl << endl;
	cerr << "     -c <file> : tabular file which contains covtigs [sequence_name start_position end_position average_coverage]." << endl;
	cerr << "     -j <file> : tabular file which contains junctions [sequence_name end_position_of_exon1 start_position_of_exon2 strand]." << endl;
	cerr << "                 strand is 1 for forward and -1 for reverse; only four first fields are taken into account." << endl;
	cerr << "     -P <file> : tabular file which contains proteic alignment phase information [sequence_name base_coordinate phase]." << endl;
	cerr << "     --annot <file> : annotation file in gff "<<endl;
	cerr << "     --rna <file> : rna file in gff "<<endl;
	cerr << "     --prot <file> : prot file in gff"<<endl;

	cerr << "  OUTPUT FILES" << endl;
	cerr << "     -o <folder> : output folder, by default (./out) " << endl;
	cerr << "     -C <file> : output extended covtigs in the given file [sequence_name start_position end_position average_coverage]" << endl;
	cerr << "     -J <file> : output validated junctions in the given file [sequence_name start_position end_position strand .... ]" << endl;
	cerr << endl;
	cerr << "     -S        : do not output single exon models." << endl;
	cerr << "     -L <int>  : reduce the exon list size to -L parameters for each covtigs, default is 100.000." << endl;
	cerr << "     -x <int>  : size of regions where to find splice site around covtigs boundaries, default is 0." << endl;
	cerr << "                 if != 0, adds complexity to graph structure and to execution time." << endl;
	cerr << "                 only works with canonical splice sites."<< endl;
	cerr << "     -e <int>  : minimal size of exons, default is 25 nucleotides." << endl;
	cerr << "     -i <int>  : minimal size of introns, default is 9 nucleotides." << endl;
	cerr << "     -m <int>  : maximal size of introns, default is 50.000 nucleotides." << endl;
	cerr << "     -p <int>  : maximal number of paths inside a connected component, default is 10,000." << endl;
	cerr << "     -b <int>  : number of nucleotides around exons boundaries where to find start and stop codons, default is 30." << endl;
	cerr << "     -n <int>  : number of neighbour covtigs to test, default is 20." << endl;
	cerr << "     -t        : gtf format annotation file - default is gff3" << endl;
	cerr << "     -u        : choose model strand according to longest ORF - only works if input junctions are non-oriented" << endl;
	cerr << "     -v        : silent mode" << endl;
	cerr << "     -h        : this help" << endl;
	cerr << "--------------------------------------------------------------------------------------------" << endl;
	exit(1);
}
