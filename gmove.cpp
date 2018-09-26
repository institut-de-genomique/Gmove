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
#include <malloc.h> //TODO do we need all of those include ?

#include "SSRContig.h"
#include "SSRContigLists.h"
#include "NetEx.h"
#include "GeneModel.h"
#include "GeneModelList.h"
#include "GffRecord.h"
#include "GffRecordList.h"


using namespace std;
using namespace boost;

void usage();
void error(string msg);
map<string, string> loadFasta(char* filename);
bool existence(char* name);
void createDirectory(string,bool, ofstream &ofs_rawModelsFilename, ofstream & ofs_filterModelsFilename);
map<string,s32>  selectPhase(map<string,map<s32,s32> > cds);

int main(int argc, char** argv) {
  
  // default values of parameters
  s32 c, extend_for_splicesites=0, min_size_exon=9, min_size_intron=9, max_size_intron=50000, max_nb_paths=10000, search_window=30, nb_neighbour=20, verbose=0, min_size_cds=100;
  char *fastaFilename = NULL, *annotationFilename = NULL, *longReadsFilename = NULL, *rnaFilename = NULL, *dataProtFilename = NULL, *abinitioFilename = NULL;
  string out = "";
  bool formatDef=true, keep_single_exon_gene=true, unoriented_junctions=false, rawData = false, ratioCdsUtr = false;

  // options given as parameters
  const struct option longopts[] =
    {
      {"rna",   required_argument,        0, 'r'},
	  {"longReads", required_argument, 0, 'l'},
	  {"prot",   required_argument,        0, 'p'},
	  {"annot",   required_argument,        0, 'a'},
	  {"out",   required_argument,        0, 'o'},
	  {"abinitio",required_argument,0,'g'},
	  {"cds", required_argument, 0 ,'c'},
	  {"raw",no_argument ,0,'z'},
	  {"ratio",no_argument,0,'w'},

       {0,0,0,0},
    };
  int index;
  cerr << "Gmove (Gene MOdeling using Various Evidence)" << endl;

  cerr <<"Command line " ;
  for (int i = 0 ; i< argc ; i++){
	  cerr << argv[i] << " ";
  }
  cerr << endl;
  time_t  t = time(0);
  cerr << ctime(&t)<<endl;

  while ((c =  getopt_long(argc, argv, "f:c:a:r:l:d:o:g:S:x:e:i:m:p:P:b:l:t:u:vzwh",longopts,&index)) != -1) {//getopt(argc, argv, "f:c:j:a:D:d:G:C:J:SL:x:e:i:m:p:P:b:n:l:t:u:v:h")) != -1) {
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
    	rnaFilename = optarg;
    	break;
    case 'l':
    	longReadsFilename = optarg;
    	break;
    case 'p':
    	dataProtFilename = optarg;
    	break;
    case 'o':
    	out = optarg;
      break;
    case 'g':
    	abinitioFilename = optarg;
    	break;
    case 'S': 
      keep_single_exon_gene = false;
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
    case 'P':
      max_nb_paths = atoi(optarg);
      break;
    case 'b':
      search_window = atoi(optarg);
      break;
    case 't':
      formatDef=false;
      break;
    case 'v': 
      verbose = 1;
      break;
    case 'z':
       	rawData = true;
       	break;
    case 'w':
    	ratioCdsUtr = true;
    	break;
    case '?':
       default:
	   /* invalid option */
	   fprintf(stderr, "%s: option '-%c' is invalid: ignored\n",argv[0], optopt);
	   return EXIT_FAILURE;
    case 'h':
      usage();
    }
  }

  if(max_nb_paths < 1){
	  cerr << "invalid value for option -p " << max_nb_paths << endl;
	  return EXIT_FAILURE;
  }

  // general constants
  SSRContigList::NBNEIGHBOUR = nb_neighbour;
  SSRContigList::MINSIZEINTRON = min_size_intron;
  SSRContigList::MAXSIZEINTRON = max_size_intron;
  SSRContig::VERBOSE= verbose;
  SSRContig::EXTEND = extend_for_splicesites;
  SSRContig::MINSIZEEXON = min_size_exon;

  // OPTION FOR BLAT MODELS
  // Check input and output options
  if(fastaFilename == NULL || !existence(fastaFilename)) error("Input genomic sequence : empty filename or file not found.");
  if(rnaFilename != NULL && !existence(rnaFilename)) error("Input rna file (gff) : empty file or file not found.");
  if(longReadsFilename != NULL && !existence(longReadsFilename)) error("Input long reads file (gff) : empty file or file not found.");
  if(dataProtFilename != NULL && !existence(dataProtFilename)) error("Input prot file (gff) : empty file or file not found.");
  if(annotationFilename != NULL && !existence(annotationFilename)) error("Input annotation file (gff) : empty file or file not found.");
  if(abinitioFilename != NULL && !existence(abinitioFilename )) error("Input abinitio file (gff) : empty file or file not found.");


  if (longReadsFilename != NULL  && ( rnaFilename != NULL || dataProtFilename != NULL || annotationFilename != NULL)  ) error("Long reads data should be alone."); //
  cerr << "======Start reading input files " << endl;
  // create network
  list<NetEx*> listReseau;
  cerr << "\tLoad fasta file " << fastaFilename << "..." << endl;
  map<string, string> fastaDB = loadFasta(fastaFilename);

  // create output files
  //Load Rna seq Data
  time_t before = time(NULL);
  map<string,map<s32,s32> > cds;
  if(!existence(rnaFilename))
	  cerr <<"\tNo rna file" << endl;
  else
	  cerr << "\tLoad data file (gff) " << rnaFilename << "..."<<endl;

  GffRecordList gffRecord(rnaFilename,fastaDB,"rna");
  if(existence(rnaFilename))
	  cerr << "\tIt load " << gffRecord.getRecords()->size() << " exons from the rna file."<<endl;

  if(existence(longReadsFilename)){
	  GffRecordList gffRecordLong(longReadsFilename,fastaDB,"rna");
	  cerr << "\tIt load " << gffRecordLong.getRecords()->size() << " exons from the longReads file."<<endl;
	  gffRecord.copyGffRecordList(gffRecordLong);
  }

  if(existence(dataProtFilename)){
	  //Load protein data
	  cerr << "\tLoad protein (gff) "<< dataProtFilename<<"..."<<endl;
	  GffRecordList gffRecord2(dataProtFilename,fastaDB,"prot");
	  cerr << "\tIt load " << gffRecord2.getRecords()->size() << " exons from the proteins file."<<endl;
	  gffRecord2.loadGff(cds);
	  gffRecord.copyGffRecordList(gffRecord2); //TODO remove copy
  }
  else
	  cerr << "\tNo protein file" << endl;

  time_t after = time(NULL);
 if(SSRContig::VERBOSE) cout <<"time loadGff " << difftime(after,before) << endl;
  if(existence(annotationFilename)){
	 before = time(NULL);
	 cerr << "\tLoad annotation file (gff) " << annotationFilename << "..."<< endl;
	 GffRecordList gffRecord3(annotationFilename,fastaDB,"annot");
	 cerr << "\tIt load " << gffRecord3.getRecords()->size() << " exons from the annotation file."<<endl;
	 bool useCds = true;
	 gffRecord3.loadAnnotation(cds,useCds);
	 gffRecord.copyGffRecordList(gffRecord3);//TODO remove copy !
	 after = time(NULL);
	 if(SSRContig::VERBOSE) cout << "time load Annotation " << difftime(after,before) << endl;
  }
  else
	  cerr<<"\tNo annotation file"<<endl;

  if(existence(abinitioFilename)){
	  before = time(NULL);
	  cerr << "\tLoad ab initio file (gff) " << abinitioFilename << "..."<< endl;
	  cout << "load abinitio "<< endl;
	  GffRecordList gffRecord4(abinitioFilename,fastaDB,"abinitio");
	  cerr << "\tIt load " << gffRecord4.getRecords()->size() << " exons from the ab initio file."<<endl;
	  bool useCds = false;
	  gffRecord4.loadAnnotation(cds,useCds);
	  gffRecord.copyGffRecordList(gffRecord4);
	  after = time(NULL);
	  if(SSRContig::VERBOSE) cout << "time load ab initio " << difftime(after,before) << endl;
  }
  else
	  cerr <<"\tNo abinitio file"<< endl;

  map<string, s32 > proteinPhaseCoord = selectPhase(cds);
  before = time(NULL);
  ofstream ofs_rawModelsFilename, ofs_filterModelsFilename;
 if(gffRecord.empty()){
	  cerr<<"\tWarning : no data load. Stop the program"<<endl;
	  exit(0);
  }
  else{
	  createDirectory(out,rawData,ofs_rawModelsFilename, ofs_filterModelsFilename);
	  cerr << "\tIt load " << gffRecord.getRecords()->size() << " exons from all input files."<<endl;

  }

//cleanMono
 cerr << "\n======Mono exonique genes"<<endl;
 gffRecord.cleanMono();

 after = time(NULL);
 if(verbose)  cout <<"time fusionMonoExons " << difftime(after,before) << endl;

	  cerr << "\n======Building graph" << endl;
  gffRecord.intron();
  SSRContigLists listeContigs(gffRecord,fastaDB);

  before = time(NULL);
  listReseau = listeContigs.buildGraph();
  after = time(NULL);
  if(verbose)  cout << " time building graph " << difftime(after,before)<<endl;

  list<NetEx*>::iterator itNet;
   s32 cb=0;
   cerr << "There are " << listReseau.size() << " to annotate " << endl;
  for (itNet = listReseau.begin() ; itNet != listReseau.end() ; itNet++) {
    NetEx* oneRes = *itNet;
    s32 num = std::distance(listReseau.begin(), itNet);
    ++num;
    cerr << "Start Scaffold " << oneRes->getSeqName()<< " scaff " << num << "/" << listReseau.size() << endl;
    before = time(NULL);
    cout << "out " << out << endl;
    oneRes->graphOut(out,"beforeClean"); //+"dot_files/"
    cerr << "\n======Clean graph" << endl;
    oneRes->cleanGraph();


    oneRes->graphOut(out,"afterClean"); //+"dot_files/"
    after = time(NULL);
    if(verbose)  cout << " time cleaning graph " << difftime(after,before) << endl;

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
    	cout<<endl;
    	cerr << "Change cc "<< c+1 <<"/"<< components.size()<<endl;
    	cout << "Change cc "<< c+1 <<"/"<< components.size()<<endl;
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
    		cerr << " nbPaths " << nbPaths <<" size component "<< component.size() << " vertices"<< endl;
    		// calculate coordinates of too furnished region
    		s32 threshold =1;

    //		while(nbPaths > max_nb_paths || nbPaths == -1) { //TODO it is not the same with pathFinderWithCondition
    		/*			  s32 begReg = -1, endReg = 0;
				  for(list<s32>::const_iterator itComp = component.begin(); itComp != component.end(); itComp++){
					  TSSRList::iterator it1 = oneRes->getVertices()->begin();
					  advance(it1,*itComp);
					  SSRContig* ctgTag = *it1;
					  if (begReg == -1 || begReg > ctgTag->start() ) begReg = ctgTag->start();
					  	  if (ctgTag->end() > endReg) endReg = ctgTag->end();
				  }
				  cerr <<"\t[Warning] Too many paths ( " << nbPaths << " ) : " << oneRes->getSeqName() << "\t" << begReg << "\t" << endReg << endl;
				  cerr <<"Simplify big graph by removing intron with occurences of " << threshold << endl;
		//		  oneRes->simplifyBigGraph(component,threshold);
			  oneRes->graphOut(out,"gifGraph"); //+"dot_files/"
				  components = oneRes->getComponents();
				  component = components[c];//replace former component with the new simplify one
				  nbPaths = oneRes->count_allPaths(component);
				  ++threshold;
	*/
	//	  }
				   oneRes->synchronisedId();
/*		  before = time(NULL);


		   s32 begReg = -1, endReg = 0;
		  for(list<s32>::const_iterator itComp = component.begin(); itComp != component.end(); itComp++){
			  TSSRList::iterator it1 = oneRes->getVertices()->begin();
			  advance(it1,*itComp);
			  SSRContig* ctgTag = *it1;
			  if (begReg == -1 || begReg > ctgTag->start() ) begReg = ctgTag->start();
				  if (ctgTag->end() > endReg) endReg = ctgTag->end();
		  }
		//  cerr <<"\t[Warning] Too many paths ( " << nbPaths << " ) : " << oneRes->getSeqName() << "\t" << begReg << "\t" << endReg << endl;

*/
		  if(longReadsFilename != NULL)
			 allpaths= oneRes->PathsFinderWithCondition(component);

		  else if(nbPaths > max_nb_paths || nbPaths == -1){
			  s32 begReg = -1, endReg = 0;
							  for(list<s32>::const_iterator itComp = component.begin(); itComp != component.end(); itComp++){
								  TSSRList::iterator it1 = oneRes->getVertices()->begin();
								  advance(it1,*itComp);
								  SSRContig* ctgTag = *it1;
								  if (begReg == -1 || begReg > ctgTag->start() ) begReg = ctgTag->start();
								  	  if (ctgTag->end() > endReg) endReg = ctgTag->end();
							  }
							  cerr <<"\t[Warning] Too many paths ( " << nbPaths << " ) : " << oneRes->getSeqName() << "\t" << begReg << "\t" << endReg << endl;

			 cout << "too many path"<< endl;
			  allpaths =  oneRes->pathWithHigherWeight( component);
		  }

		  else
			  allpaths = oneRes->allPathsFinder(component);

		  after = time(NULL);
		  if(verbose) cout << " time allPathFinder " << difftime(after,before)<< " number of paths " << allpaths.size() << endl;
		  nb_used_cc++;
		  s32 nbPath=1;
		  // find a model for each path in the connected component
		  before = time(NULL);
		  MapModel mapGeneModel;
		  cerr<<" ======Compute models "<<endl;
		  for(list<list<s32> >::iterator itPa = allpaths.begin(); itPa != allpaths.end(); itPa++) {
			  cerr << "Look at path " << distance(allpaths.begin(), itPa)<< "/" << allpaths.size() << endl;

			  GeneModel gene = GeneModel(*itPa, oneRes->getVertices(), oneRes->getSeqName(), c+1, search_window, probExons);//XXX What is probExons ?
			  std::clock_t c_start = std::clock();
			  bool contain_orf = gene.findORF(proteinPhaseCoord);
			  std::clock_t c_end = std::clock();
			  if(verbose) cout << " time findORF " <<  1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << endl;
			  // On imprime une seule CDS par transcrit
			  if(contain_orf != gene.containCDS()){
				  cerr << "Error contain_orf " << contain_orf << " gene.containCDS() " << gene.containCDS()<< gene.getCDS().first<< " end CDS " <<  gene.getCDS().second <<endl;
			//	  exit(1);
			  }
			  if(gene.containCDS()) {
				  gene.keepModel(mapGeneModel); //TODO remove this !!!
			  }
		  }


		  before = time(NULL);
		  vector<GeneModel> vectorOfGene;
		  nbPath=1;

		  //FIXME pourquoi faire tout ça ??? ***************************************
		  for(MapModel::iterator itMap = mapGeneModel.begin(); itMap != mapGeneModel.end(); ++itMap){
			  itMap->second.setID(nbPath);
			  nbPath++;
			  nb_models++;
		  }
		  cerr<<"\tClean models"<<endl;
		  std::clock_t beforeClean = std::clock();
		  GeneModelList listOfGene(mapGeneModel);//FIXME another way to convert map to list ? splice ??
		  std::clock_t afterClean = std::clock();

		  //beforeClean = time(NULL);
		  listOfGene.deleteIncludedModel();
		  if(rawData)
			  s32 nb_models_raw = listOfGene.printOut(ofs_rawModelsFilename,formatDef);



		 // afterClean = time(NULL);
		//  if(verbose)  cout << " time deleteIncludedModel " << 1000.0 * (afterClean-beforeClean) / CLOCKS_PER_SEC  << endl;

		 // beforeClean = time(NULL);
	  	  if(listOfGene.getSize() > 1)
			  listOfGene.longestCDS();
		//  afterClean = time(NULL);
		 // if(verbose)   cout << " time deleteSmallCDS and print annot " << 1000.0 * (afterClean-beforeClean) / CLOCKS_PER_SEC << endl;
		  allModelFromOneRes.insertModels(listOfGene); //FIXME What is it for ? copy de list, pb ! mapGeneModel -> listOfGene -> allModelFromOneRes Pourquoi listOfGene existe ?

    	}
      after = time(NULL);
      if(verbose)  cout << " time to loop on cc " << difftime(after,before_cc) << endl;

    }
    if(verbose) cout << "allModelfromOneRes " << allModelFromOneRes.getSize() << endl;

  	allModelFromOneRes.filter(ratioCdsUtr,min_size_cds, longReadsFilename);
    nb_models_filter = allModelFromOneRes.printOut(ofs_filterModelsFilename, formatDef);


    cerr <<"Number of usable connected component(s): " << nb_used_cc << endl;
      cerr <<"Number of gene model(s): " << nb_models << endl;
      cerr <<"Number of gene model(s) after filtering : " << nb_models_filter << endl;

     /*
      list<pair<s32,s32> > mypbCC; //TODO what is mypbCC ?
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
*/
    delete oneRes;
  }
  t = time(0);
  cerr <<"Gmove finish at " <<  ctime(&t)<<endl;
  ofs_rawModelsFilename.close();
  ofs_filterModelsFilename.close();
  return 0;
}


//---------------------------------------------------------------------------------------------------
// Functions
//---------------------------------------------------------------------------------------------------


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
    	  fastaDB.insert(make_pair(name,sequence));
    	  sequence.erase();
      }
      seqline.erase(0,1);
      name= seqline;
      fstrm.getline(ch,1000);
    }
    else {
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


void createDirectory(string out, bool rawData,ofstream& ofs_rawModelsFilename,ofstream& ofs_filterModelsFilename){
	string outRaw, outFilter;
	if(out.empty()){
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
		std::srand(std::time(0)); // use current time as seed for random generator
		outRaw += "tmpfileXXXXXX";
		char *tmpnameRaw = strdup(outRaw.c_str());
		mkstemp(tmpnameRaw);
		ofs_rawModelsFilename.open(tmpnameRaw, std::ofstream::out);
	}

	outFilter =out;
	outFilter += "filter/gmove_";
	out +=  "gmove_";
	std::srand(std::time(0)); // use current time as seed for random generator
	out += "tmpfileXXXXXX";
	string randName = "tmpfileXXXXXX";
	outFilter += randName;
	char *tmpnameFilter = strdup(outFilter.c_str());
	mkstemp(tmpnameFilter);
	ofs_filterModelsFilename.open(tmpnameFilter, std::ofstream::out);

	//create folder for dot_files
	string cmd = "mkdir -p ";
	cmd += out;
	cmd += "dot_files/";
	system(cmd.c_str());

}


void usage() {
	cerr << "--------------------------------------------------------------------------------------------" << endl;

	cerr << "gmove - Gene modelling using various evidence." << endl << endl;
	cerr << "Usage : gmove -f <reference sequence> {Options}" << endl;
	cerr << "INPUT FILES" << endl;
	cerr << "\t*-f <file>\t\t: fasta file which contains genome sequence(s)." << endl ;
	cerr << "\t--rna <file>\t\t: rna file in gff (expect tag 'exon' or 'HSP' in column 3)"<<endl;
	cerr << "\t--longReads <file>\t\t: rna file in gff from long reads sequencing (expect tag 'exon' or 'HSP' in column 3)"<<endl;
	cerr << "\t--prot <file>\t\t: prot file in gff (expect tag 'exon' or 'HSP' in column 3)"<<endl;
	cerr << "\t--annot <file>\t\t: annotation file in gff (expect tag 'CDS' or 'UTR' in column 3)"<<endl;
	cerr << "\t--abinitio <file>\t: ab initio file in gff (expect tag 'CDS' or 'UTR' in column 3)"<<endl;

	cerr <<endl;

	cerr << "OUTPUT FILES" << endl;
	cerr << "\t-o <folder>\t\t: output folder, by default ./out " << endl;
	cerr << "\t--raw\t\t\t: output raw data " << endl;

	cerr << endl;

	cerr << "\t--ratio\t\t\t: ratio CDS/UTR min 80% de CDS" << endl;
	cerr << "\t-S\t\t\t: do not output single exon models." << endl;
	cerr << "\t-e <int>\t\t: minimal size of exons, default is 9 nucleotides." << endl;
	cerr << "\t-i <int>\t\t: minimal size of introns, default is 9 nucleotides." << endl;
	cerr << "\t-m <int>\t\t: maximal size of introns, default is 50.000 nucleotides." << endl;
	cerr << "\t-P <int>\t\t: maximal number of paths inside a connected component, default is 10,000." << endl;
	cerr << "\t-x <int>\t\t: size of regions where to find splice site around covtigs boundaries, default is 0." << endl;
	cerr << "\t-b <int>\t\t: number of nucleotides around exons boundaries where to find start and stop codons, default is 30." << endl;
	cerr << "\t-t\t\t\t: gtf format annotation file - default is gff3" << endl;
//	cerr << "\t-u\t\t: choose model strand according to longest ORF - only works if input junctions are non-oriented" << endl;
	cerr << "\t--cds\t\t\t: min size CDS, by default 100 " << endl;
	cerr << "\t-v\t\t\t: verbose" << endl;
	cerr << "\t-h\t\t\t: this help" << endl;

	cerr << "--------------------------------------------------------------------------------------------" << endl;
	exit(1);
}
