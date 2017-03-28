# Gmove

# Dependencies
The script is written in c++, it uses the library boost (http://www.boost.org/) and the library seqAn (http://www.seqan.de/). You should install them and change the path in the Makefile, before compiling it.

# Running Gmove
`gmove -f <reference sequence> --rna <rna.gff> {Options}`

# Option

    gmove - Gene modelling using various evidence.

    Usage : gmove -f <reference sequence> --rna rna.gff {Options}
    !! Note : Arguments with * are required, the other are optionnal !!

    INPUT FILES
    *-f <file>          : fasta file which contains genome sequence(s).

     --annot <file>     : annotation file in gff
     --rna <file>               : rna file in gff
     --prot <file>      : prot file in gff
    OUTPUT FILES
     -o <folder>        : output folder, by default (./out)
     --raw              : output raw data

     -S                 : do not output single exon models.
     -e <int>   : minimal size of exons, default is 9 nucleotides.
     -i <int>   : minimal size of introns, default is 9 nucleotides.
     -m <int>   : maximal size of introns, default is 50.000 nucleotides.
     -p <int>   : maximal number of paths inside a connected component, default is 10,000.
     -b <int>   : number of nucleotides around exons boundaries where to find start and stop codons, default is 30.
     -t         : gtf format annotation file - default is gff3
     -u         : choose model strand according to longest ORF - only works if input junctions are non-oriented
     --cds         : min size CDS, by default 100
     -h         : this help
     --ratio                        : ratio CDS/UTR min 80% de CDS


# Input
Gmove reads **GFF2** and **GFF3** files. It recognizes some specific tags at column 3 : 
  - for files parse with the options `--rna` and `--prot`, it recognizes the tags **HSP** and **exon**
  - for files parse with the option `--annot`, it recognizes the tags **UTR** and **CDS**
  
  ## Preparing the files
    The transcript/protein's id (in the last column of the GFF files ) has to be uniq.

    You should sort your files : the exon with the smaller position on the genome should be at the first line of the transcript. In this way Gmove can reconstruct the transcript/protein by just reading once the files. So be careful, when the strand is "-", Gmove needs to read the last exon first for a corresponding transcript/protein. 
    
    for a GFF2 file : `sort -k1,1 -k10,10 -k4,4n rna.gff > sortRna.gff`
    
    for a GFF3 file :  `sort -k1,1 -k9,9 -k4,4n rna.gff > sortRna.gff`

# Output
The script will output a GFF or GTF file. 



# Installation
Install the dependencies

Download the git repository

Change the path in the Makefile

make gmove

# License
Gmove is distributed open-source under CeCILL FREE SOFTWARE LICENSE. Check out http://www.cecill.info/ for more information about the contents of this license.

# Contact
gmove [a] genoscope [.] cns [.] fr


# Acknowledgments
This work was financially supported by the Genoscope, Institut de Genomique, CEA and Agence Nationale de la Recherche (ANR), and France Génomique (ANR-10-INBS-09-08).
