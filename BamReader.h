/**
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 
 **/ 

#include <string>
#include <iostream>
#include <sstream>

#include "htslib/sam.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "samtools.h"

using namespace std;

class BamReader{
  samFile * in;
  bam1_t *b;
  hts_idx_t *idx; 
  bam_hdr_t *header;
  int nread;

public:
  ~BamReader(){
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
  }

  void setFile(string filename){
    if((in = sam_open(filename.c_str(), "r")) == 0) {
      cerr << "failed to open "<< filename << " for reading." << endl;
      exit(1);
    }
    if((header = sam_hdr_read(in)) == 0) {
      cerr << "failed to open "<< filename << " header for reading." << endl;
      exit(1);
    }
    //load the index
    idx = 0;
    idx = bam_index_load(filename.c_str());
    if(idx==0){
      cerr << "failed find index file for "<< filename << endl;
      exit(1);
    }
    b = bam_init1();
    nread=0;
  }
  
  string get_next_bad_map_seq(int flank_size, int min_gap){
    while(sam_read1(in, header, b) >= 0){
      nread++;
      if( nread % 1000000 == 0 ) cerr << nread/1000000 << " million reads processed" << endl;
      int seq_length=b->core.l_qseq;
      uint32_t *cigar = bam_get_cigar(b); //bam1_cigar(b);
      int matched=0;
      int rest=0;
      for(int k=0; k < b->core.n_cigar; k++){
	int c_oper=(cigar[k])&BAM_CIGAR_MASK;
	int c_size=(cigar[k])>>BAM_CIGAR_SHIFT;
	if(c_oper==BAM_CMATCH)
	  matched+=c_size;
	else
	  rest+=c_size;
      }
      if((seq_length-matched)<flank_size && rest < min_gap) continue ;

      //get the read sequence
      char qseq[seq_length];
      uint8_t * s = bam_get_seq(b); //bam1_seq(b);
      for(int n=0; n<seq_length; n++){
	char v = bam_seqi(s,n); //bam1_seqi(s,n);
	qseq[n] = seq_nt16_str[v]; //bam_nt16_rev_table[v];
      }
      return(string(qseq,seq_length));
    }
    return("");
  };

  bool is_current_read_unaligned(){
    return(b->core.flag & BAM_FUNMAP);
  };

  int get_coverage(string chrom, int pos){
    //have to do this formatting hack to get the chromosome index.
    stringstream formatted_pos;
    formatted_pos << chrom << ":" << pos << "-" << pos;
    int result=0;
    hts_itr_t *iter = sam_itr_querys(idx,header,formatted_pos.str().c_str());
    if(iter==NULL) return -1;
    while ( sam_itr_next(in, iter, b) >= 0) result++;
    hts_itr_destroy(iter);
    return result;
  }
} ;

