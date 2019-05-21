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

  pair<int, int> get_allele_depth(string & chrom, int & pos){
    stringstream formatted_pos;
    formatted_pos << chrom << ":" << pos << "-" << pos;
    hts_itr_t *iter = sam_itr_querys(idx,header,formatted_pos.str().c_str());
    if(iter==NULL) return make_pair(-1,-1);
    int tally[4]={0};
    while ( sam_itr_next(in, iter, b) >= 0){
      //first work out the offset between the start of the mapped position
      //and the SNP position
      uint32_t map_pos=b->core.pos;
      int diff=pos-map_pos;
      //now get the cigar string and work out which base in the read we need
      //to look at.
      uint32_t *cigar = bam_get_cigar(b);
      int read_offset=0, k=0;
      stringstream cig;
      while( k < b->core.n_cigar & diff>=0){
        int c_oper=(cigar[k])&BAM_CIGAR_MASK;
        int c_size=(cigar[k])>>BAM_CIGAR_SHIFT;
	cig << c_size ;
	switch(c_oper)
	  {
	  case BAM_CMATCH: {
	    diff-=c_size; cig << "M" ;
	    read_offset+=c_size;
	    break;
	  }
	  case BAM_CINS:{ cig << "I" ;
	      read_offset+=c_size; break ;}
	  case BAM_CDEL: { cig << "D" ;
	      diff-=c_size; break;}
	  case BAM_CREF_SKIP: { cig << "N" ;
	      diff-=c_size; break;}
	  case BAM_CSOFT_CLIP: { cig << "S" ;
	      read_offset+=c_size; break;}
	  default:
	    break;
	  }
	k++;
      }
      read_offset+=diff; //final adjustment
      if(read_offset>b->core.l_qseq | read_offset<0 | diff>0)
	continue;
      /**      cout << "  Pos="<<pos<<" MapPos=" << map_pos ;
      cout << " Diff="<<pos-map_pos;
      cout << " ReadOffset="<<read_offset ;
      cout << " ReadLength="<<b->core.l_qseq;
      cout << " Cigar=" << cig.str();**/
      uint8_t * seq = bam_get_seq(b); //bam1_seq(b);
      char base = bam_seqi(seq,read_offset); //bam1_seqi(s,n);
      tally[seq_nt16_int[base]]++;
      //      cout << " " << seq_nt16_str[base] << endl;
    }
    pair<int,int> allele_count(0,0);
    for(int t=0; t<4 ; t++){
      if(tally[t]>allele_count.first){
	allele_count.second=allele_count.first;
	allele_count.first=tally[t];
      } else if ( tally[t]>allele_count.second ){
	allele_count.second=tally[t];
      }
      //      cout << t << ":" << tally[t] << " ";
    }
    //    cout << endl;
    hts_itr_destroy(iter);
    return allele_count;
  }

  void destroy(){
    bam_destroy1(b);
    bam_hdr_destroy(header);
    sam_close(in);
  }

} ;

