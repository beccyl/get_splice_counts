/**
 ** A class with utilities for reading bam files.
 ** 
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: Rebecca Evans <rebecca.louise.evans@gmail.com>
 ** Date Modified: 19 May 2020
 **/ 

#include <string>
#include <iostream>
#include <sstream>
#include <vector>

//depends on samtools header files
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>
#include <htslib/khash.h>
#include <samtools.h>

// macros adapted from samtools/stats.c 
// From the spec
// If 0x4 is set, no assumptions can be made about RNAME, POS, CIGAR, MAPQ, bits 0x2, 0x10, 0x100 and 0x800, and the bit  0x20 of the previous read in the template.
// 0x4 is BAM_FUNMAP, 0x2 is BAM_PROPER_PAIR, 0x10 is BAM_FREVERSE, 0x100 is BAM_FSECONDARY, 0x20 BAM_FMREVERSE
#define IS_PAIRED(bam) ((bam)->core.flag&BAM_FPAIRED)
#define IS_PAIRED_AND_MAPPED(bam) (((bam)->core.flag&BAM_FPAIRED) && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))
#define IS_PROPERLYPAIRED(bam) (((bam)->core.flag&(BAM_FPAIRED|BAM_FPROPER_PAIR)) == (BAM_FPAIRED|BAM_FPROPER_PAIR) && !((bam)->core.flag&BAM_FUNMAP) && !((bam)->core.flag&BAM_FMUNMAP))
#define IS_UNMAPPED(bam) ((bam)->core.flag&BAM_FUNMAP)
#define IS_REVERSE(bam) ((bam)->core.flag&BAM_FREVERSE)
#define IS_MATE_REVERSE(bam) ((bam)->core.flag&BAM_FMREVERSE)
#define IS_READ1(bam) ((bam)->core.flag&BAM_FREAD1)
#define IS_READ2(bam) ((bam)->core.flag&BAM_FREAD2)
#define IS_DUP(bam) ((bam)->core.flag&BAM_FDUP)
#define IS_ORIGINAL(bam) (((bam)->core.flag&(BAM_FSECONDARY|BAM_FSUPPLEMENTARY)) == 0)

// codes for determining which case arose to cause non-linear sequence
// withh be helpful for second stage of project
// ideally use syslog.h TODO

using namespace std;

class BamReader {
    samFile * in;
    bam1_t *b;
    hts_idx_t *idx; 
    bam_hdr_t *header;
    int nread;

    // copy constructor
    BamReader(const BamReader& other) { }

    // copy assignment
    BamReader& operator=(const BamReader& other)
    {
        if (this != &other)
        {
            *this = BamReader(other);
        }
        return *this;
    }

    // move constructor
    BamReader(BamReader && other) noexcept {}

    // move assignment
    BamReader& operator=(BamReader&& other) noexcept { return *this;}

    public:

    // constructor
    BamReader() { }

    ~BamReader() {
        // could also call "destroy()" 
        if (NULL != b) {
           bam_destroy1(b);
        }
        b = NULL;

        if (NULL != header) {
            bam_hdr_destroy(header);
        }
        header = NULL;

        if (NULL != idx) {
            hts_idx_destroy(idx);
        }
        idx = NULL;

        // close the samFile TODO: open in constructor, instead of setFile(?)
        if (NULL != in) {
            sam_close(in);
        }
        in = NULL;

        nread=0;
    }

    // TODO: ideally this should be a private function, and called by constructor
    bool setFile(string filename){
        if((in = sam_open(filename.c_str(), "r")) == 0) {
            cerr << "failed to open "<< filename << " for reading." << endl;
            return false;
        }
        if((header = sam_hdr_read(in)) == 0) {
            cerr << "failed to open "<< filename << " header for reading." << endl;
            sam_close(in); // close the samFile
            return false;
        }
        //load the index
        idx = 0;
        idx = bam_index_load(filename.c_str());
        if(idx==0){
            cerr << "failed find index file for "<< filename << endl;
            hts_idx_destroy(idx);
            sam_close(in); // close the samFile
            return false;
        }
        b = bam_init1();
        nread=0;
        return true;
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
    }

    bool is_current_read_unaligned(){
        return(b->core.flag & BAM_FUNMAP);
    }

    int get_coverage(string chrom, int pos1, int pos2){
        //have to do this formatting hack to get the chromosome index.
        stringstream formatted_pos;
        formatted_pos << chrom << ":" << pos1 << "-" << pos2;
        int result=0;
        hts_itr_t *iter = sam_itr_querys(idx,header,formatted_pos.str().c_str());
        if(iter==NULL) return -1;
        while ( sam_itr_next(in, iter, b) >= 0) result++;
        hts_itr_destroy(iter);
        return result;
    }

    vector< pair<int,int> >  get_bad_pair_positions(string chrom,int start,int end){
        const char* const dirn[] = { "+", "-" }; // 0 = forward, 1 = reverse
        vector< pair<int,int> > result;
        stringstream formatted_pos;
        formatted_pos << chrom << ":" <<  start << "-" << end;
        hts_itr_t *iter = sam_itr_querys(idx,header,formatted_pos.str().c_str());
        if(iter==NULL) return result;
        while(sam_itr_next(in, iter, b) >= 0){
            // continue if fqcfail or fdup ...
            //
            // if paired AND not a proper pair
            // (filter out single reads!, as can cause a core-dump)
            // BAM files may contain both single reads and paired-end reads. (they should have different RG?)
            // same direction read mappings, should not be a proper pair??? why are these occurring
            //
            if ( (b)->core.flag & BAM_FQCFAIL) { continue; } // skip if quality check failed
            if ( (b)->core.flag & BAM_FDUP) { continue; } // skip if pcr/optical duplicate
            if ( !((b)->core.flag & BAM_FPAIRED)) { continue; } // skip if not paired

            // TODO: test this use case, currently has not occurred on sampled data
            // IF both are set,  or both are unset -- read sam spec (currently skipping)
            if ( (((b)->core.flag & BAM_FREAD1) != 0) && (((b)->core.flag & BAM_FREAD2) != 0)) { continue; }
            if ( (((b)->core.flag & BAM_FREAD1) == 0) && (((b)->core.flag & BAM_FREAD2) == 0)) { continue; }

            // Capture "exceptional cases" in properly paired reads
            if (IS_PROPERLYPAIRED(b)) {
                // may occur if there is a supplementary (?chimeric) alignment
                // these may indicate a fusion(?)
                // From sam spec, if supplementary, inherit PROPER_PAIR flag from primary alignment
                if (bam_is_rev(b) == bam_is_mrev(b) || (b->core.tid != b->core.mtid)) {
                    //TODO:
                    cout << header->target_name[b->core.tid] << "\t" << b->core.pos
                        << "\t" << dirn[bam_is_rev(b)] << "\t"
                        << header->target_name[b->core.mtid] << "\t" << b->core.mpos
                        << "\t" << dirn[bam_is_mrev(b)] << "\t"
                        << "[" << (b)->core.flag << "]\t[1]"
                        << endl;
                } else if (bam_is_rev(b) && b->core.pos < b->core.mpos) {
                    // complementary orientation, same strand, proper pair
                    // check for outward facing -- can indicate possible duplication
                    cout << header->target_name[b->core.tid] << "\t" << b->core.pos
                       << "\t" << dirn[bam_is_rev(b)] << "\t"
                       << header->target_name[b->core.mtid] << "\t" << b->core.mpos
                       << "\t" << dirn[bam_is_mrev(b)] << "\t"
                       << "[" << b->core.flag << "]\t[2]"
                       << endl;
                }
            } else if (IS_PAIRED_AND_MAPPED(b) && !IS_PROPERLYPAIRED(b)) {
                bool diff_chrom = (b->core.tid != b->core.mtid);
                // if on different chromosomes, or same stranded -- non-linear
                if (diff_chrom || bam_is_rev(b) == bam_is_mrev(b)) {
                    // non linear if mapped to different chromosomes or same strand
                    cout << header->target_name[b->core.tid] << "\t" << b->core.pos
                        << "\t" << dirn[bam_is_rev(b)] << "\t"
                        << header->target_name[b->core.mtid] << "\t" << b->core.mpos
                        << "\t" << dirn[bam_is_mrev(b)]<< "\t"
                        << "[" << b->core.flag << "]\t[3]\t"
                        << endl;
                } else if (!bam_is_rev(b) && b->core.pos < b->core.mpos && b->core.mpos < end) {
                    // same chromosome and different strand (!(A||B) == !A && !B)
                    // case (C) forward facing downstream read, and mate within range
                    result.push_back( make_pair(b->core.pos+b->core.l_qseq,b->core.mpos));
                } else if (!bam_is_rev(b) || b->core.pos < b->core.mpos) {
                    // else  (not case (C)**) 
                    // either a forward facing read (infer must be upstream)**
                    // OR a reverse facing read which is downstream
                    // Essentially this means an outward facing pair.
                    cout << header->target_name[b->core.tid] << "\t" << b->core.pos
                        << "\t" << dirn[bam_is_rev(b)] << "\t"
                        << header->target_name[b->core.mtid] << "\t" << b->core.mpos
                        << "\t" << dirn[bam_is_mrev(b)]<< "\t"
                        << "[" << b->core.flag << "]\t[4]\t"
                        << endl;
                }
            }
        }
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

    int get_nreads(){
        return nread;
    }

    void destroy(){
        bam_destroy1(b);
        bam_hdr_destroy(header);
        sam_close(in);
    }

} ;
