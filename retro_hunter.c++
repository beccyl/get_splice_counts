/**
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 
 **/ 

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <experimental/string_view>
#include <algorithm>

#include "BamReader.h"
#include "ReferenceReader.h"

//#include <gperftools/profiler.h>
using namespace std;
using namespace std::experimental;

static const bool ALLOW_MISMATCH=false;
static const int MIN_GAP=50; //200000;
static const int MIN_COUNTS=2;

static int n_first_match=0;
static int n_perfect_match=0;

static BamReader bam_reader;
static ReferenceReader ref_reader;

//unordered_map<string_view,string>::iterator find_non_exact_match(string_view & orig_kmer,JunctionSeq& ref_reader);

//complimenting code taken from stackoverflow
char compliment(char& c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

class RetroJunction {
public:
  string gene;
  string chrom;
  int split_reads;
  int split_pairs;
  pair<int,int> position;
  pair<int,int> exon_coverage;
  pair<int,int> intron_coverage;
  void print(ofstream & ofs){
    ofs << chrom << ":" << position.first << "-" << position.second 
	<< "," << split_reads << "," << split_pairs
	<< "," << intron_coverage.first << "/" << exon_coverage.first 
	<< "," << intron_coverage.second << "/" << exon_coverage.second;
      }
};




static class Counts {
  unordered_map< string, pair<int,int> > _counts;

  void get_interesting_junctions(string name, int & read_support, int & pair_support, 
				 unordered_set<string> & black_list,
				 unordered_map< string, vector < RetroJunction > > & retro_results
				 ){
    //get positions and genes from the exon junction ids
    stringstream ss(name);
    string field;
    vector<string> gene_info;
    while(getline(ss,field, EXON_ID_DELIM)) {
      gene_info.push_back(field);
    } //in case the format looks wrong.
    if(gene_info.size()!=8){
      cerr << "Issue with exon sequence IDs." 
	   << "Format should be Gene:Chrom:Start:END/START(+/-)" << endl;
      exit(1);
    }

    //otherwise fill in the gene info
    vector<string> gene{gene_info.at(0),gene_info.at(4)};
    vector<string> chrom{gene_info.at(1),gene_info.at(5)};
    vector<int>pos{atoi(gene_info.at(2).c_str()), 
	atoi(gene_info.at(6).c_str())};
    vector<char> strand{gene_info.at(3).at(gene_info.at(3).size()-2),
	gene_info.at(7).at(gene_info.at(7).size()-2)}; //get the second last character
    
    bool enough_support = (read_support >= MIN_COUNTS) && (pair_support >= MIN_COUNTS);
    stringstream junc_pos_formatted;
    junc_pos_formatted << gene[0] << "\t" << chrom[0] << "\t" << pos[0] << "\t"
		       << gene[1] << "\t" << chrom[1] << "\t" << pos[1] ;
    bool not_in_black_list = black_list.find(junc_pos_formatted.str())==black_list.end();
      
    int idx1=0; int idx2=1;
    if( strand[0]=='-') {idx1=1 ; idx2=0; }

    if( gene[0]==gene[1] & enough_support & not_in_black_list){
      RetroJunction rj;
      rj.gene=gene[0];
      rj.chrom=chrom[0];
      rj.position=make_pair(pos[idx1],pos[idx2]);
      rj.split_reads = read_support;
      rj.split_pairs = pair_support;
      rj.intron_coverage=make_pair(
				   bam_reader.get_coverage(chrom[idx1],pos[idx1]+4,pos[idx1]+4),//FIX
				   bam_reader.get_coverage(chrom[idx2],pos[idx2]-4,pos[idx2]-4)
				   );
      rj.exon_coverage=make_pair(
				 bam_reader.get_coverage(chrom[idx1],pos[idx1]-1,pos[idx1]-1),//FIX
				 bam_reader.get_coverage(chrom[idx2],pos[idx2]+1,pos[idx2]+1)
				 );
      
      retro_results[gene[0]].push_back(rj);
    }
  };
  
public:

  vector<string> get_gene_short_list(){
    //return the names of all genes with junctions
    //loop over the counts and get gene names
    vector<string> genes;
    unordered_map<string, pair<int,int> >::iterator counts_itr=_counts.begin();
    for(;counts_itr!=_counts.end(); counts_itr++){
      //if the event is not in the black list check if it's interesting..
      size_t pos = counts_itr->first.find_first_of(EXON_ID_DELIM);
      genes.push_back(counts_itr->first.substr(0,pos));
    }
    sort( genes.begin(), genes.end() );
    genes.erase( unique( genes.begin(), genes.end() ), genes.end() );
    return genes;
  };

  void increment_split_read_count(string& end,string& start){
    string pair = end + EXON_ID_DELIM + start;
    _counts[pair].first++;
  };

  void increment_split_pair_count(string& pair){
    _counts[pair].second++;
  };

  void print_table(unordered_set<string> & black_list, string outfile ){
    cerr << "Reporting counts..."<< endl;

    //loop through junction counts and apply some filtering
    unordered_map<string,pair<int,int>>::iterator counts_itr=_counts.begin();
    unordered_map< string, vector < RetroJunction > > retro_results ;
    for(;counts_itr!=_counts.end(); counts_itr++){
      //if the event is not in the black list check if it's interesting..
      get_interesting_junctions(counts_itr->first,
				counts_itr->second.first,
				counts_itr->second.second,
				black_list,retro_results);
    }
    ofstream ofs(outfile);
    unordered_map< string , vector<RetroJunction> >::iterator ri;  

    for(ri=retro_results.begin() ; ri !=retro_results.end() ; ri++){
      ofs << ri->first << "\t" ;
      vector<RetroJunction> rj=ri->second;
      ofs << rj.size() << "\t" ;
      
      //sort the junctions by genomic position
      sort(rj.begin(), rj.end(), [=](RetroJunction& a, RetroJunction& b){
	return a.position.first < b.position.first;
      });      

      //get coverage over introns and exons
      int total_read_support=0;
      for(int j=0; j< rj.size() ; j++ )
	total_read_support += rj.at(j).split_reads + rj.at(j).split_pairs ;
      ofs << total_read_support << "\t" << bam_reader.get_nreads() << "\t" ;

      for(int j=0; j< rj.size() ; j++ ){
	rj.at(j).print(ofs);
	ofs << "\t";
      }
      ofs << endl;
    }
    ofs.close();
  };

} counts;

//static unordered_map< string, vector < pair <int, int > > > exon_positions;


unordered_map<string_view,string>::iterator 
find_non_exact_match(string_view & orig_kmer,string type){
  string kmer(orig_kmer);
  for(int base=0; base < FLANK_SIZE ; base++){
    vector<string> nuc{"A","G","C","T"};
    for(int n=0; n< nuc.size(); n++){
      kmer.replace(base,1,nuc.at(n));
      string_view sv_kmer=kmer;
      unordered_map<string_view,string>::iterator match=ref_reader.find(sv_kmer,type);
      if(match!=ref_reader.end(type))
	return match;
    }
  }
  return ref_reader.end(type);
}

//loop through the sequence to find a match
bool get_match(string & seq){
  if(seq.size()<(2*FLANK_SIZE)) return false;
  unordered_map<string_view,string>::iterator end; //end of exon1
  unordered_map<string_view,string>::iterator start; //joins to start of exon2
  string_view sv_seq = seq;
  string_view kmer1,kmer2;
  //search in the forward direction
  for(int pos=0; pos < (seq.size()-FLANK_SIZE) ; pos++){
    kmer1=sv_seq.substr(pos,FLANK_SIZE);
    end=ref_reader.find(kmer1,END_LABEL);
    //if a match is found. look for other side of the junction
    if(end!=ref_reader.end(END_LABEL)){
      n_first_match++;
      kmer2=sv_seq.substr(pos+FLANK_SIZE,FLANK_SIZE);
      start=ref_reader.find(kmer2,START_LABEL);
      //if other end is found
      if(start!=ref_reader.end(START_LABEL)){
	n_perfect_match++;
	//cout << "FOUND perfect match" << endl;
	counts.increment_split_read_count(end->second,start->second);
	return true;
      }
      //start not found. Try permutating the bases to account for 1 mismatch
      if(ALLOW_MISMATCH){
	start=find_non_exact_match(kmer2,START_LABEL);
	if(start!=ref_reader.end(START_LABEL)){
	  counts.increment_split_read_count(end->second,start->second);
	  return true;
	}
      }
    }
  }
  if(ALLOW_MISMATCH){
    //check again in reverse, permutating the end bases:
    for(int pos=seq.size()-FLANK_SIZE-1; pos >= FLANK_SIZE ; pos--){
      kmer1=sv_seq.substr(pos,FLANK_SIZE);
      start=ref_reader.find(kmer1,START_LABEL);
      //if a match is found. look for other side of the junction
      if(start!=ref_reader.end(START_LABEL)){
	kmer2=sv_seq.substr(pos-FLANK_SIZE,FLANK_SIZE);
	end=find_non_exact_match(kmer2,END_LABEL);
	if(end!=ref_reader.end(END_LABEL)){
	  counts.increment_split_read_count(end->second,start->second);
	  return true;
	}
      }
    }
  }
  return false;
}

int main(int argc, char *argv[]){
  if(!(argc==4 | argc==5)){
    cerr << "Usage: retro_hunter <exon_flanking_seq.fasta> <output prefix> <in.bam> [black_list]" << endl;
    exit(1);
  }
  string flank_fasta=argv[1];
  string in_filename=argv[3];
  string out_prefix=argv[2];
  string fusions_out_file=out_prefix+".ret";

  unordered_set<string> black_list;
  if(argc==5){ // read the black list
    ifstream black_list_stream;
    black_list_stream.open(argv[4]);
    if(!(black_list_stream.good())){ //check it opens
      cout << "Unable to open file " << argv[4] << endl;
      exit(1);
    }
    string line;
    while(getline (black_list_stream,line)){
      black_list.insert(line);
    }
  }
  
  //Read the fasta reference file
  cerr << "Reading fasta file of junction sequences: " << flank_fasta << endl;
  unordered_map<string,string> seqs;
  ifstream file;
  file.open(flank_fasta);
  if(!(file.good())){ //check it opens
    cout << "Unable to open file " << flank_fasta << endl;
    exit(1);
  } 
  string id="";
  string line;
  while ( getline (file,line) ){
    int start=line.find(">")+1;
    if(start==1){ //if this is the ID line...
        int end=line.find_first_of("\t\n ")-1;
        id=line.substr(start,end);
    } else {
      seqs[id]=seqs[id]+line;
    }
  }
  ref_reader.read_fasta(seqs);
  cerr << "Done reading fasta" << endl;
  //ProfilerStart("prof.out");

  //Read the bam file (using htslib API)
  bam_reader.setFile(in_filename);  

  int r;
  int i=0;
  int nread_processed=0;
  int f_count=0;
  int r_count=0;
  int nread=0;
  int unmapped=0;

  string seq;
  while((seq=bam_reader.get_next_bad_map_seq(FLANK_SIZE,MIN_GAP))!=""){
    nread_processed++;
    //loop through the sequence to find the first match
    if(get_match(seq)){
      f_count++;
      //also check the reverse compliment if the read is unmapped.
    } else { //if (bam_reader.is_current_read_unaligned()) {
      //unmapped++;
      reverse(seq.begin(),seq.end());
      transform(seq.begin(),seq.end(),seq.begin(),compliment);
      if(get_match(seq)){
	r_count++;
	//	cout << seq << endl;
      }
    }
    //find all matches
  }

  vector<string> genes = counts.get_gene_short_list();
  for(int g=0; g<genes.size(); g++){
    int start, end;
    string chrom;
    ref_reader.get_gene_coords(genes.at(g),chrom,start,end);
    cout << genes.at(g) << " " << chrom << ":" << start << "-" << end << endl;
    
    vector< pair< int, int> > pairs=bam_reader.get_bad_pair_positions(chrom,start,end);
    for(int i=0; i< pairs.size(); i++){
      int jend,jstart;
      string junc_string=ref_reader.get_closest_junction(genes.at(g), chrom, pairs.at(i).first, pairs.at(i).second);
      if(junc_string!="")
	counts.increment_split_pair_count(junc_string);
      //	cout << junc_string << endl;
    }
  }

  //  ProfilerStop();
  /**
  cerr << "Reads Total=" << nread << endl;
  cerr << "Reads Processed=" << nread_processed << endl;
  cerr << "One match=" << n_first_match << endl;
  cerr << "Junction counts=" << f_count << "   " << r_count << endl;
  cerr << "Perfect matches=" << n_perfect_match << endl;
  cerr << "Unmapped=" << unmapped << endl;
  **/

  //print out the table of counts
  //bam_reader.setFile(in_filename);
  counts.print_table(black_list,fusions_out_file);

  //Get the SNPs
  file.close();

}
