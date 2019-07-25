/**
 ** Author: Nadia Davidson, nadia.davidson@mcri.edu.au
 ** Modified: 
 **/ 

#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <experimental/string_view>

using namespace std;
using namespace std::experimental;

static const int FLANK_SIZE=20;
static const string END_LABEL=":END";
static const string START_LABEL=":START";
static const char EXON_ID_DELIM=':';

class ExonPositions {
  string _chrom="";
  char _strand='?'; 
  vector<int> starts;
  vector<int> ends;

 public:

  /** check that this isn't a duplicate genes and add exon boundary to database */
  bool add_boundary(string chrom, int pos, string type, char strand){
    if(_chrom=="") _chrom=chrom;
    if(_strand=='?') _strand=strand;
    else if (_chrom!=chrom | _strand!=strand)
      return false;
    if(type==START_LABEL)
      starts.push_back(pos);
    else
      ends.push_back(pos);
    return true;
  };
  
  /** get the start and end position of the gene */
  void get_position(string & chrom, int & start, int &end){
    if(starts.size()==0 | ends.size()==0){
      cerr << "Gene with no exons?" << endl; exit(1);
    }
    chrom=_chrom;
    start=*min_element(starts.begin(),starts.end());
    end=*max_element(ends.begin(),ends.end());
  };

  /** get which junction is spanned by paired end reads
      return false if pair doesn't span a junction */
  bool get_closest_junction(string chrom, int r1_end, int r2_start, int &jend, int&jstart){
    //    cout << r1_end << "\t" << r2_start << endl;
    if(_chrom!=chrom) return false;
    //find the closes end
    jend=-1;
    for(int i=0; i< ends.size(); i++){
      int best_diff=jend-r1_end;
      int this_diff=ends.at(i)-r1_end;
      if(this_diff>0 & (this_diff<best_diff | best_diff<0) )
	jend=ends.at(i);
    }
    jstart=-1;
    for(int i=0; i< starts.size(); i++){
      int best_diff=r2_start-jstart;
      int this_diff=r2_start-starts.at(i);
      if(this_diff>0 & (this_diff<best_diff | best_diff<0) )
	jstart=starts.at(i);
    }
    if(jend<0 | jstart<0) return false;
    if(jstart<jend && _strand=='+') return false;
    if(jend<jstart && _strand=='-') return false;
    else return true;
  };

};

class ReferenceReader { //read the fasta
  unordered_map<string,ExonPositions> ref_pos;
  unordered_map<string_view,string > junc_seq_start;
  unordered_map<string_view,string > junc_seq_end;

  void read_fasta( unordered_map<string,string> & _seqs, const string type){
    unordered_map<string_view,string > junc_seq;
    //loop through the sequences and sort into start or end flanking sequence
    //mark any duplicate sequences for later removal
    vector<string_view> to_erase; //list of junction sequences that aren't unique.
    unordered_map<string,string>::iterator seq_itr=_seqs.begin();
    for(; seq_itr!=_seqs.end(); seq_itr++){
      //if more than one junction with this sequence
      //will need to remove later.
      string id=seq_itr->first;
      if(seq_itr->second.size()<FLANK_SIZE){
	cerr << "Reference sequence too short for flank length" << endl;
	exit(1);
      }
      bool is_type = id.find(type)!=string::npos;
      
      if(is_type){// & right_size){
	if(type==START_LABEL)
	  seq_itr->second=seq_itr->second.substr(0,FLANK_SIZE);
	else
	  seq_itr->second=seq_itr->second.substr(seq_itr->second.size()-FLANK_SIZE);
	string_view seq=seq_itr->second;
	  if(junc_seq.find(seq)!=junc_seq.end()){
	  to_erase.push_back(seq); //remove edges that differ by only one.
	  }
	  //add exon position information; 
	  stringstream ss(id);
	  string gene,chrom,pos_str,extra_info;
	  getline(ss,gene, EXON_ID_DELIM);
	  getline(ss,chrom, EXON_ID_DELIM);
	  getline(ss,pos_str, EXON_ID_DELIM);
	  getline(ss,extra_info, EXON_ID_DELIM);

	  if(ref_pos[gene].add_boundary(chrom,atoi(pos_str.c_str()),type,extra_info.at(extra_info.size()-2)))
	    junc_seq[seq]=id; //add junction sequence to database
      }
    }
    if(junc_seq.size()==0){
      cerr << "Found no compatible sequences in exon reference fasta file"<< endl;
      exit(1);
    }

    sort(to_erase.begin(),to_erase.end());
    to_erase.erase(unique(to_erase.begin(),to_erase.end()),to_erase.end());
    //now loop again and remove all the black listed junctions
    cerr << "Removing " << to_erase.size() 
	 << " exon edge sequences for being non-unique" << endl;
    for(int i=0; i<to_erase.size(); i++)
      junc_seq.erase(to_erase.at(i));
    cerr << junc_seq.size() << " edges remaining" << endl;

    if(type==START_LABEL)
      junc_seq_start=junc_seq;
    else
      junc_seq_end=junc_seq;
  };

 public:

  inline unordered_map<string_view,string>::iterator find(string_view & key, string type){
    if(type==START_LABEL)
      return junc_seq_start.find(key);
    else
      return junc_seq_end.find(key);
  };
  inline unordered_map<string_view,string>::iterator end(string type){
    if(type==START_LABEL)
      return junc_seq_start.end();
    else
      return junc_seq_end.end();
  };

  void read_fasta( unordered_map<string,string> & seqs){
    read_fasta(seqs,START_LABEL);
    read_fasta(seqs,END_LABEL);
  }

  void get_gene_coords(string gene, string & chrom, int & start, int &end){
    ref_pos[gene].get_position(chrom, start, end);
  }

  string get_closest_junction(string gene, string chrom, int r1_end, int r2_start){
    int jend,jstart;
    if(ref_pos[gene].get_closest_junction(chrom,r1_end,r2_start,jend,jstart)){
      stringstream junc_ss;
      string strand="(+)";
      if(jstart<jend) strand="(-)";
      junc_ss << gene << EXON_ID_DELIM << chrom << EXON_ID_DELIM << jend << END_LABEL << strand ;
      junc_ss << EXON_ID_DELIM ;
      junc_ss << gene << EXON_ID_DELIM << chrom << EXON_ID_DELIM << jstart << START_LABEL << strand ;
      return junc_ss.str();
    }
    else
      return "";
  };
  //loop over the junctions and read the IDS:
    /**unordered_map<string,string >::iterator ji = ref_seq->begin();
    bool found=false;
    for(;ji!=ref_seq->end(); ji++){
      stringstream ss(ji->first);
      string field;
      getline(ss,field, EXON_ID_DELIM);
      if(field!=gene) continue;
      getline(ss,chrom, EXON_ID_DELIM);
      getline(ss,field, EXON_ID_DELIM);
      int pos=atoi(field.c_str());
      if(!found){
	found=true;
	start=pos;
	end=pos;
      }
      if(pos<start) start=pos;
      if(pos>end) end=pos;
    } //in case the format looks wrong.
    if(!found){
      cerr << "Could not find gene: " << gene << " in reference." << endl;
      exit(1);
    }
    };**/

  /**pair<int,int> get_junction_from_pair(string gene, int start, int end){
     };**/

};


