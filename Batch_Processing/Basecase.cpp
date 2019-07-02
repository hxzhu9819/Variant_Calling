#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <bitset>
using namespace std;

int wordSize = 288;
int BatchSize = 100*8*1024;

class Words{
	public:
		string header;
		vector<string> wordr;
		vector<string> wordh;
		int curr_r;
		int curr_h;
};
vector<Words> words;


class Read{
	public:
		int l;
		string bp;
		vector<int> q;
		Read(int,string);
};


class Haplotype{
	public:
		int l;
		string bp;
		int log2Init;
		int init;
		Haplotype(int,string);
};


class Region{
	public:
		Words words;
		int numReads;
		int numHap;
		vector<Read> reads;
		vector<Haplotype> haplotypes;
		Region(int);
		void generate_read_bin();
		void generate_hap_bin();
};

Region::Region(int n){
	numReads = n;
}

int getDigits(string& str){
	stringstream ss(str);
	string trash;
	int digit;
	ss>>trash>>digit;
	return digit;	
}

Read::Read(int length, string basePair){
	l = length;
	bp = basePair;
	q.reserve(l);
}

Haplotype::Haplotype(int length, string basePair){
	l = length;
	bp = basePair;
}

Read get_read(ifstream& infile, string& txt, int l){
	string bp;
	string line = txt;
	auto index = line.find(' ');
	bp = line.substr(++index);
//	cout<<bp<<endl;
	Read read(l,bp);
	
	//get quality
	getline(infile,line);
	for(int i = 0; i < l; i++){
		getline(infile,line);
		read.q.push_back(stoi(line) > 64 ? 64 : stoi(line));
	} //for
	
//	for(auto a : read.q){
//		cout<<a<<" ";
//	}
//	cout<<endl;
	return read;
}

void decToBinary(int n){ 
	int binaryNum[32]; 
	int i = 0; 
	while (n > 0) { 
		binaryNum[i] = n % 2; 
		n = n / 2; 
		i++; 
	} 
	for (int j = i - 1; j >= 0; j--) 
		cout << binaryNum[j]; 
} 

void Region::generate_read_bin(){
	for(auto read: reads){
		string out = "01";
		out.append(bitset<14>(read.l).to_string());
		for(int i = 0; i < read.l;i++){
			string base;
			if(read.bp[i] == 'A'){
				base = "00";
			}
			else if(read.bp[i] == 'T'){
				base = "01";
			}
			else if(read.bp[i] == 'G'){
				base = "10";
			}
			else if(read.bp[i] == 'C'){
				base = "11";
			}
			else{
				base ="ER";
			}
			string quality = bitset<6>(read.q[i]).to_string();
			base.append(quality);
			out.append(base);
		}//for per read
		words.wordr.push_back(out);
		cout<<out<<endl;
		cout<<"-------"<<endl;
	}//for per region
}

void Region::generate_hap_bin(){
	for(auto hap : haplotypes){
		string out = "10";
		out.append(bitset<14>(hap.l).to_string());
		string log2Value = bitset<32>(hap.log2Init).to_string();
		string initValue = bitset<32>(hap.init).to_string();
		out.append(log2Value);
		out.append(initValue);
		for(int i = 0; i < hap.l;i++){
			string base;
			if(hap.bp[i] == 'A'){
				base = "00";
			}
			else if(hap.bp[i] == 'T'){
				base = "01";
			}
			else if(hap.bp[i] == 'G'){
				base = "10";
			}
			else if(hap.bp[i] == 'C'){
				base = "11";
			}
			else{
				base ="ER";
			}
			out.append(base);
		}//for per hap
		cout<<out<<endl;
		words.wordh.push_back(out);
		cout<<"-------"<<endl;
	}//for per region
}

int hap_size(Region& region){
	int total = 0;
	for(auto hap:region.haplotypes){
		total += 80;
		total += hap.l * 2;
	}
	return total;
}

int read_size(Region& region){
	int total = 0;
	for(auto read:region.reads){
		total += 16;
		total += read.l * 8;
	}
	return total;
}


int main(int argc, char *argv[]) {
	ifstream infile("test");
	string line;
	if(infile.is_open()){
		while(getline(infile,line)){ if(line[0] == '#') break;}
		int numReads = getDigits(line);		
		Region region(numReads);
			
		//import reads
		while(getline(infile,line)){
			if(line.find('l') != 0) break;
			int l = getDigits(line);
			getline(infile,line);
			// if bp contains N, then skip and look for the next read
			if(line.find('N')){
				for(int i = 0; i < l + 1; i++){
					getline(infile,line);
				} 
				continue;
			}
			region.reads.push_back(get_read(infile, line, l));
		}
		region.generate_read_bin();
		
		//import haplotype
		if(line.find("#h") == 0){
			region.numHap = getDigits(line);
			for(int i = 0; i < region.numHap; i++){
				int l,trash,log2init,init;
				string bp;
				getline(infile,line);
				auto index = line.find(' ');
				bp = line.substr(++index);
				getline(infile,line);
				l = getDigits(line);
				getline(infile,line);
				
				//TODO:
				init = getDigits(line);
				init = 7;
				log2init = 378358;
				//END TODO
				Haplotype haplotype(l,bp);
				haplotype.init = init;
				haplotype.log2Init = log2init;
				region.haplotypes.push_back(haplotype);	
			}//for
			region.generate_hap_bin();
		}//if
		
		//generate header
		string header = "0100000000000000";
		header.append(bitset<16>(region.numReads).to_string());
		header.append(bitset<16>(region.numHap).to_string());
		
		//Generate header and Batch data
		int hpStartAddr;
		int addrLastRow;
		if(read_size(region) + hap_size(region) + wordSize < BatchSize){
			cout<<read_size(region)<<"good"<<endl;
			hpStartAddr = 1 + region.numReads;
			addrLastRow = 1 + region.numReads + region.numHap;
		}
		
	} //if is_open
	
}


