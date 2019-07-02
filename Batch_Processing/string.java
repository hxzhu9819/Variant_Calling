import java.math.*;
import java.util.*;
import javax.swing.text.*;

class Region{
	int numRead;
	int numHaplotype;
	List<Read> reads = new ArrayList<Read>();
	List<Hap> haps= new ArrayList<Hap>();
	List<String> readsbin= new ArrayList<String>();
	List<Integer> numliner= new ArrayList<Integer>();
	List<String> hapsbin= new ArrayList<String>();
	List<Integer> numlineh= new ArrayList<Integer>();
	String bin;
	
	/**
	 * @input: num: value to be converted to binary string.
	 * @input: digits: the number of digits 
     * @return: binary string
	 */
	public static String toBinary(int num, int digits) {
		char[] chars = new char[digits];
		Arrays.fill(chars, '0');
		String cover = new String(chars);
		String s = Integer.toBinaryString(num);
		return s.length() < digits ? cover.substring(s.length()) + s : s;
	}

	/**
	 * Generate binary string for all reads and haplotypes in List reads and List haps
	 * Call parse() to parse based on word size
	 * Generate Header and arrange the words
	*/
	public void generate_bin(){
		for(int i = 0; i < reads.size(); i++){
			reads.get(i).generate_read_bin();
			parse(reads.get(i).bin,0);
		}
		for(int i = 0; i < haps.size(); i++){
			haps.get(i).generate_hap_bin();
			parse(haps.get(i).bin,1);
		}
		//check size TODO
		int batchSize = 20;
		int sramSize = 1000;
		int outputSize = numRead * numHaplotype;
		int totNumLine = numliner.stream().mapToInt(Integer::intValue).sum() + numlineh.stream().mapToInt(Integer::intValue).sum();
		if(totNumLine< sramSize - outputSize){
			System.out.println("Sufficient Space!");
			System.out.println(numliner.stream().mapToInt(Integer::intValue).sum()  + " totNumLiner <--> totNumLineh " + numlineh.stream().mapToInt(Integer::intValue).sum());
			//Generate Header
			bin = "0000000000000000" + toBinary(numRead,16) + toBinary(numHaplotype,16) + toBinary(numliner.stream().mapToInt(Integer::intValue).sum(),16) + toBinary(numliner.stream().mapToInt(Integer::intValue).sum(),288 - 16 * 4) ;
			//System.out.println("-----header-----\n" + bin + "\n");
			bin = bin + "\n";
			for(String s : readsbin){
				bin += s + "\n";
			}
			for(String s : hapsbin){
				bin += s + "\n";
			}
			System.out.println(bin);
		}
		else{
			//TODO
			System.out.println("Insufficient Size");
		}
		
	}
	/**
	 * Parse binary based on word size
	 * Store the word in readsbin, and hapsbin accordingly
	 * Record the number of lines taken by each read and haplotype
	*/
	public void parse(String bin, int type){
		int i = 0;
		int numLine = 0;
		while(i < bin.length()){
			if(type == 0){
				readsbin.add( bin.substring(i,i+288 > bin.length() ? bin.length() : i+288) );
//				System.out.print("read!");
			}
			else{
				hapsbin.add( bin.substring(i,i+288 > bin.length() ? bin.length() : i+288) );
//				System.out.print("hap!");
			}
//			System.out.print( bin.substring(i,i+288 > bin.length() ? bin.length() : i+288) + "\n-----\n" );
			i += 288;
			numLine++;
		}
		
		if(type == 0){
			numliner.add(numLine);
		}
		else{
			numlineh.add(numLine);
		}
//		System.out.println("Line num:" + numLine);
	}
	
}

class Read {
	/**
	 * @input: num: value to be converted to binary string.
	 * @input: digits: the number of digits 
     * @return: binary string
	*/
	public static String toBinary(int num, int digits) {
			char[] chars = new char[digits];
			Arrays.fill(chars, '0');
			String cover = new String(chars);
			String s = Integer.toBinaryString(num);
			return s.length() < digits ? cover.substring(s.length()) + s : s;
	}
	
	public int l;
	public String bp;
	public int[] quality;
	public String bin;
	public String isN = " ";
	
	/**
	  * Convert a read to bin.
	  * bp(3 bits) + Quality(6 bits) per base pair 
	  * bin: bp1(9 bits) bp2(9 bits) bp3 bp4 .....
	  * Modifies Read.bin
	*/
	public void generate_read_bin(){
		bin = "01" + toBinary(l,14);
		for(int i = 0; i < l; i++){
			// bp+Q
			String base;
			if(bp.charAt(i) == 'A'){
				base = "00";
				isN = isN + "0";
			}
			else if(bp.charAt(i) == 'T'){
				base = "01";
				isN = isN + "0";
			}
			else if(bp.charAt(i) == 'G'){
				base = "10";
				isN = isN + "0";
			}
			else if(bp.charAt(i) == 'C'){
				base = "11";
				isN = isN + "0";
			}
			else{
				base ="00";
				isN = isN + "1";
			}
			base = base + toBinary(quality[i] > 64 ? 64 : quality[i], 6);
			//add to bin
			bin = bin + base;
			//System.out.print(base + '\n');
		}//for
		bin = bin.concat(isN.trim());
//		System.out.print(bin + "\n-----\n");
		//System.out.print(isN);
		
//		//拆成288
//		int i = 0;
//		int numLine = 0;
//		while(i < bin.length()){
//			System.out.print( bin.substring(i,i+288 > bin.length() ? bin.length() : i+288) + "\n-----\n" );
//			i += 288;
//			numLine++;
//		}
//		System.out.println("Line num:" + numLine);
		
	}
}

class Hap {
	public static String toBinary(int num, int digits) {
		char[] chars = new char[digits];
		Arrays.fill(chars, '0');
		String cover = new String(chars);
		String s = Integer.toBinaryString(num);
		return s.length() < digits ? cover.substring(s.length()) + s : s;
	}
	
	public int l;
	public String bp;
	public int init;
	public int loginit;
	public String bin;
	
	/**
	  * Convert a hap to bin.
	  * bp(3 bits) + Quality(6 bits) per base pair 
	  * bin: bp1(9 bits) bp2(9 bits) bp3 bp4 .....
	  * Modifies Read.bin
	*/
	public void generate_hap_bin(){
		//overhead
		bin = "10" + toBinary(l,14) + toBinary(loginit,32) + toBinary(init,32);
		for(int i = 0; i < l; i++){
			// bp+Q
			String base;
			if(bp.charAt(i) == 'A'){
				base = "00";
			}
			else if(bp.charAt(i) == 'T'){
				base = "01";
			}
			else if(bp.charAt(i) == 'G'){
				base = "10";
			}
			else if(bp.charAt(i) == 'C'){
				base = "11";
			}
			else{
				base ="00";
			}
			//add to bin
			bin = bin + base;
			//System.out.print(base + '\n');
		}//for
		
//		//拆成288
//		int i = 0;
//		int numLine = 0;
//		while(i < bin.length()){
//			System.out.print( bin.substring(i,i+288 > bin.length() ? bin.length() : i+288) + "\n-----\n" );
//			i += 288;
//			numLine++;
//		}
//		System.out.print("Line num:" + numLine);
		
	}
}



class Main {
	public static void main(String[] args) {
		Read r1 = new Read();
		r1.l = 93;
		r1.bp = "NNNNNNNNNNAGTCATNNNNNNNNNNNNNNNNNNGGTATCTTTACTATAAAAGCTATTGTGTAAGCTAGTCATANTNNNNNGTTGGCTCAGGA";
		r1.quality = new int[]{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 18, 17, 32, 34, 32, 36, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 12, 12, 20, 12, 33, 37, 39, 31, 36, 39, 39, 40, 38, 40, 38, 39, 37, 39, 39, 39, 37, 35, 26, 26, 30, 31, 31, 34, 35, 36, 34, 36, 32, 35, 34, 34, 29, 34, 35, 2, 11, 2, 2, 2, 2, 2, 11, 11, 17, 30, 30, 34, 31, 32, 25, 30, 34, 19};
		
		Read r2 = new Read();
		r2.l = 92;
		r2.bp = "NNNNNNNNNAGTCATNNNNNNNNNNNNNNNNNNGGTATCTTTACTATAAAAGCTATTGTGTAAGCTAGTCATANTNNNNNGTTGGCTCAGGA";
		r2.quality = new int[]{2, 2, 2, 2, 2, 2, 2, 2, 2, 18, 17, 32, 34, 32, 36, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 12, 12, 20, 12, 33, 37, 39, 31, 36, 39, 39, 40, 38, 40, 38, 39, 37, 39, 39, 39, 37, 35, 26, 26, 30, 31, 31, 34, 35, 36, 34, 36, 32, 35, 34, 34, 29, 34, 35, 2, 11, 2, 2, 2, 2, 2, 11, 11, 17, 30, 30, 34, 31, 32, 25, 30, 34, 19};

		Hap h1 = new Hap();
		h1.l = 86;
		h1.bp = "GGCTTTAGGGAGCCATAGGTGGAGACCGTAAAGAGGTATCTTTACTATAAAAGCTATTGTGTAAGCTAGTCATATTAAGTTATTGGCTCAGGAGTTTGATAGTTCTTGGGCAGTAAGAGTGAGTAATAGAATATTCAGTGAGCCTAGGGTGTTGTGAGTGTAAATTAGTGCGATGAGTAGGGGAAGGGAGCCTACTAGGGTGTAGAATAGGAAGTATGTGCCTGCGTTCAGGCGTTCTGG";
		h1.init = 7;
		h1.loginit = 378358;
		
		Region reg1 = new Region();
		reg1.reads.add(r1);
		reg1.reads.add(r2);
		reg1.haps.add(h1);
		reg1.numRead = 2;
		reg1.numHaplotype = 1;
//	r1.generate_read_bin();
//	h1.generate_hap_bin();
	
		reg1.generate_bin();
		
		
	}
}