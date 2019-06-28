import java.math.*;


class Read {
	public int l;
	public String bp;
	public int[] quality;
	public long bin;
	public string
	
	/**
	  * Convert a read to bin.
	  * bp(3 bits) + Quality(6 bits) per base pair 
	  * bin: bp1(9 bits) bp2(9 bits) bp3 bp4 .....
	  * Modifies Read.bin
	*/
	public void generate_read_bin(){
		bin = (0b01 << 14) + l;
		for(int i = 0; i < l; i++){
			// bp+Q
			int base;
			if(bp.charAt(i) == 'A'){
				base = 0b000;
			}
			else if(bp.charAt(i) == 'T'){
				base = 0b001;
			}
			else if(bp.charAt(i) == 'G'){
				base = 0b010;
			}
			else if(bp.charAt(i) == 'C'){
				base = 0b011;
			}
			else{
				base =0b100;
			}
			base = base << 6;
			base += (quality[i] > 64 ? 64 : quality[i]);
			//add to bin
			bin = (bin << 9) + base;
			System.out.print(Integer.toBinaryString(base) + '\n');
		}//for
		System.out.print(Long.toBinaryString(bin) + '\n');
		System.out.print(bin);
	}
}


class Untitled {
	public static void main(String[] args) {
		Read r1 = new Read();
		r1.l = 93;
		r1.bp = "NNNNNNNNNNAGTCATNNNNNNNNNNNNNNNNNNGGTATCTTTACTATAAAAGCTATTGTGTAAGCTAGTCATANTNNNNNGTTGGCTCAGGA";
		r1.quality = new int[]{2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 18, 17, 32, 34, 32, 36, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 12, 12, 12, 20, 12, 33, 37, 39, 31, 36, 39, 39, 40, 38, 40, 38, 39, 37, 39, 39, 39, 37, 35, 26, 26, 30, 31, 31, 34, 35, 36, 34, 36, 32, 35, 34, 34, 29, 34, 35, 2, 11, 2, 2, 2, 2, 2, 11, 11, 17, 30, 30, 34, 31, 32, 25, 30, 34, 19};
		/////
		
		r1.generate_read_bin();
		

	
	}
}