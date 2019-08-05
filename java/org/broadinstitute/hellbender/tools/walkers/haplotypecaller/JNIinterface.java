class JNIinterface {
	// Load Library
	static{
		System.loadLibrary("JNIinterface");
	}

	//Functions that will apply c functions

	/* Send string to c
	   @parameters: String a
	   @output: void
	*/
	public native void JNI_send_string_to_c(String str);

	/* Receive string from c
		   @parameters: void
		   @output: a String object from c
	*/
	public native String JNI_request_from_c();

	/* Test Main Entry
	*/
	//public static void main(String[] args) {
		//JNIinterface test = new JNIinterface();
		//test.JNI_send_string_to_c("This is from JAVA");

		//String rec_str = test.JNI_request_from_c();
		//System.out.println(rec_str);
	//}
}
