# Variant_Calling Hardware Version

This is the hardware version of variant_calling. It applies multi-thread and have intergrated data-tranfer interface and binary converter.



## Added files

* `JNIinterface.java`:

  Provides the JNI interface class. Notice that all the member functions are native.

  * `JNI_send_string_to_c`: take a Java String as an argument and send it to a C function
  * `JNI_request_from_c`: Ask the corresponding C function for a String

* `slave.c`:

  C implementations of the native functions decleared in `JNIinterface.java`

  * `Java_JNIinterface_JNI_1send_1string_1to_1c`: C implementation for `JNI_send_string_to_c` in `JNIinterface.java`. It receives the String, converts it to `char*` and call `dma_send`.
  * `dma_send`:  This function fills a buffer with input data and then uses DMA to copy that buffer into each of the 4 DDR DIMMS.
  * `Java_JNIinterface_JNI_1request_1from_1c`: C implementation for `JNI_request_from_c` in `JNIinterface.java`. It calls `dma_read` to get `char*` and send it back to Java.
  * `dma_read`: This function fetch data from DDR to buffer, and store it to a global variable.

* `jni.h`

  JNI requires this file to be moved to the working directory. 

* `jni_md.h`

  JNI requires this file to be moved to the working directory. 

## Modified files

* `PairHMMLikelihoodCalculationEngine.java`:
  * `toBinary`: helper function for `hardware_process`. It takes values and desired number of digits as input, and output the value in binary in string  with required length.
  * `hardware_process`: See JavaDoc for detail. This function takes all the required data for PairHMM as input, convert them into binary according to the assigned bitmap, and call `JNI_send_string_to_c` to send the binary to FPGA.

## TODO

* Currently, the data-receiving part of `hardware_process` has not been integrated yet since we currently have not decided the output format of data from FPGA.
* The code has not been tested with FPGA yet. Need to figure out how we can write stuff to the buffer in Host

## Version Info

Commit: [a5b7215](https://github.com/hxzhu9819/Variant_Calling/commit/a5b72156379eae9ec8967bd98092c8e1f225ad8b): Integrate JNI interface. The code has removed exact

Commit: [dc79b08](https://github.com/hxzhu9819/Variant_Calling/commit/dc79b085c64616616dbe874c1996a2cc296c2e01): Add converter to the program
