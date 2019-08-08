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

||||||| merged common ancestors
# Variant_Calling
This is the software implementation after PairHMM 
=======
# GATK prune method
This is the prune method software implementation of after PairHMM. The prune method is based on [GATK 4.0.11](https://github.com/broadinstitute/gatk/tags).

## Changes after PairHMM
The prune method applies bound values for genotype likelihoods calculations, and recompute the exact likelihoods when bound check fails.

### Filter part
The main changes happen in the functions `filterPoorlyModeledReads` and `readIsPoorlyModelled` in the file `/java/org/broadinstitute/hellbender/utils/genotyper/ReadLikelihoods.java`.

When bound check fails use `recomputeOneRead` function to calculate the exact likelihoods values for failed read against all haplotypes.

### Assign genotype likelihoods
The recompute happens in the functions `marginalLikelihoods` in the file `/java/org/broadinstitute/hellbender/utils/genotyper/ReadLikelihoods.java`.

For current version, there are two methods that decide which reads to recompute.
* recompute the reads with worst bound gap
* recompute the reads in the region one by one

The second version has an additional input parameter `readPT` to determine the read point.

When recompute, the algorithm mainly do the following steps:
* determine the best haplotype index according to the allele-haplotype mapping table.

* determine the haplotype index which will be chosen as cap value in the normalize step.

* recompute the likelihoods values for the read against the haplotypes chosen.

* re-check the best haplotype index for each allele, and re-check the cap haplotype index for the read. If either one is changed, go back to step 1 and do the recompute.

* normalize (cap the likelihoods) after all index keep the same.

### Multi-threads for hardware version
To combine with hardware, multi-threads is used. For now the main progress is divided into three threads

* Handle sanity check, assemble the haplotypes

* Handle the PairHMM calculation

* Assign the genotype likelihood and pass-emit

## Important functions added/modified
* **`/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/PairHMMLikelihoodCalculationEngine.java`**

  * `computeReadLikelihoods`: divided into three sub functions. The first two generate the gap penalties and processed read lists, and the third one do the compute likelihoods.

  * `getPenaltyMap`:  generate the gap penalties. It will be kept in the `ReadLikelihoods` class for recomputation.

* **`/java/org/broadinstitute/hellbender/utils/genotyper/ReadLikelihoods.java`**

  The following variables are added:
  * `double cap_difference`: used for normalize likelihoods

  * `PairHMM recompute_tool`: for use of recompute.

  * `List<Map<GATKRead, byte[]>> gapPenalties`: gap penalties got from likelihoods calculation.

  * `List<List<GATKRead>> processedReadsList`: processed reads used for compute likelihoods

  The following functions are added or modified:

  * Constructor `ReadLikelihoods`: in this Constructor, add the initialization of the PairHMM tool `recompute_tool` with the PairHMM tool in the likelihoodCalculationEngine.

  * `getBestIndex`: the function is used to find the best alternative haplotype index. Only used in recompute part of the function `marginalLikelihoods`

  * `marginalLikelihoods`: two versions of this function, one for reads with worst bound gap, and one for recompute reads one by one. Recomputation happens here.

  * `filterPoorlyModeledReads`: add a parameter `cap_difference`. After determine the removed indices, apply `removeSampleProcessedReads` to remove the reads in the `processedReadsList`

  * `readIsPoorlyModelled`: recompute one read against all haplotypes if bound check fails. Then normalize (cap) this read.

  * `recomputeOneRead`: recompute One read against all haplotypes.

  * `recomputeOneHap`: recompute one read against one haplotype.

  * `removeSampleProcessedReads`: given the removed read indices, remove these reads in the `processedReadsList`.

  * `get_statics`: print out the workload for original method and recompute steps.

* **`/java/org/broadinstitute/hellbender/utils/pairhmm/PairHMM.java`**

  * `computeOneReadLikelihoodGivenHaplotypeLog10`: only calculate the likelihood for one Read against one Haplotype

* **`/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCaller.java`**
  * `applyAllRegion`: the function that handle three threads that complete a progress for an active region

  * `ThreadOne` to `ThreadFour`: different threads that handle different steps of GATK algorithm

* **`/java/org/broadinstitute/hellbender/tools/walkers/haplotypecaller/HaplotypeCallerEngine.java`**
  * Line 148 to line 169: lists are containers for outputs of different threads, and key objects are for locking threads when there are conflicts.

  * `callRegionStepOne` to `callRegionStepFour`: divide the `callRegion` function into four minor functions, and each one is called by different threads.

## TODO

When the hardware interface is ready, test the code. There may still exist bugs or conflicts when applying the multi-threads version.

## Version Description
The hardware branch originates from tag v1.0 (commit `4e085f3`)

* commit `a420463`: generate three threads for callRegion. hardware interface is not injected yet. exact values are still kept but not used.

* commit `70fc18c`: hardware interface is injected. The exact value calculations are removed, except in PairHMM step. Since PairHMM is due to hardware, it should not be used anyway.
>>>>>>> d6f8cf29e45a276c48ecbf54a145323da48d2ad1
