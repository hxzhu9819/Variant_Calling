package org.broadinstitute.hellbender.utils.genotyper;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.util.Locatable;
import htsjdk.variant.variantcontext.Allele;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.objects.Object2IntMap;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;
import org.apache.commons.collections.ListUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.downsampling.AlleleBiasedDownsamplingUtils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
// added by Chenhao
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


/**
 * Read-likelihoods container implementation based on integer indexed arrays.
 *
 * @param <A> the type of the allele the likelihood makes reference to.
 *
 * Note: this class uses FastUtil collections for speed.
 *
 * @author Valentin Ruano-Rubio &lt;valentin@broadinstitute.org&gt;
 */
public class ReadLikelihoods<A extends Allele> implements SampleList, AlleleList<A> {

    /**
     * Index indicaintg that the reference allele is missing.
     */
    private static final int MISSING_REF = -1;

    /**
     * Reads by sample index. Each sub array contains reference to the reads of the ith sample.
     */
    protected final GATKRead[][] readsBySampleIndex;

    /**
     * Indexed per sample, allele and finally read (within sample).
     * <p>
     *     valuesBySampleIndex[s][a][r] == lnLk(R_r | A_a) where R_r comes from Sample s.
     * </p>
     */
    //protected final double[][][] valuesBySampleIndex;
    public final double[][][] valuesBySampleIndex;
    //Xiao added this: keeps record of which valuesBySampleIndex has been recomputed
    public boolean[][][] recomputeRecord;
    //Xiao added this: keeps record of old allele index chosen in marginalization
    private int [][][] bestHapMapLO;
    private int [][][] bestHapMapUP;

    /**
     * Sample list
     */
    protected final SampleList samples;

    /**
     * Allele list
     */
    protected AlleleList<A> alleles;

    /**
     * Cached allele list.
     */
    private List<A> alleleList;

    /**
     * Cached sample list.
     */
    private List<String> sampleList;

    //Xiao added this
    private List<Haplotype> haplotypes;

    /**
     * Maps from each read to its index within the sample.
     *
     * <p>In order to save CPU time the indices contained in this array (not the array itself) is
     * lazily initialized by invoking {@link #readIndexBySampleIndex(int)}.</p>
     */
    private final Object2IntMap<GATKRead>[] readIndexBySampleIndex;

    /**
     * Index of the reference allele if any, otherwise {@link #MISSING_REF}.
     */
    private int referenceAlleleIndex = MISSING_REF;

    /**
     * Caches the read-list per sample list returned by {@link #sampleReads}
     */
    private final List<GATKRead>[] readListBySampleIndex;

    /**
     * Sample matrices lazily initialized (the elements not the array) by invoking {@link #sampleMatrix(int)}.
     */
    private final LikelihoodMatrix<A>[] sampleMatrices;

    //added by Chenhao: create a PairHMM tool for recompute
    private PairHMM recompute_tool = null;

    /**
     * Is this container expected to have the per-allele liklihoods calculations filled in.
     */
    public boolean hasFilledLikelihoods() {
        return true;
    }

    /**
     * Constructs a new read-likelihood collection.
     *
     * <p>
     *     The initial likelihoods for all allele-read combinations are
     *     0.
     * </p>
     *
     * @param samples all supported samples in the collection.
     * @param alleles all supported alleles in the collection.
     * @param reads reads stratified per sample.
     *
     * @throws IllegalArgumentException if any of {@code allele}, {@code samples}
     * or {@code reads} is {@code null},
     *  or if they contain null values.
     */
    // added by Chenhao: new constructor
    @SuppressWarnings({"rawtypes", "unchecked"})
    public ReadLikelihoods(final SampleList samples,
                           final AlleleList<A> alleles,
                           final Map<String, List<GATKRead>> reads,
                           final PairHMM pairHMM) {
        Utils.nonNull(alleles, "allele list cannot be null");
        Utils.nonNull(samples, "sample list cannot be null");
        Utils.nonNull(reads, "read map cannot be null");

        // added by Chenhao: initialize the pairHMM tool
        this.recompute_tool = pairHMM;

        this.samples = samples;
        this.alleles = alleles;

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        readsBySampleIndex = new GATKRead[sampleCount][];
        readListBySampleIndex = (List<GATKRead>[])new List[sampleCount];
        valuesBySampleIndex = new double[sampleCount][][];
        //System.err.print("ReadLikelihoods: ReadLikelihoods0\n");
        recomputeRecord = new boolean[sampleCount][][];
        bestHapMapLO = new int[sampleCount][][];
        bestHapMapUP = new int[sampleCount][][];
        referenceAlleleIndex = findReferenceAllele(alleles);

        readIndexBySampleIndex = new Object2IntMap[sampleCount];

        setupIndexes(reads, sampleCount, alleleCount);

        sampleMatrices = (LikelihoodMatrix<A>[]) new LikelihoodMatrix[sampleCount];
    }

    @SuppressWarnings({"rawtypes", "unchecked"})
    public ReadLikelihoods(final SampleList samples,
                           final AlleleList<A> alleles,
                           final Map<String, List<GATKRead>> reads) {
        Utils.nonNull(alleles, "allele list cannot be null");
        Utils.nonNull(samples, "sample list cannot be null");
        Utils.nonNull(reads, "read map cannot be null");

        this.samples = samples;
        this.alleles = alleles;

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        readsBySampleIndex = new GATKRead[sampleCount][];
        readListBySampleIndex = (List<GATKRead>[])new List[sampleCount];
        valuesBySampleIndex = new double[sampleCount][][];
        //System.err.print("ReadLikelihoods: ReadLikelihoods0\n");
        recomputeRecord = new boolean[sampleCount][][];
        bestHapMapLO = new int[sampleCount][][];
        bestHapMapUP = new int[sampleCount][][];
        referenceAlleleIndex = findReferenceAllele(alleles);

        readIndexBySampleIndex = new Object2IntMap[sampleCount];

        setupIndexes(reads, sampleCount, alleleCount);

        sampleMatrices = (LikelihoodMatrix<A>[]) new LikelihoodMatrix[sampleCount];
    }

    // Internally used constructor.
    @SuppressWarnings({"unchecked", "rawtypes"})
    ReadLikelihoods(final AlleleList alleles,
                    final SampleList samples,
                    final GATKRead[][] readsBySampleIndex,
                    final Object2IntMap<GATKRead>[] readIndex,
                    final double[][][] values) {
        this.samples = samples;
        this.alleles = alleles;
        this.readsBySampleIndex = readsBySampleIndex;
        this.valuesBySampleIndex = values;

        final int sampleCount = samples.numberOfSamples();
        recomputeRecord = new boolean[sampleCount][][];
        bestHapMapLO = new int[sampleCount][][];
        bestHapMapUP = new int[sampleCount][][];
        //System.err.print("ReadLikelihoods:1\n");
        for(int s=0;s<sampleCount;s++){
            recomputeRecord[s] = new boolean[values[s].length][values[s][0].length];
            bestHapMapLO[s] = new int[values[s].length][values[s][0].length];
            bestHapMapUP[s] = new int[values[s].length][values[s][0].length];
            for(int a=0;a<values[s].length;a++){
                for(int r=0;r<values[s][a].length;r++){
                    //System.err.printf("ReadLikelihoods: ReadLikelihoods1 s=%d a=%d r=%d\n",s,a,r);
                    recomputeRecord[s][a][r]=false;
                    bestHapMapLO[s][a][r] = -1;
                    bestHapMapUP[s][a][r] = -1;
                }
            }
        }
        this.readIndexBySampleIndex = readIndex;
        this.readListBySampleIndex = (List<GATKRead>[])new List[sampleCount];

        referenceAlleleIndex = findReferenceAllele(alleles);
        sampleMatrices = (LikelihoodMatrix<A>[]) new LikelihoodMatrix[sampleCount];
    }


    // Internally used constructor.
    @SuppressWarnings({"unchecked", "rawtypes"})
    ReadLikelihoods(final List<Haplotype> haplotypes,
                    final AlleleList alleles,
                    final SampleList samples,
                    final GATKRead[][] readsBySampleIndex,
                    final Object2IntMap<GATKRead>[] readIndex,
                    final double[][][] values,
                    boolean [][][] recomputeRecord,
                    int [][][] bestHapMapLO,
                    int [][][] bestHapMapUP) {
        this.haplotypes = haplotypes;
        this.samples = samples;
        this.alleles = alleles;
        this.readsBySampleIndex = readsBySampleIndex;
        this.valuesBySampleIndex = values;
        this.recomputeRecord = recomputeRecord;
        this.bestHapMapLO = bestHapMapLO;
        this.bestHapMapUP = bestHapMapUP;
        //System.err.print("ReadLikelihoods:2\n");
        this.readIndexBySampleIndex = readIndex;
        final int sampleCount = samples.numberOfSamples();
        this.readListBySampleIndex = (List<GATKRead>[])new List[sampleCount];

        referenceAlleleIndex = findReferenceAllele(alleles);
        sampleMatrices = (LikelihoodMatrix<A>[]) new LikelihoodMatrix[sampleCount];
    }

    // Add all the indices to alleles, sample and reads in the look-up maps.
    private void setupIndexes(final Map<String, List<GATKRead>> reads, final int sampleCount, final int alleleCount) {
        for (int i = 0; i < sampleCount; i++) {
            setupSampleData(i, reads, alleleCount);
        }
    }

    // Assumes that {@link #samples} has been initialized with the sample names.
    private void setupSampleData(final int sampleIndex,
                                 final Map<String, List<GATKRead>> readsBySample,
                                 final int alleleCount) {
        final String sample = samples.getSample(sampleIndex);

        final List<GATKRead> reads = readsBySample.get(sample);
        readsBySampleIndex[sampleIndex] = reads == null
                ? new GATKRead[0]
                : reads.toArray(new GATKRead[reads.size()]);
        final int sampleReadCount = readsBySampleIndex[sampleIndex].length;

        final double[][] sampleValues = new double[alleleCount][sampleReadCount];
        valuesBySampleIndex[sampleIndex] = sampleValues;
        //System.err.print("ReadLikelihoods: setupSampleData\n");
        recomputeRecord[sampleIndex]=new boolean[alleleCount][sampleReadCount];
        bestHapMapLO[sampleIndex] = new int[alleleCount][sampleReadCount];
        bestHapMapUP[sampleIndex] = new int[alleleCount][sampleReadCount];
        for(int a=0;a<alleleCount;a++){
            for(int r=0;r<sampleReadCount;r++){
                recomputeRecord[sampleIndex][a][r]=false;
                bestHapMapLO[sampleIndex][a][r] = -1;
                bestHapMapUP[sampleIndex][a][r] = -1;
            }
        }
    }

    /**
     * Create an independent copy of this read-likelihoods collection
     */
    @VisibleForTesting
    ReadLikelihoods<A> copy() {

        final int sampleCount = samples.numberOfSamples();
        final int alleleCount = alleles.numberOfAlleles();

        final double[][][] newLikelihoodValues = new double[sampleCount][alleleCount][];

        @SuppressWarnings({"unchecked", "rawtypes"})
        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];

        for (int s = 0; s < sampleCount; s++) {
            newReadsBySampleIndex[s] = readsBySampleIndex[s].clone();
            for (int a = 0; a < alleleCount; a++) {
                newLikelihoodValues[s][a] = valuesBySampleIndex[s][a].clone();
            }
        }
        boolean [][][] newrecomputeRecord = new boolean[sampleCount][][];
        int [][][] newbestHapMapLO = new int[sampleCount][][];
        int [][][] newbestHapMapUP = new int[sampleCount][][];
        for(int s=0;s<sampleCount;s++){
            for(int a=0;a<alleleCount;a++){
                newrecomputeRecord[s][a]=recomputeRecord[s][a].clone();
                newbestHapMapLO[s][a] = bestHapMapLO[s][a].clone();
                newbestHapMapUP[s][a] = bestHapMapUP[s][a].clone();
            }
        }
        // Finally we create the new read-likelihood
        return new ReadLikelihoods<>(
                haplotypes,
                alleles,
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues,
                newrecomputeRecord,
                newbestHapMapLO,
                newbestHapMapUP);
    }


    // Search for the reference allele, if not found the index is {@link MISSING_REF}.
    private static int findReferenceAllele(final AlleleList<?> alleles) {
        return IntStream.range(0, alleles.numberOfAlleles()).filter(i -> alleles.getAllele(i).isReference()).findAny().orElse(MISSING_REF);
    }

    /**
     * Returns the index of a sample within the likelihood collection.
     *
     * @param sample the query sample.
     *
     * @throws IllegalArgumentException if {@code sample} is {@code null}.
     * @return -1 if the allele is not included, 0 or greater otherwise.
     */
    @Override
    public int indexOfSample(final String sample) {
        return samples.indexOfSample(sample);
    }

    /**
     * Number of samples included in the likelihood collection.
     * @return 0 or greater.
     */
    @Override
    public int numberOfSamples() {
        return samples.numberOfSamples();
    }

    /**
     * Returns sample name given its index.
     *
     * @param sampleIndex query index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is negative.
     *
     * @return never {@code null}.
     */
    @Override
    public String getSample(final int sampleIndex) {
        return samples.getSample(sampleIndex);
    }

    /**
     * Returns the index of an allele within the likelihood collection.
     *
     * @param allele the query allele.
     *
     * @throws IllegalArgumentException if {@code allele} is {@code null}.
     *
     * @return -1 if the allele is not included, 0 or greater otherwise.
     */
    @Override
    public int indexOfAllele(final A allele) {
        return alleles.indexOfAllele(allele);
    }

    /**
     * Returns number of alleles in the collection.
     * @return 0 or greater.
     */
    @Override
    public int numberOfAlleles() {
        return alleles.numberOfAlleles();
    }

    /**
     * Returns the allele given its index.
     *
     * @param alleleIndex the allele index.
     *
     * @throws IllegalArgumentException the allele index is {@code null}.
     *
     * @return never {@code null}.
     */
    @Override
    public A getAllele(final int alleleIndex) {
        return alleles.getAllele(alleleIndex);
    }

    /**
     * Returns the reads that belong to a sample sorted by their index (within that sample).
     *
     * @param sampleIndex the requested sample.
     * @return never {@code null} but perhaps a zero-length array if there is no reads in sample. No element in
     *   the array will be null.
     */
    public List<GATKRead> sampleReads(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        final List<GATKRead> extantList = readListBySampleIndex[sampleIndex];
        if (extantList == null) {
            return readListBySampleIndex[sampleIndex] = Collections.unmodifiableList(Arrays.asList(readsBySampleIndex[sampleIndex]));
        } else {
            return extantList;
        }
    }

    //Xiao added this
    public List<Haplotype> haplotypes() {
        return haplotypes;
    }
    /**
     * Returns a read vs allele likelihood matrix corresponding to a sample.
     *
     * @param sampleIndex target sample.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not null.
     *
     * @return never {@code null}
     */
    public LikelihoodMatrix<A> sampleMatrix(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        final LikelihoodMatrix<A> extantResult = sampleMatrices[sampleIndex];
        if (extantResult == null) {
            return sampleMatrices[sampleIndex] = new SampleMatrix(sampleIndex);
        } else {
            return extantResult;
        }
    }
    /**
     * print for debug.
     *
     */
    public void printlikelihoods(){
        System.err.print("printlikelihoods\n");
        for(int i=0;i<numberOfSamples();i++){
            for(int a=0;a<alleles.numberOfAlleles();a++){
                for(int r=0;r<readsBySampleIndex[i].length;r++){
                    System.err.printf("a=%d r=%d values=%f\n",a,r,valuesBySampleIndex[i][a][r]);
                }
            }
        }

    }

    /**
     * Adjusts likelihoods so that for each read, the best allele likelihood is 0 and caps the minimum likelihood
     * of any allele for each read based on the maximum alternative allele likelihood.
     *
     * @param maximumLikelihoodDifferenceCap maximum difference between the best alternative allele likelihood
     *                                           and any other likelihood.
     *
     * @throws IllegalArgumentException if {@code maximumDifferenceWithBestAlternative} is not 0 or less.
     */
    public void normalizeLikelihoods(final double maximumLikelihoodDifferenceCap) {
        Utils.validateArg(maximumLikelihoodDifferenceCap < 0.0 && !Double.isNaN(maximumLikelihoodDifferenceCap),
                "the minimum reference likelihood fall must be negative");

        if (maximumLikelihoodDifferenceCap == Double.NEGATIVE_INFINITY) {
            return;
        }

        final int alleleCount = alleles.numberOfAlleles();
        if (alleleCount == 0){ // trivial case there is no alleles.
            return;
        } else if (alleleCount == 1) {
            return;
        }

        for (int s = 0; s < valuesBySampleIndex.length; s++) {
            final double[][] sampleValues = valuesBySampleIndex[s];
            final int readCount = readsBySampleIndex[s].length;
            for (int r = 0; r < readCount; r++) {
                normalizeLikelihoodsPerRead(maximumLikelihoodDifferenceCap, sampleValues, s, r);
            }
        }
    }

    // Does the normalizeLikelihoods job for each read.
    private void normalizeLikelihoodsPerRead(final double maximumBestAltLikelihoodDifference,
                                             final double[][] sampleValues, final int sampleIndex, final int readIndex) {

        final BestAllele bestAlternativeAllele = searchBestAllele(sampleIndex,readIndex,false, false);

        final double worstLikelihoodCap = bestAlternativeAllele.likelihood + maximumBestAltLikelihoodDifference;

        final int alleleCount = alleles.numberOfAlleles();

        // Guarantee to be the case by enclosing code.
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][readIndex] < worstLikelihoodCap) {
                sampleValues[a][readIndex] = worstLikelihoodCap;
            }
        }

    }

    /**
     * Returns the samples in this read-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable view on the read-likelihoods sample list.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<String> samples() {
        return Collections.unmodifiableList(sampleList == null ? sampleList = samples.asListOfSamples() : sampleList);
    }

    /**
     * Returns the samples in this read-likelihood collection.
     * <p>
     *     Samples are sorted by their index in the collection.
     * </p>
     *
     * <p>
     *     The returned list is an unmodifiable. It will not be updated if the collection
     *     allele list changes.
     * </p>
     *
     * @return never {@code null}.
     */
    public List<A> alleles() {
        return Collections.unmodifiableList(alleleList == null ? alleleList = alleles.asListOfAlleles() : alleleList);
    }


    /**
     * Search the best allele for a read.
     *
     * @param sampleIndex including sample index.
     * @param readIndex  target read index.
     *
     * @param useReferenceIfUninformative
     * @return never {@code null}, but with {@link BestAllele#allele allele} == {@code null}
     * if non-could be found.
     */
    private BestAllele searchBestAllele(final int sampleIndex, final int readIndex, final boolean canBeReference, boolean useReferenceIfUninformative) {
        final int alleleCount = alleles.numberOfAlleles();
        if (alleleCount == 0 || (alleleCount == 1 && referenceAlleleIndex == 0 && !canBeReference)) {
            return new BestAllele(sampleIndex, readIndex, -1, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY);
        }

        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        int bestAlleleIndex = canBeReference || referenceAlleleIndex != 0 ? 0 : 1;

        double bestLikelihood = sampleValues[bestAlleleIndex][readIndex];
        double secondBestLikelihood = Double.NEGATIVE_INFINITY;
        for (int a = bestAlleleIndex + 1; a < alleleCount; a++) {
            if (!canBeReference && referenceAlleleIndex == a) {
                continue;
            }
            final double candidateLikelihood = sampleValues[a][readIndex];
            if (candidateLikelihood > bestLikelihood) {
                bestAlleleIndex = a;
                secondBestLikelihood = bestLikelihood;
                bestLikelihood = candidateLikelihood;
            } else if (candidateLikelihood > secondBestLikelihood) {
                secondBestLikelihood = candidateLikelihood;
            }
        }

        // if our read is not informative against the ref we set the ref as the best allele.  This is so that bamouts don't
        // spuriously show deletions in ref reads that end in STRs
        if (useReferenceIfUninformative && canBeReference && referenceAlleleIndex != MISSING_REF && bestAlleleIndex != referenceAlleleIndex) {
            final double referenceLikelihood = sampleValues[referenceAlleleIndex][readIndex];
            if ( bestLikelihood - referenceLikelihood < BestAllele.INFORMATIVE_THRESHOLD ) {
                secondBestLikelihood = bestLikelihood;
                bestAlleleIndex = referenceAlleleIndex;
                bestLikelihood = referenceLikelihood;
            }
        }
        return new BestAllele(sampleIndex,readIndex,bestAlleleIndex,bestLikelihood,secondBestLikelihood);
    }

    public void changeReads(final Map<GATKRead, GATKRead> readRealignments) {
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final Object2IntMap<GATKRead> readIndex = readIndexBySampleIndex[s];
            final int sampleReadCount = sampleReads.length;
            for (int r = 0; r < sampleReadCount; r++) {
                final GATKRead read = sampleReads[r];
                final GATKRead replacement = readRealignments.get(read);
                if (replacement == null) {
                    continue;
                }
                sampleReads[r] = replacement;
                if (readIndex != null) {
                    readIndex.remove(read);
                    readIndex.put(replacement, r);
                }
            }
        }
    }

    /**
     * Add alleles that are missing in the read-likelihoods collection giving all reads a default
     * likelihood value.
     * @param candidateAlleles the potentially missing alleles.
     * @param defaultLikelihood the default read likelihood value for that allele.
     *
     * @return {@code true} iff the the read-likelihood collection was modified by the addition of the input alleles.
     *  So if all the alleles in the input collection were already present in the read-likelihood collection this method
     *  will return {@code false}.
     *
     * @throws IllegalArgumentException if {@code candidateAlleles} is {@code null} or there is more than
     * one missing allele that is a reference or there is one but the collection already has
     * a reference allele.
     */
    public boolean addMissingAlleles(final Collection<A> candidateAlleles, final double defaultLikelihood) {
        Utils.nonNull(candidateAlleles, "the candidateAlleles list cannot be null");
        if (candidateAlleles.isEmpty()) {
            return false;
        }
        final List<A> allelesToAdd = candidateAlleles.stream().filter(allele -> !alleles.containsAllele(allele)).collect(Collectors.toList());

        if (allelesToAdd.isEmpty()) {
            return false;
        }

        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = alleles.numberOfAlleles() + allelesToAdd.size();

        alleleList = null;
        int referenceIndex = this.referenceAlleleIndex;

        @SuppressWarnings("unchecked")
        final List<A> newAlleles = ListUtils.union(alleles.asListOfAlleles(), allelesToAdd);
        alleles = new IndexedAlleleList<>(newAlleles);

        // if we previously had no reference allele, update the reference index if a reference allele is added
        // if we previously had a reference and try to add another, throw an exception
        final OptionalInt indexOfReferenceInAllelesToAdd = IntStream.range(0, allelesToAdd.size())
                .filter(n -> allelesToAdd.get(n).isReference()).findFirst();
        if (referenceIndex != MISSING_REF) {
            Utils.validateArg(!indexOfReferenceInAllelesToAdd.isPresent(), "there can only be one reference allele");
        } else if (indexOfReferenceInAllelesToAdd.isPresent()){
            referenceAlleleIndex = oldAlleleCount + indexOfReferenceInAllelesToAdd.getAsInt();
        }

        //copy old allele likelihoods and set new allele likelihoods to the default value
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final int sampleReadCount = readsBySampleIndex[s].length;
            final double[][] newValuesBySampleIndex = Arrays.copyOf(valuesBySampleIndex[s], newAlleleCount);
            for (int a = oldAlleleCount; a < newAlleleCount; a++) {
                newValuesBySampleIndex[a] = new double[sampleReadCount];
                if (defaultLikelihood != 0.0) {
                    Arrays.fill(newValuesBySampleIndex[a], defaultLikelihood);
                }
            }
            valuesBySampleIndex[s] = newValuesBySampleIndex;
        }
        return true;
    }

    /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each read in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @return never {@code null}. The result will have the requested set of new alleles (keys in {@code newToOldAlleleMap}, and
     * the same set of samples and reads as the original.
     *
     * @throws IllegalArgumentException is {@code newToOldAlleleMap} is {@code null} or contains {@code null} values,
     *  or its values contain reference to non-existing alleles in this read-likelihood collection. Also no new allele
     *  can have zero old alleles mapping nor two new alleles can make reference to the same old allele.
     */
    @SuppressWarnings({"unchecked", "rawtypes"})
    //Original
    public <B extends Allele> ReadLikelihoods<B> marginalize(final Map<B, List<A>> newToOldAlleleMap) {
        Utils.nonNull(newToOldAlleleMap);

        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final int[] oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        // We calculate the marginal likelihoods.
        final double[][][] newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, oldToNewAlleleIndexMap, null);

        final int sampleCount = samples.numberOfSamples();

        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];

        for (int s = 0; s < sampleCount; s++) {
            newReadsBySampleIndex[s] = readsBySampleIndex[s].clone();
        }

        // Finally we create the new read-likelihood
        //System.err.print("ReadLikelihoods: marginalize original\n");
        return new ReadLikelihoods<>(
                new IndexedAlleleList(newAlleles),
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex, newLikelihoodValues);
    }
    //Prune method
    @SuppressWarnings({"unchecked", "rawtypes"})
    public <B extends Allele> List<ReadLikelihoods<B>> marginalize(final List<Haplotype> haplotypes,ReadLikelihoods<Haplotype> readLikelihoods_upperbound,ReadLikelihoods<Haplotype> readLikelihoods_exact, final Map<B, List<A>> newToOldAlleleMap,final boolean moreRecompute, int[] readPT,  int [] redoSquares ) {
        Utils.nonNull(newToOldAlleleMap);
        // 这里传进来的newToOldAlleleMap应该是Allele-haplotype的map
        // ??? 看后面的用法 newToOldAlleleMap应该是一个Allele-Allelelist的map???

        // newAlleles是mapper中的key 也就是allele的数组
        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final int[] oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        // We calculate the marginal likelihoods.
        final double[][][][] newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount,oldToNewAlleleIndexMap, null,readLikelihoods_upperbound,readLikelihoods_exact,moreRecompute, readPT, redoSquares);

        final int sampleCount = samples.numberOfSamples();

        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];
        boolean[][][] newrecomputeRecord = new boolean[sampleCount][][];
        int [][][] newbestHapMapLO = new int[sampleCount][][];
        int [][][] newbestHapMapUP = new int[sampleCount][][];
        for (int s = 0; s < sampleCount; s++) {
            newReadsBySampleIndex[s] = readsBySampleIndex[s].clone();
            newrecomputeRecord[s] = recomputeRecord[s].clone();
            newbestHapMapLO[s] = bestHapMapLO[s].clone();
            newbestHapMapUP[s] = bestHapMapUP[s].clone();

        }
        final List<ReadLikelihoods<B>> result = new ArrayList<ReadLikelihoods<B>>();
        result.add(new ReadLikelihoods<>(
                haplotypes,
                new IndexedAlleleList(newAlleles),
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues[0],
                newrecomputeRecord,
                newbestHapMapLO,
                newbestHapMapUP));
        result.add(new ReadLikelihoods<>(
                haplotypes,
                new IndexedAlleleList(newAlleles),
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues[1],
                newrecomputeRecord,
                newbestHapMapLO,
                newbestHapMapUP));
        result.add(new ReadLikelihoods<>(
                haplotypes,
                new IndexedAlleleList(newAlleles),
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues[2],
                newrecomputeRecord,
                newbestHapMapLO,
                newbestHapMapUP));



        // Finally we create the new read-likelihood
        return result;
    }


    /**
     * Perform marginalization from an allele set to another (smaller one) taking the maximum value
     * for each read in the original allele subset.
     *
     * @param newToOldAlleleMap map where the keys are the new alleles and the value list the original
     *                          alleles that correspond to the new one.
     * @return never {@code null}. The result will have the requested set of new alleles (keys in {@code newToOldAlleleMap}, and
     * the same set of samples and reads as the original.
     *
     * @param overlap if not {@code null}, only reads that overlap the location (with unclipping) will be present in
     *                        the output read-collection.
     *
     * @throws IllegalArgumentException is {@code newToOldAlleleMap} is {@code null} or contains {@code null} values,
     *  or its values contain reference to non-existing alleles in this read-likelihood collection. Also no new allele
     *  can have zero old alleles mapping nor two new alleles can make reference to the same old allele.
     */
    //Original
    public <B extends Allele> ReadLikelihoods<B> marginalize(final Map<B, List<A>> newToOldAlleleMap, final Locatable overlap) {
        Utils.nonNull(newToOldAlleleMap, "the input allele mapping cannot be null");
        if (overlap == null) {
            System.err.print("ReadLikelihoods: marginalize original 0\n");
            return marginalize(newToOldAlleleMap);
        }

        @SuppressWarnings("unchecked")
        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final int[] oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        final int[][] readsToKeep = overlappingReadIndicesBySampleIndex(overlap);
        // We calculate the marginal likelihoods.

        final double[][][] newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount, oldToNewAlleleIndexMap, readsToKeep);

        final int sampleCount = samples.numberOfSamples();

        @SuppressWarnings({"rawtypes","unchecked"})
        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = (Object2IntMap<GATKRead>[])new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];

        for (int s = 0; s < sampleCount; s++) {
            final int[] sampleReadsToKeep = readsToKeep[s];
            final GATKRead[] oldSampleReads = readsBySampleIndex[s];
            final int oldSampleReadCount = oldSampleReads.length;
            final int newSampleReadCount = sampleReadsToKeep.length;
            if (newSampleReadCount == oldSampleReadCount) {
                newReadsBySampleIndex[s] = oldSampleReads.clone();
            } else {
                newReadsBySampleIndex[s] = new GATKRead[newSampleReadCount];
                for (int i = 0; i < newSampleReadCount; i++) {
                    newReadsBySampleIndex[s][i] = oldSampleReads[sampleReadsToKeep[i]];
                }
            }
        }

        // Finally we create the new read-likelihood
        return new ReadLikelihoods<>(new IndexedAlleleList<>(newAlleles), samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex, newLikelihoodValues);
    }
    //Prune method
    public <B extends Allele> List<ReadLikelihoods<B>> marginalize(final List<Haplotype> haplotypes, final ReadLikelihoods<Haplotype> readLikelihoods_upperbound,final ReadLikelihoods<Haplotype> readLikelihoods_exact,final Map<B, List<A>> newToOldAlleleMap, final Locatable overlap,final boolean moreRecompute, int[] readPT,  int [] redoSquares) {
        Utils.nonNull(newToOldAlleleMap, "the input allele mapping cannot be null");
        if (overlap == null) {
            return marginalize(haplotypes,readLikelihoods_upperbound,readLikelihoods_exact,newToOldAlleleMap,moreRecompute, readPT, redoSquares);
        }

        @SuppressWarnings("unchecked")
        // haplotype和allele的map 互相转换
        // newToOldAlleleMap是 allele和hap-list的map，keyset就是allele
        final B[] newAlleles = newToOldAlleleMap.keySet().toArray((B[]) new Allele[newToOldAlleleMap.size()]);
        final int oldAlleleCount = alleles.numberOfAlleles();
        final int newAlleleCount = newAlleles.length;

        // we get the index correspondence between new old -> new allele, -1 entries mean that the old
        // allele does not map to any new; supported but typically not the case.
        final int[] oldToNewAlleleIndexMap = oldToNewAlleleIndexMap(newToOldAlleleMap, oldAlleleCount, newAlleles);

        final int[][] readsToKeep = overlappingReadIndicesBySampleIndex(overlap);

        // We calculate the marginal likelihoods.
        // 计算marginal的likelihood
        final double[][][][] newLikelihoodValues = marginalLikelihoods(oldAlleleCount, newAlleleCount,oldToNewAlleleIndexMap,readsToKeep,readLikelihoods_upperbound,readLikelihoods_exact,moreRecompute, readPT, redoSquares);
        // 现在newLikelihoodValues是新的allele-read直接的likelihood

        final int sampleCount = samples.numberOfSamples();

        @SuppressWarnings({"rawtypes","unchecked"})
        final Object2IntMap<GATKRead>[] newReadIndexBySampleIndex = (Object2IntMap<GATKRead>[])new Object2IntMap[sampleCount];
        final GATKRead[][] newReadsBySampleIndex = new GATKRead[sampleCount][];

        // sample count是输入DNA的数量
        // newSampleReadCount是保留下来的read的数量
        // 这一步大概只是在赋值？
        for (int s = 0; s < sampleCount; s++) {
            final int[] sampleReadsToKeep = readsToKeep[s];
            final GATKRead[] oldSampleReads = readsBySampleIndex[s];
            final int oldSampleReadCount = oldSampleReads.length;
            final int newSampleReadCount = sampleReadsToKeep.length;
            if (newSampleReadCount == oldSampleReadCount) {
                newReadsBySampleIndex[s] = oldSampleReads.clone();
            } else {
                newReadsBySampleIndex[s] = new GATKRead[newSampleReadCount];
                for (int i = 0; i < newSampleReadCount; i++) {
                    newReadsBySampleIndex[s][i] = oldSampleReads[sampleReadsToKeep[i]];
                }
            }
        }
        //-------

        // 这是一个三维矩阵 记录是否需要重算（猜测是 [sample index][allele index][read index]）
        boolean[][][] newrecomputeRecord = new boolean[sampleCount][][];
        // 一个三维矩阵，记录了某个sample的 hap-read矩阵，每个元素是该allele-read最高的下界
        int [][][] newbestHapMapLO = new int[sampleCount][][];
        // 一个三维矩阵，记录了某个sample的 hap-read矩阵，每个元素是该allele-read最高的上界
        int [][][] newbestHapMapUP = new int[sampleCount][][];

        // 给三个矩阵赋值
        for (int s = 0; s < sampleCount; s++) {
            newbestHapMapLO[s] = bestHapMapLO[s].clone();
            newbestHapMapUP[s] = bestHapMapUP[s].clone();
            newrecomputeRecord[s] = recomputeRecord[s].clone();
        }

        // 下面三个ReadLikelihoods变量分别由lo up exact生成 为什么依然要三个？
        // 最后的genotype表里的值是否要为精确值？bound-method只能生成估计值，为什么还能用？

        // Finally we create the new read-likelihood
        //System.err.print("leaving marginalize\n");
        final  List<ReadLikelihoods<B>> result= new ArrayList<ReadLikelihoods<B>>();
        result.add(new ReadLikelihoods<>(
                haplotypes,
                new IndexedAlleleList<>(newAlleles),
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues[0],
                newrecomputeRecord,
                newbestHapMapLO,
                newbestHapMapUP));
        result.add(new ReadLikelihoods<>(
                haplotypes,
                new IndexedAlleleList<>(newAlleles),
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues[1],
                newrecomputeRecord,
                newbestHapMapLO,
                newbestHapMapUP));
        result.add(new ReadLikelihoods<>(
                haplotypes,
                new IndexedAlleleList<>(newAlleles),
                samples,
                newReadsBySampleIndex,
                newReadIndexBySampleIndex,
                newLikelihoodValues[2],
                newrecomputeRecord,
                newbestHapMapLO,
                newbestHapMapUP));

        return result;
    }

    private int[][] overlappingReadIndicesBySampleIndex(final Locatable overlap) {
        if (overlap == null) {
            return null;
        }
        final int sampleCount = samples.numberOfSamples();
        final int[][] result = new int[sampleCount][];
        final IntArrayList buffer = new IntArrayList(200);

        for (int s = 0; s < sampleCount; s++) {
            buffer.clear();
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final int sampleReadCount = sampleReads.length;
            buffer.ensureCapacity(sampleReadCount);
            for (int r = 0; r < sampleReadCount; r++) {
                if (sampleReads[r].overlaps(overlap)) {
                    buffer.add(r);
                }
            }
            result[s] = buffer.toIntArray();
        }
        return result;
    }

    // Calculate the marginal likelihoods considering the old -> new allele index mapping.
    //Original
    private double[][][] marginalLikelihoods(final int oldAlleleCount, final int newAlleleCount, final int[] oldToNewAlleleIndexMap, final int[][] readsToKeep) {

        final int sampleCount = samples.numberOfSamples();
        final double[][][] result = new double[sampleCount][][];

        for (int s = 0; s < sampleCount; s++) {
            final int sampleReadCount = readsBySampleIndex[s].length;
            final double[][] oldSampleValues = valuesBySampleIndex[s];
            final int[] sampleReadToKeep = readsToKeep == null || readsToKeep[s].length == sampleReadCount ? null : readsToKeep[s];
            final int newSampleReadCount = sampleReadToKeep == null ? sampleReadCount : sampleReadToKeep.length;
            final double[][] newSampleValues = result[s] = new double[newAlleleCount][newSampleReadCount];
            // We initiate all likelihoods to -Inf.
            for (int a = 0; a < newAlleleCount; a++) {
                Arrays.fill(newSampleValues[a], Double.NEGATIVE_INFINITY);
            }
            // For each old allele and read we update the new table keeping the maximum likelihood.
            for (int r = 0; r < newSampleReadCount; r++) {
                for (int a = 0; a < oldAlleleCount; a++) {
                    final int oldReadIndex = newSampleReadCount == sampleReadCount ? r : sampleReadToKeep[r];
                    final int newAlleleIndex = oldToNewAlleleIndexMap[a];
                    if (newAlleleIndex == -1) {
                        continue;
                    }
                    final double likelihood = oldSampleValues[a][oldReadIndex];
                    if (likelihood > newSampleValues[newAlleleIndex][r]) {
                        newSampleValues[newAlleleIndex][r] = likelihood;
                    }
                }
            }
        }
        return result;
    }
    // Calculate the marginal likelihoods considering the old -> new allele index mapping.
    //Prune method
    private double[][][][] marginalLikelihoods(final int oldAlleleCount, final int newAlleleCount, final int[] oldToNewAlleleIndexMap, final int[][] readsToKeep, final ReadLikelihoods<Haplotype> readLikelihoods_upperbound, final ReadLikelihoods<Haplotype> readLikelihoods_exact, final boolean moreRecompute, int[] readPT, int [] redoSquares) {

        final int sampleCount = samples.numberOfSamples();
        final double[][][][] result = new double[3][sampleCount][][];//[0] is lowerbound, [1] is upperbound,[3] is exact (for debug)
        final int readIncrease = 1;
        redoSquares[0]=0;
        for (int s = 0; s < sampleCount; s++) {
            final int sampleReadCount = readsBySampleIndex[s].length;
            final double[][] oldSampleValues = valuesBySampleIndex[s];
            final int[] sampleReadToKeep = readsToKeep == null || readsToKeep[s].length == sampleReadCount ? null : readsToKeep[s];
            final int newSampleReadCount = sampleReadToKeep == null ? sampleReadCount : sampleReadToKeep.length;
            //
            final double[][] newSampleValues_lo = result[0][s] = new double[newAlleleCount][newSampleReadCount];
            final double[][] newSampleValues_up = result[1][s] = new double[newAlleleCount][newSampleReadCount];
            final double[][] newSampleValues_exact = result[2][s] = new double[newAlleleCount][newSampleReadCount];
            // 上面是在初始化矩阵

            //We recompute here if instructed by outer loop
            //final int readPTEnd = Math.min(newSampleReadCount,readPT[0]+readIncrease);
            //System.err.print("printlowerbound\n");
            //printlikelihoods();
            //System.err.print("printupperbound\n");
            //readLikelihoods_upperbound.printlikelihoods();

            // 如果需要重算
            if(moreRecompute){
                //for (int readCount=readPT[0];readCount<readPTEnd;readCount++){
                //System.err.printf("moreRecompute r=%d\n",r);
                boolean nextRound = true;
                //find the read value with largest upperbound-lowerbound gap
                int [] worstReadIndex = new int[newAlleleCount];
                double [] worstBoundGap = new double[newAlleleCount];
                Arrays.fill(worstReadIndex, 0);
                Arrays.fill(worstBoundGap, 0);
                //exact_only[0] = true;
                for(int a=0; a<newAlleleCount; a++){
                    for(int r = 0; r<newSampleReadCount; r++){
                        final int oldReadIndex = newSampleReadCount == sampleReadCount ? r : sampleReadToKeep[r];
                        final int bestOldAlleleIndexLO = bestHapMapLO[s][a][r];
                        final int bestOldAlleleIndexUP = bestHapMapUP[s][a][r];
                        double boundGap = readLikelihoods_upperbound.valuesBySampleIndex[s][bestOldAlleleIndexLO][oldReadIndex]-valuesBySampleIndex[s][bestOldAlleleIndexUP][oldReadIndex];
                        //System.err.printf("ReadLikelihoods: a=%d r =%d gap=%f\n",a,r,boundGap);
                        // 这一步得到了最糟糕的gap所对应的read的index
                        if(boundGap>worstBoundGap[a]){
                            worstBoundGap[a] = boundGap;
                            worstReadIndex[a] = r;
                        }
                    }
                    //System.err.printf("ReadLikelihoods: worstReadIndex[%d]=%d gap=%f\n", a,worstReadIndex[a],worstBoundGap[a] );
                    //exact_only[0] &= (worstBoundGap[a]==0);
                }

                //System.err.printf("ReadLikelihoods: exact_only=%b\n", exact_only[0] );
                do{
                    nextRound = true;
                    for(int newA=0; newA<newAlleleCount; newA++){
                        // r是最糟糕的那个read的index
                        int r = worstReadIndex[newA];
                        final int oldReadIndex = newSampleReadCount == sampleReadCount ? r : sampleReadToKeep[r];
                        final int bestOldAlleleIndexLO = bestHapMapLO[s][newA][r];
                        final int bestOldAlleleIndexUP = bestHapMapUP[s][newA][r];

                        //do the haplotype corresponding to upperbound
                        if(!recomputeRecord[s][bestOldAlleleIndexUP][oldReadIndex]){
                            // 如果没有重算过 就重算
                            // TODO 这里应该要真实地重算 把exact的部分代替掉
                            //System.err.printf("ReadLikelihoods: UP: newA=%d readIndex=%d hapIndex=%d totHap=%d totRead=%d gap=%f\n",newA,oldReadIndex,bestOldAlleleIndexUP,oldAlleleCount,newSampleReadCount,readLikelihoods_upperbound.valuesBySampleIndex[s][bestOldAlleleIndexUP][oldReadIndex]-valuesBySampleIndex[s][bestOldAlleleIndexLO][oldReadIndex]);
                            valuesBySampleIndex[s][bestOldAlleleIndexUP][oldReadIndex]
                                    =readLikelihoods_exact.valuesBySampleIndex[s][bestOldAlleleIndexUP][oldReadIndex];
                            readLikelihoods_upperbound.valuesBySampleIndex[s][bestOldAlleleIndexUP][oldReadIndex]
                                    =readLikelihoods_exact.valuesBySampleIndex[s][bestOldAlleleIndexUP][oldReadIndex];
                            redoSquares[0]+=readsBySampleIndex[s][oldReadIndex].getLength()*alleles.getAllele(bestOldAlleleIndexUP).getBases().length;
                            recomputeRecord[s][bestOldAlleleIndexUP][oldReadIndex]=true;
                            readLikelihoods_upperbound.recomputeRecord[s][bestOldAlleleIndexUP][oldReadIndex]=true;
                            readLikelihoods_exact.recomputeRecord[s][bestOldAlleleIndexUP][oldReadIndex]=true;
                        }

                    }
                    for (int a = 0; a < newAlleleCount; a++) {
                        int r = worstReadIndex[a];
                        newSampleValues_lo[a][r]=Double.NEGATIVE_INFINITY;
                        newSampleValues_up[a][r]=Double.NEGATIVE_INFINITY;
                        newSampleValues_exact[a][r]=Double.NEGATIVE_INFINITY;
                    }
                    //System.err.printf("done this round moreRecompute\n");
                    for(int allele = 0; allele<newAlleleCount; allele ++){
                        int r = worstReadIndex[allele];
                        final int oldReadIndex = newSampleReadCount == sampleReadCount ? r : sampleReadToKeep[r];
                        int newbestHapMapLO=0;
                        int newbestHapMapUP=0;
                        //check for bestAllele
                        // 找到现在的likelihood中最高值对应的haplotype的index
                        for (int a = 0; a < oldAlleleCount; a++) {
                            final int newAlleleIndex = oldToNewAlleleIndexMap[a];
                            if (newAlleleIndex !=allele) {
                                continue;
                            }
                            double likelihood = valuesBySampleIndex[s][a][oldReadIndex]; //this is lowerbound
                            if (likelihood > newSampleValues_lo[newAlleleIndex][r]) {
                                newSampleValues_lo[newAlleleIndex][r] = likelihood;
                                newbestHapMapLO = a;
                            }
                            likelihood = readLikelihoods_upperbound.valuesBySampleIndex[s][a][oldReadIndex]; //this is upperbound
                            if (likelihood > newSampleValues_up[newAlleleIndex][r]) {
                                newSampleValues_up[newAlleleIndex][r] = likelihood;
                                newbestHapMapUP = a;
                            }
                        }
                        // 检查得到的index和原本的index是否一样
                        boolean isExact = newbestHapMapLO==bestHapMapLO[s][allele][r] && newbestHapMapUP==bestHapMapUP[s][allele][r];
                        nextRound &= isExact;
                        bestHapMapLO[s][allele][r] = newbestHapMapLO;
                        bestHapMapUP[s][allele][r] = newbestHapMapUP;
                        //System.err.printf("moreRecompute allele=%d r=%d nextRound=%b\n",allele,r,nextRound);
                    }
                }while(!nextRound);
            }

            // 如果不重算直接跳到这里
            // We initiate all likelihoods to -Inf.
            for (int a = 0; a < newAlleleCount; a++) {
                Arrays.fill(newSampleValues_lo[a], Double.NEGATIVE_INFINITY);
                Arrays.fill(newSampleValues_up[a], Double.NEGATIVE_INFINITY);
                Arrays.fill(newSampleValues_exact[a], Double.NEGATIVE_INFINITY);
            }

            // 从hap表格提到allele表格时
            // lo和up都取最大的那个值
            // For each old allele and read we update the new table keeping the maximum likelihood.
            ////WE take maximum of upperbound and lowerbound now
            // 这里的两个矩阵是index的映射
            int [][] newbestHapMapLO = new int[newAlleleCount][newSampleReadCount];
            int [][] newbestHapMapUP = new int[newAlleleCount][newSampleReadCount];
            for (int r = 0; r < newSampleReadCount; r++) {
                final int oldReadIndex = newSampleReadCount == sampleReadCount ? r : sampleReadToKeep[r];
                for (int a = 0; a < oldAlleleCount; a++) {
                    final int newAlleleIndex = oldToNewAlleleIndexMap[a];
                    if (newAlleleIndex == -1) {
                        continue;
                    }
                    double likelihood = readLikelihoods_exact.valuesBySampleIndex[s][a][oldReadIndex]; //this is exact
                    if (likelihood > newSampleValues_exact[newAlleleIndex][r]) {
                        newSampleValues_exact[newAlleleIndex][r] = likelihood;
                    }
                    likelihood = valuesBySampleIndex[s][a][oldReadIndex]; //this is lowerbound
                    if (likelihood > newSampleValues_lo[newAlleleIndex][r]) {
                        newSampleValues_lo[newAlleleIndex][r] = likelihood;
                        newbestHapMapLO[newAlleleIndex][r] = a;
                    }
                    likelihood = readLikelihoods_upperbound.valuesBySampleIndex[s][a][oldReadIndex]; //this is upperbound
                    if (likelihood > newSampleValues_up[newAlleleIndex][r]) {
                        newSampleValues_up[newAlleleIndex][r] = likelihood;
                        newbestHapMapUP[newAlleleIndex][r] = a;
                    }
                }

            }
            // 这两个变量是 最好的那个值所对应的原本的hap的index
            bestHapMapLO[s] = newbestHapMapLO.clone();
            bestHapMapUP[s] = newbestHapMapUP.clone();
        }
        return result;
    }

    // calculates an old to new allele index map array.
    // 首先old和new都是allele的东西 它们之间需要一个map
    // newToOldAlleleMap: allele -> hap-list 的map
    private <B extends Allele> int[] oldToNewAlleleIndexMap(final Map<B, List<A>> newToOldAlleleMap, final int oldAlleleCount, final B[] newAlleles) {
        Arrays.stream(newAlleles).forEach(Utils::nonNull);
        Utils.containsNoNull(newToOldAlleleMap.values(), "no new allele list can be null");

        newToOldAlleleMap.values().stream().forEach(oldList -> Utils.containsNoNull(oldList,"old alleles cannot be null"));

        final int[] oldToNewAlleleIndexMap = new int[oldAlleleCount];
        Arrays.fill(oldToNewAlleleIndexMap, -1);  // -1 indicate that there is no new allele that make reference to that old one.

        for (int newIndex = 0; newIndex < newAlleles.length; newIndex++) {
            // 对newAllele进行遍历(newAlleles就是mapper的key_set)
            final B newAllele = newAlleles[newIndex];
            for (final A oldAllele : newToOldAlleleMap.get(newAllele)) {
                // 这一步有疑问??
                // oldAllele是mapper中key为newAllele的那个list中的每一个元素
                // 这一步得到了该元素在alleles对象中的index
                final int oldAlleleIndex = indexOfAllele(oldAllele);
                if (oldAlleleIndex == -1) {
                    throw new IllegalArgumentException("missing old allele " + oldAllele + " in likelihood collection ");
                }
                if (oldToNewAlleleIndexMap[oldAlleleIndex] != -1) {
                    throw new IllegalArgumentException("collision: two new alleles make reference to the same old allele");
                }
                oldToNewAlleleIndexMap[oldAlleleIndex] = newIndex;
            }
        }
        return oldToNewAlleleIndexMap;
    }

    /**
     * Removes those read that the best possible likelihood given any allele is just too low.
     *
     * <p>
     *     This is determined by a maximum error per read-base against the best likelihood possible.
     * </p>
     *
     * @param maximumErrorPerBase the minimum acceptable error rate per read base, must be
     *                            a positive number.
     *
     * @throws IllegalStateException is not supported for read-likelihood that do not contain alleles.
     *
     * @throws IllegalArgumentException if {@code maximumErrorPerBase} is negative.
     */
    //Prune method
    // Chenhao: new parameter
    public void filterPoorlyModeledReads(final double maximumErrorPerBase, ReadLikelihoods<Haplotype> result_upperbound, ReadLikelihoods<Haplotype> result_exact, List<Map<GATKRead, byte[]>> gapPenalties, List<List<GATKRead>> processedReadList) {
        Utils.validateArg(alleles.numberOfAlleles() > 0, "unsupported for read-likelihood collections with no alleles");
        Utils.validateArg(!Double.isNaN(maximumErrorPerBase) && maximumErrorPerBase > 0.0, "the maximum error per base must be a positive number");

        new IndexRange(0, samples.numberOfSamples()).forEach(s -> {
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            // added by Chenhao
            Map<GATKRead, byte[]> penalty = gapPenalties.get(s);
            List<GATKRead> processedReads = processedReadList.get(s);
            // Debug
            System.out.println("gapPenalty: " + gapPenalties);
            System.out.println("penalty: " + penalty);
            // added by Chenhao: new parameter penalty
            final List<Integer> removeIndices = new IndexRange(0, sampleReads.length).filter(r -> readIsPoorlyModelled(s, r, sampleReads[r], maximumErrorPerBase, result_upperbound,result_exact, penalty, processedReads.get(r)));

            removeSampleReads(s, removeIndices, alleles.numberOfAlleles());
            result_upperbound.removeSampleReads(s, removeIndices, alleles.numberOfAlleles());
            result_exact.removeSampleReads(s, removeIndices, alleles.numberOfAlleles());

            // added by Chenhao
            System.out.println("--------------------");
            System.out.println("Removed index: " + removeIndices);
        });
    }
    //Original
    public void filterPoorlyModeledReads(final double maximumErrorPerBase) {
        Utils.validateArg(alleles.numberOfAlleles() > 0, "unsupported for read-likelihood collections with no alleles");
        Utils.validateArg(!Double.isNaN(maximumErrorPerBase) && maximumErrorPerBase > 0.0, "the maximum error per base must be a positive number");

        new IndexRange(0, samples.numberOfSamples()).forEach(s -> {
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final List<Integer> removeIndices = new IndexRange(0, sampleReads.length)
                    .filter(r -> readIsPoorlyModelled(s, r, sampleReads[r], maximumErrorPerBase));
            removeSampleReads(s, removeIndices, alleles.numberOfAlleles());

            // added by Chenhao
            System.out.println("--------------------");
            System.out.println("Region index: " + s);
            System.out.println("Removed index: " + removeIndices);
        });

    }

    //Original
    private boolean readIsPoorlyModelled(final int sampleIndex, final int readIndex, final GATKRead read, final double maxErrorRatePerBase) {
        final double maxErrorsForRead = Math.min(2.0, Math.ceil(read.getLength() * maxErrorRatePerBase));
        final double log10QualPerBase = -4.0;
        final double log10MaxLikelihoodForTrueAllele = maxErrorsForRead * log10QualPerBase;

        final int alleleCount = alleles.numberOfAlleles();
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][readIndex] >= log10MaxLikelihoodForTrueAllele) {
                return false;
            }
        }
        return true;
    }

    /// added by Chenhao
    long num_reads = 0;
    long saved_by_lower = 0;
    long delete_by_upper = 0;
    long unknown = 0;
    long original_work = 0;
    long filter_work = 0;

    //Prune method
    // added by Chenhao: new parameter
    private boolean readIsPoorlyModelled(final int sampleIndex, final int readIndex, final GATKRead read, final double maxErrorRatePerBase,ReadLikelihoods<Haplotype> result_upperbound, ReadLikelihoods<Haplotype> result_exact, final Map<GATKRead, byte[]> penalty, final GATKRead processedRead) {
        final double maxErrorsForRead = Math.min(2.0, Math.ceil(read.getLength() * maxErrorRatePerBase));
        final double log10QualPerBase = -4.0;
        final double log10MaxLikelihoodForTrueAllele = maxErrorsForRead * log10QualPerBase;

        final int alleleCount = alleles.numberOfAlleles();
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        boolean rm_lowerbound=true;
        boolean rm_upperbound=true;
        boolean rm_exact = true;

        // this part for debug
        if (true){
            int recompute_squares = 0;
            final boolean isFirstHaplotype = true;
            for (int a = 0; a < alleleCount; a++){
                System.out.println("recompute for hap index: " + a);
                final List<Haplotype> my_alleles = result_upperbound.sampleMatrix(sampleIndex).alleles();
                final byte[] nextAlleleBases = a == my_alleles.size() - 1 ? null : my_alleles.get(a + 1).getBases();
                // check penalty
                System.out.println("penalty: " + penalty.get(processedRead));
                double exact_readLikelihood = recompute_tool.recomputeExactLog10Likelihoods(sampleIndex, processedRead, a, result_upperbound, penalty, isFirstHaplotype, nextAlleleBases);
                // Debug: check the result first
                if (exact_readLikelihood != result_exact.valuesBySampleIndex[sampleIndex][a][readIndex]){
                    System.out.println("some error");
                }

                if (exact_readLikelihood >= log10MaxLikelihoodForTrueAllele){
                    rm_exact = false;
                }
                recompute_squares += read.getLength()*alleles.getAllele(a).getBases().length;
                filter_work += recompute_squares;
                valuesBySampleIndex[sampleIndex][a][readIndex] = exact_readLikelihood;
                result_upperbound.valuesBySampleIndex[sampleIndex][a][readIndex] = exact_readLikelihood;
            }
        }


        for (int a = 0; a < alleleCount; a++) {
            if (sampleValues[a][readIndex] >= log10MaxLikelihoodForTrueAllele) {
                rm_lowerbound=false;
            }
        }
        for (int a = 0; a < alleleCount; a++) {
            if (result_upperbound.valuesBySampleIndex[sampleIndex][a][readIndex] >= log10MaxLikelihoodForTrueAllele) {
                rm_upperbound=false;
            }
        }

        // added by Chenhao
        num_reads = num_reads + 1;
        for (int a = 0; a < alleleCount; a++) {
            original_work += read.getLength()*alleles.getAllele(a).getBases().length;
        }

        if(rm_upperbound) {
            // added by Chenhao
            delete_by_upper = delete_by_upper + 1;

            //System.err.printf("Xiao: yes remove read from utils/genotyper/ReadLikelihoods.java/readIsPoorlyModelled index=%d\n",readIndex);
            return true;
        }else if(!rm_lowerbound) {
            // added by Chenhao
            saved_by_lower = saved_by_lower + 1;

            //System.err.printf("Xiao: no remove read from utils/genotyper/ReadLikelihoods.java/readIsPoorlyModelled index=%d\n",readIndex);
            return false;
        }else{
            // added by Chenhao
            unknown = unknown + 1;
            System.err.printf("Xiao: redo from utils/genotyper/ReadLikelihoods.java/readIsPoorlyModelled index=%d\n",readIndex);
            int recompute_squares = 0;
            // added by Chenhao: the following is the original prune version
            //
            //for (int a = 0; a < alleleCount; a++){
            //    if (result_exact.valuesBySampleIndex[sampleIndex][a][readIndex] >= log10MaxLikelihoodForTrueAllele) {
            //        rm_exact = false;
            //    }
            //    recompute_squares += read.getLength()*alleles.getAllele(a).getBases().length;
            //    valuesBySampleIndex[sampleIndex][a][readIndex] = result_exact.valuesBySampleIndex[sampleIndex][a][readIndex];
            //    result_upperbound.valuesBySampleIndex[sampleIndex][a][readIndex] = result_exact.valuesBySampleIndex[sampleIndex][a][readIndex];
            //}

            // added by Chenhao: the following is the recompute version
            final boolean isFirstHaplotype = true;
            for (int a = 0; a < alleleCount; a++){
                final List<Haplotype> my_alleles = result_upperbound.sampleMatrix(sampleIndex).alleles();
                final byte[] nextAlleleBases = a == my_alleles.size() - 1 ? null : my_alleles.get(a + 1).getBases();
                double exact_readLikelihood = recompute_tool.recomputeExactLog10Likelihoods(sampleIndex, processedRead, a, result_upperbound, penalty, isFirstHaplotype, nextAlleleBases);
                // Debug: check the result first
                if (exact_readLikelihood >= log10MaxLikelihoodForTrueAllele){
                    rm_exact = false;
                }
                recompute_squares += read.getLength()*alleles.getAllele(a).getBases().length;
                filter_work += recompute_squares;
                valuesBySampleIndex[sampleIndex][a][readIndex] = exact_readLikelihood;
                result_upperbound.valuesBySampleIndex[sampleIndex][a][readIndex] = exact_readLikelihood;
            }
            // added by Chenhao: end of recompute version
            System.err.printf("Xiao:redo total squares = %d reason = filterReads\n", recompute_squares);
            return rm_exact;
        }
    }


    /**
     * Add more reads to the collection.
     *
     * @param readsBySample reads to add.
     * @param initialLikelihood the likelihood for the new entries.
     *
     * @throws IllegalArgumentException if {@code readsBySample} is {@code null} or {@code readsBySample} contains
     *  {@code null} reads, or {@code readsBySample} contains read that are already present in the read-likelihood
     *  collection.
     */
    public void addReads(final Map<String,List<GATKRead>> readsBySample, final double initialLikelihood) {
        for (final Map.Entry<String,List<GATKRead>> entry : readsBySample.entrySet()) {
            final String sample = entry.getKey();
            final List<GATKRead> newSampleReads = entry.getValue();
            final int sampleIndex = samples.indexOfSample(sample);

            if (sampleIndex == -1) {
                throw new IllegalArgumentException("input sample " + sample +
                        " is not part of the read-likelihoods collection");
            }

            if (newSampleReads == null || newSampleReads.isEmpty()) {
                continue;
            }

            final int sampleReadCount = readsBySampleIndex[sampleIndex].length;
            final int newSampleReadCount = sampleReadCount + newSampleReads.size();

            appendReads(newSampleReads, sampleIndex, sampleReadCount, newSampleReadCount);
            extendsLikelihoodArrays(initialLikelihood, sampleIndex, sampleReadCount, newSampleReadCount);
        }
    }

    // Extends the likelihood arrays-matrices.
    private void extendsLikelihoodArrays(final double initialLikelihood, final int sampleIndex, final int sampleReadCount, final int newSampleReadCount) {
        final double[][] sampleValues = valuesBySampleIndex[sampleIndex];
        final int alleleCount = alleles.numberOfAlleles();
        for (int a = 0; a < alleleCount; a++) {
            sampleValues[a] = Arrays.copyOf(sampleValues[a], newSampleReadCount);
        }
        if (initialLikelihood != 0.0) // the default array new value.
        {
            for (int a = 0; a < alleleCount; a++) {
                Arrays.fill(sampleValues[a], sampleReadCount, newSampleReadCount, initialLikelihood);
            }
        }
    }

    // Append the new read reference into the structure per-sample.
    private void appendReads(final List<GATKRead> newSampleReads, final int sampleIndex,
                             final int sampleReadCount, final int newSampleReadCount) {
        final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex] =
                Arrays.copyOf(readsBySampleIndex[sampleIndex], newSampleReadCount);

        int nextReadIndex = sampleReadCount;
        final Object2IntMap<GATKRead> sampleReadIndex = readIndexBySampleIndex[sampleIndex];
        for (final GATKRead newRead : newSampleReads) {
            //    if (sampleReadIndex.containsKey(newRead)) // might be worth handle this without exception (ignore the read?) but in practice should never be the case.
            //        throw new IllegalArgumentException("you cannot add reads that are already in read-likelihood collection");
            if (sampleReadIndex != null ) {
                sampleReadIndex.put(newRead, nextReadIndex);
            }
            sampleReads[nextReadIndex++] = newRead;
        }
    }

    /**
     * Adds the non-reference allele to the read-likelihood collection setting each read likelihood to the second
     * best found (or best one if only one allele has likelihood).
     *
     * <p>Nothing will happen if the read-likelihoods collection already includes the non-ref allele</p>
     *
     * <p>
     *     <i>Implementation note: even when strictly speaking we do not need to demand the calling code to pass
     *     the reference the non-ref allele, we still demand it in order to lead the
     *     the calling code to use the right generic type for this likelihoods
     *     collection {@link Allele}.</i>
     * </p>
     *
     * @param nonRefAllele the non-ref allele.
     *
     * @throws IllegalArgumentException if {@code nonRefAllele} is anything but the designated &lt;NON_REF&gt;
     * symbolic allele {@link Allele#NON_REF_ALLELE}.
     */
    public void addNonReferenceAllele(final A nonRefAllele) {
        Utils.nonNull(nonRefAllele, "non-ref allele cannot be null");
        if (!nonRefAllele.equals(Allele.NON_REF_ALLELE)) {
            throw new IllegalArgumentException("the non-ref allele is not valid");
        } else if (alleles.containsAllele(nonRefAllele)) {
            return;
        } else if (addMissingAlleles(Collections.singleton(nonRefAllele), Double.NEGATIVE_INFINITY)) {
            updateNonRefAlleleLikelihoods();
        }
    }

    /**
     * Updates the likelihoods of the non-ref allele, if present, considering all non-symbolic alleles avaialble.
     */
    public void updateNonRefAlleleLikelihoods() {
        updateNonRefAlleleLikelihoods(alleles);
    }

    /**
     * Updates the likelihood of the NonRef allele (if present) based on the likelihoods of a set of non-symbolic
     * <p>
     *     This method does
     * </p>
     *
     *
     * @param allelesToConsider
     */
    @SuppressWarnings("unchecked")  // for the cast (A) Allele.NON_REF_ALLELE below
    public void updateNonRefAlleleLikelihoods(final AlleleList<A> allelesToConsider) {
        final int nonRefAlleleIndex = indexOfAllele((A) Allele.NON_REF_ALLELE);
        if ( nonRefAlleleIndex < 0) {
            return;
        }
        final int alleleCount = alleles.numberOfAlleles();
        final int nonSymbolicAlleleCount = alleleCount - 1;
        // likelihood buffer reused across reads:
        final double[] qualifiedAlleleLikelihoods = new double[nonSymbolicAlleleCount];
        final Median medianCalculator = new Median();
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final double[][] sampleValues = valuesBySampleIndex[s];
            final int readCount = sampleValues[0].length;
            for (int r = 0; r < readCount; r++) {
                final BestAllele bestAllele = searchBestAllele(s, r, true, false);
                int numberOfQualifiedAlleleLikelihoods = 0;
                for (int i = 0; i < alleleCount; i++) {
                    final double alleleLikelihood = sampleValues[i][r];
                    if (i != nonRefAlleleIndex && alleleLikelihood < bestAllele.likelihood
                            && !Double.isNaN(alleleLikelihood) && allelesToConsider.indexOfAllele(alleles.getAllele(i)) != -1) {
                        qualifiedAlleleLikelihoods[numberOfQualifiedAlleleLikelihoods++] = alleleLikelihood;
                    }
                }
                final double nonRefLikelihood = medianCalculator.evaluate(qualifiedAlleleLikelihoods, 0, numberOfQualifiedAlleleLikelihoods);
                // when the median is NaN that means that all applicable likekihoods are the same as the best
                // so the read is not informative at all given the existing alleles. Unless there is only one (or zero) concrete
                // alleles with give the same (the best) likelihood to the NON-REF. When there is only one (or zero) concrete
                // alleles we set the NON-REF likelihood to NaN.
                sampleValues[nonRefAlleleIndex][r] = !Double.isNaN(nonRefLikelihood) ? nonRefLikelihood
                        : nonSymbolicAlleleCount <= 1 ? Double.NaN : bestAllele.likelihood;
            }
        }
    }



    /**
     * Downsamples reads based on contamination fractions making sure that all alleles are affected proportionally.
     *
     * @param perSampleDownsamplingFraction contamination sample map where the sample name are the keys and the
     *                                       fractions are the values.
     *
     * @throws IllegalArgumentException if {@code perSampleDownsamplingFraction} is {@code null}.
     */
    public void contaminationDownsampling(final Map<String, Double> perSampleDownsamplingFraction) {
        Utils.nonNull(perSampleDownsamplingFraction);

        final int alleleCount = alleles.numberOfAlleles();
        for (int s = 0; s < samples.numberOfSamples(); s++) {
            final String sample = samples.getSample(s);
            final Double fractionDouble = perSampleDownsamplingFraction.get(sample);
            if (fractionDouble == null) {
                continue;
            }
            final double fraction = fractionDouble;
            if (Double.isNaN(fraction) || fraction <= 0.0) {
                continue;
            }
            if (fraction >= 1.0) {
                final List<Integer> removeIndices = IntStream.range(0, readsBySampleIndex[s].length).boxed().collect(Collectors.toList());
                removeSampleReads(s, removeIndices, alleleCount);
            } else {
                final Map<A,List<GATKRead>> readsByBestAllelesMap = readsByBestAlleleMap(s);
                removeSampleReads(s, AlleleBiasedDownsamplingUtils.selectAlleleBiasedReads(readsByBestAllelesMap, fraction),alleleCount);
            }
        }
    }

    /**
     * Returns the collection of best allele estimates for the reads based on the read-likelihoods.
     * "Ties" where the ref likelihood is within {@code ReadLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken in favor of the reference.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per read in the read-likelihoods collection.
     */
    public Collection<BestAllele> bestAllelesBreakingTies() {
        return IntStream.range(0, numberOfSamples()).boxed().flatMap(n -> bestAllelesBreakingTies(n).stream()).collect(Collectors.toList());
    }

    /**
     * Returns the collection of best allele estimates for one sample's reads based on the read-likelihoods.
     * "Ties" where the ref likelihood is within {@code ReadLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken in favor of the reference.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per read in the read-likelihoods collection.
     */
    public Collection<BestAllele> bestAllelesBreakingTies(final String sample) {
        final int sampleIndex = indexOfSample(sample);
        return bestAllelesBreakingTies(sampleIndex);
    }

    /**
     * Returns the collection of best allele estimates for one sample's reads reads based on the read-likelihoods.
     * "Ties" where the ref likelihood is within {@code ReadLikelihoods.INFORMATIVE_THRESHOLD} of the greatest likelihood
     * are broken in favor of the reference.
     *
     * @throws IllegalStateException if there is no alleles.
     *
     * @return never {@code null}, one element per read in the read-likelihoods collection.
     */
    private Collection<BestAllele> bestAllelesBreakingTies(final int sampleIndex) {
        Utils.validIndex(sampleIndex, numberOfSamples());

        final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
        final int readCount = sampleReads.length;
        final List<BestAllele> result = new ArrayList<>(readCount);
        for (int r = 0; r < readCount; r++) {
            result.add(searchBestAllele(sampleIndex, r, true, true));
        }

        return result;
    }


    /**
     * Returns reads stratified by their best allele.
     * @param sampleIndex the target sample.
     * @return never {@code null}, perhaps empty.
     */
    private Map<A,List<GATKRead>> readsByBestAlleleMap(final int sampleIndex) {
        Utils.validIndex(sampleIndex, numberOfSamples());
        final int alleleCount = alleles.numberOfAlleles();
        final int sampleReadCount = readsBySampleIndex[sampleIndex].length;
        final Map<A,List<GATKRead>> result = new LinkedHashMap<>(alleleCount);
        for (int a = 0; a < alleleCount; a++) {
            result.put(alleles.getAllele(a), new ArrayList<>(sampleReadCount));
        }
        readsByBestAlleleMap(sampleIndex,result);
        return result;
    }

    /**
     * Returns reads stratified by their best allele.
     * @return never {@code null}, perhaps empty.
     */
    @VisibleForTesting
    Map<A,List<GATKRead>> readsByBestAlleleMap() {
        final int alleleCount = alleles.numberOfAlleles();
        final Map<A,List<GATKRead>> result = new LinkedHashMap<>(alleleCount);
        final int totalReadCount = readCount();
        for (int a = 0; a < alleleCount; a++) {
            result.put(alleles.getAllele(a), new ArrayList<>(totalReadCount));
        }
        final int sampleCount = samples.numberOfSamples();
        for (int s = 0; s < sampleCount; s++) {
            readsByBestAlleleMap(s, result);
        }
        return result;
    }

    private void readsByBestAlleleMap(final int sampleIndex, final Map<A, List<GATKRead>> result) {
        final GATKRead[] reads = readsBySampleIndex[sampleIndex];
        final int readCount = reads.length;

        for (int r = 0; r < readCount; r++) {
            final BestAllele bestAllele = searchBestAllele(sampleIndex,r,true, false);
            if (!bestAllele.isInformative()) {
                continue;
            }
            result.get(bestAllele.allele).add(bestAllele.read);
        }
    }

    /**
     * Returns the index of a read within a sample read-likelihood sub collection.
     * @param sampleIndex the sample index.
     * @param read the query read.
     * @return -1 if there is no such read in that sample, 0 or greater otherwise.
     */
    @VisibleForTesting
    int readIndex(final int sampleIndex, final GATKRead read) {
        final Object2IntMap<GATKRead> readIndex = readIndexBySampleIndex(sampleIndex);
        if (readIndex.containsKey(read)) {
            return readIndexBySampleIndex(sampleIndex).getInt(read);
        } else {
            return -1;
        }
    }

    /**
     * Returns the total number of reads in the read-likelihood collection.
     *
     * @return never {@code null}
     */
    public int readCount() {
        int sum = 0;
        final int sampleCount = samples.numberOfSamples();
        for (int i = 0; i < sampleCount; i++) {
            sum += readsBySampleIndex[i].length;
        }
        return sum;
    }

    /**
     * Returns the number of reads that belong to a sample in the read-likelihood collection.
     * @param sampleIndex the query sample index.
     *
     * @throws IllegalArgumentException if {@code sampleIndex} is not a valid sample index.
     * @return 0 or greater.
     */
    public int sampleReadCount(final int sampleIndex) {
        Utils.validIndex(sampleIndex, samples.numberOfSamples());
        return readsBySampleIndex[sampleIndex].length;
    }

    /**
     * Remove those reads that do not overlap certain genomic location.
     *
     * <p>
     *     This method modifies the current read-likelihoods collection.
     * </p>
     *
     * @param location the target location.
     *
     * @throws IllegalArgumentException the location cannot be {@code null} nor unmapped.
     */
    public void filterToOnlyOverlappingReads(final SimpleInterval location) {
        Utils.nonNull(location, "the location cannot be null");

        final int sampleCount = samples.numberOfSamples();

        final int alleleCount = alleles.numberOfAlleles();
        for (int s = 0; s < sampleCount; s++) {
            final GATKRead[] sampleReads = readsBySampleIndex[s];
            final List<Integer> removeIndices = new IndexRange(0, sampleReads.length)
                    .filter(r -> !sampleReads[r].overlaps(location));
            removeSampleReads(s, removeIndices, alleleCount);
        }
    }

    /**
     * Contains information about the best allele for a read search result.
     */
    public final class BestAllele {
        public static final double INFORMATIVE_THRESHOLD = 0.2;

        /**
         * Null if there is no possible match (no allele?).
         */
        public final A allele;

        /**
         * The containing sample.
         */
        public final String sample;

        /**
         * The query read.
         */
        public final GATKRead read;

        /**
         * If allele != null, the indicates the likelihood of the read.
         */
        public final double likelihood;

        /**
         * Confidence that the read actually was generated under that likelihood.
         * This is equal to the difference between this and the second best allele match.
         */
        public final double confidence;

        private BestAllele(final int sampleIndex, final int readIndex, final int bestAlleleIndex,
                           final double likelihood, final double secondBestLikelihood) {
            allele = bestAlleleIndex == -1 ? null : alleles.getAllele(bestAlleleIndex);
            this.likelihood = likelihood;
            sample = samples.getSample(sampleIndex);
            read = readsBySampleIndex[sampleIndex][readIndex];
            confidence = likelihood == secondBestLikelihood ? 0 : likelihood - secondBestLikelihood;
        }

        public boolean isInformative() {
            return confidence > INFORMATIVE_THRESHOLD;
        }
    }

    private void removeSampleReads(final int sampleIndex, final List<Integer> removeIndices, final int alleleCount) {
        if (removeIndices.isEmpty()) {
            return;
        }

        final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
        final int sampleReadCount = sampleReads.length;

        final Object2IntMap<GATKRead> indexByRead = readIndexBySampleIndex[sampleIndex];
        if (indexByRead != null) {
            removeIndices.stream().forEach(n -> indexByRead.remove(sampleReads[n]));
        }
        final boolean[] removeIndex = new boolean[sampleReadCount];
        final int firstDeleted = removeIndices.get(0);
        removeIndices.stream().forEach(n -> removeIndex[n] = true);

        final int newSampleReadCount = sampleReadCount - removeIndices.size();

        // Now we skim out the removed reads from the read array.
        final GATKRead[] oldSampleReads = readsBySampleIndex[sampleIndex];
        final GATKRead[] newSampleReads = new GATKRead[newSampleReadCount];

        System.arraycopy(oldSampleReads, 0, newSampleReads, 0, firstDeleted);
        Utils.skimArray(oldSampleReads,firstDeleted, newSampleReads, firstDeleted, removeIndex, firstDeleted);

        // Then we skim out the likelihoods of the removed reads.
        final double[][] oldSampleValues = valuesBySampleIndex[sampleIndex];
        final double[][] newSampleValues = new double[alleleCount][newSampleReadCount];
        for (int a = 0; a < alleleCount; a++) {
            System.arraycopy(oldSampleValues[a], 0, newSampleValues[a], 0, firstDeleted);
            Utils.skimArray(oldSampleValues[a], firstDeleted, newSampleValues[a], firstDeleted, removeIndex, firstDeleted);
        }
        valuesBySampleIndex[sampleIndex] = newSampleValues;
        readsBySampleIndex[sampleIndex] = newSampleReads;
        readListBySampleIndex[sampleIndex] = null; // reset the unmodifiable list.
    }


    // Requires that the collection passed iterator can remove elements, and it can be modified.
    public void removeSampleReads(final int sampleIndex, final Collection<GATKRead> readsToRemove, final int alleleCount) {
        final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
        final int sampleReadCount = sampleReads.length;

        final Object2IntMap<GATKRead> indexByRead = readIndexBySampleIndex(sampleIndex);
        // Count how many we are going to remove, which ones (indexes) and remove entry from the read-index map.
        final boolean[] removeIndex = new boolean[sampleReadCount];
        int removeCount = 0; // captures the number of deletions.
        int firstDeleted = sampleReadCount;    // captures the first position that was deleted.

        final Iterator<GATKRead> readsToRemoveIterator = readsToRemove.iterator();
        while (readsToRemoveIterator.hasNext()) {
            final GATKRead read = readsToRemoveIterator.next();
            if (indexByRead.containsKey(read)) {
                final int index = indexByRead.getInt(read);
                if (firstDeleted > index) {
                    firstDeleted = index;
                }
                removeCount++;
                removeIndex[index] = true;
                readsToRemoveIterator.remove();
                indexByRead.remove(read);
            }
        }

        // Nothing to remove we just finish here.
        if (removeCount == 0) {
            return;
        }

        final int newSampleReadCount = sampleReadCount - removeCount;

        // Now we skim out the removed reads from the read array.
        final GATKRead[] oldSampleReads = readsBySampleIndex[sampleIndex];
        final GATKRead[] newSampleReads = new GATKRead[newSampleReadCount];

        System.arraycopy(oldSampleReads,0,newSampleReads,0,firstDeleted);
        Utils.skimArray(oldSampleReads,firstDeleted, newSampleReads, firstDeleted, removeIndex, firstDeleted);

        // Update the indices for the extant reads from the first deletion onwards.
        for (int r = firstDeleted; r < newSampleReadCount; r++) {
            indexByRead.put(newSampleReads[r], r);
        }

        // Then we skim out the likelihoods of the removed reads.
        final double[][] oldSampleValues = valuesBySampleIndex[sampleIndex];
        final double[][] newSampleValues = new double[alleleCount][newSampleReadCount];
        for (int a = 0; a < alleleCount; a++) {
            System.arraycopy(oldSampleValues[a],0,newSampleValues[a],0,firstDeleted);
            Utils.skimArray(oldSampleValues[a], firstDeleted, newSampleValues[a], firstDeleted, removeIndex, firstDeleted);
        }
        valuesBySampleIndex[sampleIndex] = newSampleValues;
        readsBySampleIndex[sampleIndex] = newSampleReads;
        readListBySampleIndex[sampleIndex] = null; // reset the unmodifiable list.
    }


    private Object2IntMap<GATKRead> readIndexBySampleIndex(final int sampleIndex) {
        if (readIndexBySampleIndex[sampleIndex] == null) {
            final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
            final int sampleReadCount = sampleReads.length;
            readIndexBySampleIndex[sampleIndex] = new Object2IntOpenHashMap<>(sampleReadCount);
            for (int r = 0; r < sampleReadCount; r++) {
                readIndexBySampleIndex[sampleIndex].put(sampleReads[r], r);
            }
        }
        return readIndexBySampleIndex[sampleIndex];
    }


    /**
     * Collect a map stratified per-sample of the base pileups at the provided Location
     * NOTE: Since we shouldn't need to use the pileup if we have more reliable liklihoods, we want to discourage their use
     *
     * @param loc reference location to construct pileups for
     * @return
     */
    public Map<String, List<PileupElement>> getStratifiedPileups(final Locatable loc) {
        return null;
    }

    /**
     * Implements a likelihood matrix per sample given its index.
     */
    private final class SampleMatrix implements LikelihoodMatrix<A> {

        private final int sampleIndex;

        private SampleMatrix(final int sampleIndex) {
            this.sampleIndex = sampleIndex;
        }

        //@Override
        public List<Haplotype> haplotypes() {
            return haplotypes;
        }
        @Override
        public List<GATKRead> reads() {
            return sampleReads(sampleIndex);
        }

        @Override
        public List<A> alleles() {
            return ReadLikelihoods.this.alleles();
        }

        @Override
        public void set(final int alleleIndex, final int readIndex, final double value) {
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            Utils.validIndex(readIndex, valuesBySampleIndex[sampleIndex][alleleIndex].length);
            valuesBySampleIndex[sampleIndex][alleleIndex][readIndex] = value;
        }

        @Override
        public double get(final int alleleIndex, final int readIndex) {
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            Utils.validIndex(readIndex, valuesBySampleIndex[sampleIndex][alleleIndex].length);
            return valuesBySampleIndex[sampleIndex][alleleIndex][readIndex];
        }

        @Override
        public int indexOfAllele(final A allele) {
            Utils.nonNull(allele);
            return ReadLikelihoods.this.indexOfAllele(allele);
        }

        @Override
        public int indexOfRead(final GATKRead read) {
            Utils.nonNull(read);
            return ReadLikelihoods.this.readIndex(sampleIndex, read);
        }

        @Override
        public int numberOfAlleles() {
            return alleles.numberOfAlleles();
        }

        @Override
        public int numberOfReads() {
            return readsBySampleIndex[sampleIndex].length;
        }

        @Override
        public A getAllele(final int alleleIndex) {
            return ReadLikelihoods.this.getAllele(alleleIndex);
        }

        @Override
        public GATKRead getRead(final int readIndex) {
            final GATKRead[] sampleReads = readsBySampleIndex[sampleIndex];
            Utils.validIndex(readIndex, sampleReads.length);
            return sampleReads[readIndex];
        }

        @Override
        public void copyAlleleLikelihoods(final int alleleIndex, final double[] dest, final int offset) {
            Utils.nonNull(dest);
            Utils.validIndex(alleleIndex, valuesBySampleIndex[sampleIndex].length);
            System.arraycopy(valuesBySampleIndex[sampleIndex][alleleIndex], 0, dest, offset, numberOfReads());
        }
    }
}
