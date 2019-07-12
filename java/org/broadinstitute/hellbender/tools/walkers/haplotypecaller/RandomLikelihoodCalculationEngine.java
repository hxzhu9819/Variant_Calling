package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.util.List;
import java.util.Map;
import java.util.Random;

/**
 * Random likelihoods generator, used for testing/benchmarking purposes.
 */
public final class RandomLikelihoodCalculationEngine implements ReadLikelihoodCalculationEngine {
    @Override
    public List<ReadLikelihoods<Haplotype>> computeReadLikelihoods(AssemblyResultSet assemblyResultSet, SampleList samples,
                                                             Map<String, List<GATKRead>> perSampleReadList){
        return null;
    }

    public List<ReadLikelihoods<Haplotype>> hardware_compute(final AssemblyResultSet assemblyResultSet,
                                                             final SampleList samples,
                                                             final Map<String, List<GATKRead>> perSampleReadList,
                                                             final List<Integer> log2InitialValues,
                                                             final List<Float> realInitialValues){
        return null;
    }

    // added by Chenhao: just for pass compile
    public void printInitialValues(final byte[] haplotypeBases){
        return;
    }

    // added by Chenhao: just for pass compile
    public int getInitialUpper(final byte[] haplotypeBases){ return 0; }

    // added by Chenhao: just for pass compile
    public float getInitialLower(final byte[] haplotypeBases){ return 0; }
                                                             
    //@Override
    //public ReadLikelihoods<Haplotype> computeReadLikelihoods(final AssemblyResultSet assemblyResultSet,
    //                                              final SampleList samples,
    //                                              final Map<String, List<GATKRead>> reads) {
    //    Utils.nonNull(assemblyResultSet, "assemblyResultSet is null");
    //    Utils.nonNull(samples, "samples is null");
    //    Utils.nonNull(reads, "perSampleReadList is null");
    //    final AlleleList<Haplotype> haplotypes = new IndexedAlleleList<>(assemblyResultSet.getHaplotypeList());
    //    final ReadLikelihoods<Haplotype> result = new ReadLikelihoods<>(samples, haplotypes, reads);
    //    final Random rnd = Utils.getRandomGenerator();
    //    final int sampleCount = samples.numberOfSamples();
    //    final int alleleCount = haplotypes.numberOfAlleles();
    //    for (int i = 0; i < sampleCount; i++)  {
    //        final List<GATKRead> sampleReads = result.sampleReads(i);
    //        final int readCount = sampleReads.size();
    //        final LikelihoodMatrix<Haplotype> sampleLikelihoods = result.sampleMatrix(i);
    //        for (int a = 0; a < alleleCount; a++) {
    //            for (int r = 0; r < readCount; r++) {
    //                sampleLikelihoods.set(a, r, -Math.abs(rnd.nextDouble()));
    //            }
    //        }
    //    }
    //    return result;
    //}

    @Override
    public void close() {
    }

}
