/*

Copyright (c) 2006, the Board of Trustees of Leland
Stanford Junior University.

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the
distribution.

    * Neither the name of Stanford University nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

/*
 Author: Debashis Sahoo <sahoo@stanford.edu>
 */

package tools.microarray.StepMiner;

import tools.microarray.GeneData;
import java.util.HashSet;

public class OneStep extends Step {

  public OneStep() {
    super();
  }

  public double getStatistic() throws StepException { 
    if (num_ < 0) {
        return 0.0;
    }
    // Degrees of freedom = 4
    if (num_ > 4) {
      return (sstot_ - sse_)/3/(sse_/(num_ - 4)); 
    }
    else {
      return (sstot_ - sse_)/2/sse_; 
    }
  }

  public double getPvalue() throws StepException {
    if (num_ < 0) {
        return 1.0;
    }
    double f = getStatistic();
    double p;
    if (num_ > 4) {
      p = 1 - Utils.Fisher(f, 3, num_ - 4);
    }
    else {
      p = 1 - Utils.Fisher(f, 2, 1);
    }
    // System.out.println("f1 = " + f + " p1 = " + p);
    return p;
  }

  public double getCenter() { return (means_[0] + means_[1])/2; }

  public HashSet<String> getTags(int numArrayHeader) {
    HashSet<String> set = super.getTags(numArrayHeader);
    set.add("Step");
    if (means_[0] <  means_[1]) {
      set.add("Up");
      set.add("Up-" + (steps_[0]-numArrayHeader));
      set.add("UpOnly-" + (steps_[0]-numArrayHeader));
    }
    else {
      set.add("Down");
      set.add("Down-" + (steps_[0]-numArrayHeader));
      set.add("DownOnly-" + (steps_[0]-numArrayHeader));
    }
    return set;
  }

  public static Step fitStep(GeneData gene, AnalysisMetaData meta) 
    throws StepException {
      Step result = new OneStep();
      try {
      Object[] data = gene.getData();
      meta.sanitize(data);
      data = meta.permute(data);
      double[] sseArray = new double[meta.getNum()];
      for (int i = 0; i < sseArray.length; i++) {
        sseArray[i] = 0.0;
      }
      double sum = GeneData.getSum(data, meta.getStart(), meta.getEnd());
      int count = GeneData.getCount(data, meta.getStart(), meta.getEnd());
      result.num_ = count;
      result.numSteps_ = 1;
      double mean = GeneData.getMean(data, meta.getStart(), meta.getEnd());
      double sstot = GeneData.getSquareError(data, meta.getStart(), meta.getEnd());
      result.sstot_ = sstot;
      double sum1 = 0.0;
      int count1 = 0;
      double m1 = 0.0;
      double sum2 = sum;
      int count2 = count;
      double m2 = (sum/count);
      double sum1sq = 0.0;
      double sum2sq = GeneData.getSquareError(data, meta.getStart(), meta.getEnd());
      double sse = sum1sq + sum2sq;

      for (int i = 0; i < sseArray.length; i++) {
        Double entryVal = (Double) data[i + meta.getStart()];
        if (entryVal == null) {
          sseArray[i] = sse;
          continue;
        }
        double entry = entryVal.doubleValue();
        count1 ++;
        count2 --;
        if (count2 == 0) {
          sseArray[i] = sstot;
          continue;
        }
        double tmp = (mean - (entry + sum1)/count1);
        sum1sq = sum1sq + (entry-mean) * (entry-mean) - tmp * tmp * count1
            + (count1 - 1) * (mean - m1) * (mean - m1);
        tmp = (mean - (sum2 - entry)/count2);
        sum2sq = sum2sq - (entry-mean) * (entry-mean) - tmp * tmp * count2
            + (count2 + 1) * (mean - m2) * (mean - m2);
        sum1 += entry;
        sum2 -= entry;
        m1 = sum1/count1;
        m2 = sum2/count2;
        sse = sum1sq + sum2sq;
        sseArray[i] = sse;
      }

      /*
      for (int i = 0; i < sseArray.length; i++) {
        System.out.println(sseArray[i]);
      }
      System.out.println("------------");
      */

      double bestSse = Double.MAX_VALUE;
      int bestIndex = 0;
      int[] stepSearch = meta.getStepSearch();
      for (int i = 0; i < stepSearch.length; i++) {
        int index = stepSearch[i];
        /* 
         * Debug code
        Double entry = (Double) data[index];
        sum1sq = GeneData.getSquareError(data, meta.getStart(), index);
        sum2sq = GeneData.getSquareError(data, index + 1, meta.getEnd());
        sse = sum1sq + sum2sq;
        System.out.println(sseArray[index-meta.getStart()] + "\t" + entry +
            "\t" + sse + "\t" + sum1sq);
        */
        if (sseArray[index-meta.getStart()] < bestSse) {
            bestSse = sseArray[index-meta.getStart()];
            bestIndex = index;
        }
      }
      m1 = GeneData.getMean(data, meta.getStart(), bestIndex);
      m2 = GeneData.getMean(data, bestIndex + 1, meta.getEnd());
      /*
      sum1sq = GeneData.getSquareError(data, meta.getStart(), bestIndex);
      sum2sq = GeneData.getSquareError(data, bestIndex + 1, meta.getEnd());
      sse = sum1sq + sum2sq;
      if (sse != bestSse) {
        throw new StepException("SSE calculation is wrong :" + sse + ":"+bestSse);
      }
      */
      result.sse_ = bestSse;
      result.steps_ = new int[1];
      result.steps_[0] = bestIndex;
      result.means_ = new double[2];
      result.means_[0] = m1;
      result.means_[1] = m2;
      if (m1 < m2) {
        result.label_ = 1;
      }
      else {
        result.label_ = 2;
      }
      }
      catch(Exception e) {
        // Exception in fitting step.
        // e.printStackTrace();
        result.numSteps_ = -1;
      }
      return result;
    }

};

