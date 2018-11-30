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

public class TwoStep extends Step {

  public TwoStep() {
    super();
  }

  public double getStatistic() throws StepException { 
    if (num_ < 0) {
        return 0.0;
    }
    // Degrees of freedom = 5
    if (num_ > 5) {
      return (sstot_ - sse_)/4/(sse_/(num_ - 5)); 
    }
    else {
      return (sstot_ - sse_)/3/sse_; 
    }
  }

  public double getPvalue() throws StepException {
    if (num_ < 0) {
        return 1.0;
    }
    double f = getStatistic();
    double p;
    if (num_ > 5) {
      p = 1 - Utils.Fisher(f, 4, num_ - 5);
    }
    else {
      p = 1 - Utils.Fisher(f, 3, 1);
    }
    //System.out.println("f2 = " + f + " p2 = " + p);
    return p;
  }

  public double getCenter() { return (means_[0] + means_[1])/2; }

  public HashSet<String> getTags(int numArrayHeader) {
    HashSet<String> set = super.getTags(numArrayHeader);
    set.add("Step");
    if (means_[0] <  means_[1]) {
      set.add("UpDown");
      set.add("Up-" + (steps_[0]-numArrayHeader));
      set.add("Up-" + (steps_[0]-numArrayHeader) + "-Down");
      set.add("Down-" + (steps_[1]-numArrayHeader));
      set.add("Up-Down-" + (steps_[1]-numArrayHeader));
      set.add("Up-" + (steps_[0]-numArrayHeader) + "-Down-" + (steps_[1]-numArrayHeader));
    }
    else {
      set.add("DownUp");
      set.add("Down-" + (steps_[0]-numArrayHeader));
      set.add("Down-" + (steps_[0]-numArrayHeader) + "-Up");
      set.add("Up-" + (steps_[1]-numArrayHeader));
      set.add("Down-Up-" + (steps_[1]-numArrayHeader));
      set.add("Down-" + (steps_[0]-numArrayHeader) + "-Up-" + (steps_[1]-numArrayHeader));
    }
    return set;
  }

  public static Step fitStep(GeneData gene, AnalysisMetaData meta) 
    throws StepException {
      //meta.print();
      Step result = new TwoStep();
      try {
      double thr = meta.getPvalueThr();
      Object[] data = gene.getData();
      meta.sanitize(data);
      data = meta.permute(data);
      data = meta.clipData(data);
      data = meta.concatData(data, data);
      GeneData stepdata = new GeneData(data);
      int num = meta.getNum();
      int[] stepSearch = meta.getStepSearch();
      int[] stepSearch1 = meta.getStepSearch1();
      int count = GeneData.getCount(data, 0, num-1);
      result.num_ = count;
      double mean = GeneData.getMean(data, 0, num-1);
      double sstot = GeneData.getSquareError(data, 0, num-1);
      result.sstot_ = sstot;
      result.numSteps_ = 2;
      result.steps_ = new int[2];
      result.steps_[0] = -1;
      result.steps_[0] = num - 1;
      result.means_ = new double[3];
      result.means_[0] = 0;
      result.means_[1] = mean;
      result.means_[2] = 0;
      result.sse_ = sstot;
      for (int i = 0; i < stepSearch.length; i++) {
        int index = stepSearch[i] - meta.getStart();
        Double entry = (Double) data[index];
        
        AnalysisMetaData meta1 = new AnalysisMetaData(index+1, index+num, null, stepSearch1, thr);
        //meta1.print();
        //System.out.println(index + " " + num);
        meta1.intersectStepSearch(index+2, num-2);
        //meta1.print();
        Step onestep = OneStep.fitStep(stepdata, meta1);
        if (onestep.getNumSteps() == 1 && onestep.getSse() < result.sse_) {
            result.sse_ = onestep.getSse();
            result.steps_[0] = index;
            result.steps_[1] = onestep.getStep(0);
            result.means_[0] = onestep.getMean(1);
            result.means_[1] = onestep.getMean(0);
            result.means_[2] = result.means_[0];
        }
      }

      result.steps_[0] += meta.getStart();
      result.steps_[1] += meta.getStart();
      double m1 = result.means_[0];
      double m2 = result.means_[1];
      if (m1 < m2) {
        result.label_ = 3;
      }
      else {
        result.label_ = 4;
      }
      }
      catch(Exception e) {
        // Exceptions raised during fitting steps
        result.numSteps_ = -1;
      }
      //System.exit(0);
      //result.print();
      return result;
    }

};

