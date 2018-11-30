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

public class ZeroStep extends Step {

  public double getStatistic() throws StepException { 
    return 0;
  }

  public double getPvalue() throws StepException {
    return 1.0;
  }

  public double getCenter() { return means_[0]; }

  public HashSet<String> getTags(int numArrayHeader) {
    HashSet<String> set = super.getTags(numArrayHeader);
    set.add("Rest");
    return set;
  }

  public static Step fitStep(GeneData gene, AnalysisMetaData meta) 
    throws StepException {
      Step result = new ZeroStep();
      try {
      Object[] data = gene.getData();
      meta.sanitize(data);
      data = meta.permute(data);
      int count = GeneData.getCount(data, meta.getStart(), meta.getEnd());
      result.num_ = count;
      result.numSteps_ = 0;
      double mean = GeneData.getMean(data, meta.getStart(), meta.getEnd());
      double sstot = GeneData.getSquareError(data, meta.getStart(), meta.getEnd());
      result.sstot_ = sstot;
      result.sse_ = sstot;
      result.means_ = new double[1];
      result.means_[0] = mean;
      }
      catch(Exception e) {
        // Bad count
        result.numSteps_ = -1;
      }
      return result;
    }

};

