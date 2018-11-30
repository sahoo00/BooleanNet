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

public class BestStep extends Step {
  Step[] fittedSteps_;
  AnalysisMetaData meta_;

  public BestStep() {
    super();
    fittedSteps_ = null;
    meta_ = null;
  }

  public double getStatistic() throws StepException { 
    return 0;
  }

  public double getPvalue() throws StepException {
    return 1.0;
  }

  public Step findStep(int index) {
    return fittedSteps_[index];
  }

  public double getF12() throws StepException {
    Step step1 = findStep(1); // get oneStep
    Step step2 = findStep(2); // get twoStep
    double sse1 = step1.getSse();
    double sse2 = step2.getSse();
    double f12;
    int num = step1.getNum();
    if (num > 5) {
      f12 = (sse1 - sse2)/(sse2/(num - 5));
    }
    else {
      f12 = (sse1 - sse2)/sse2;
    }
    return f12;
  }

  public double getP12() throws StepException {
    Step step1 = findStep(1); // get oneStep
    Step step2 = findStep(2); // get twoStep
    double sse1 = step1.getSse();
    double sse2 = step2.getSse();
    double p12;
    double f12;
    int num = step1.getNum();
    if (num > 5) {
      f12 = (sse1 - sse2)/(sse2/(num - 5));
      p12 = 1 - Utils.Fisher(f12, 1, num - 5);
    }
    else {
      f12 = (sse1 - sse2)/sse2;
      p12 = 1 - Utils.Fisher(f12, 1, 1);
    }
    //step1.print();
    //step2.print();
    //System.out.println("f12 = " + f12 + " p12 = " + p12);
    return p12;
  }

  public Step findBestSingleStep() throws StepException {
    if (fittedSteps_ == null) {
      return null;
    }
    // StepMiner Algorithm
    double thr = meta_.getPvalueThr();
    double p1 = fittedSteps_[1].getPvalue();

    if (p1 < thr) {
      return fittedSteps_[1];
    }
    return fittedSteps_[0];
  }

  public Step findBestBothStep() throws StepException {
    if (fittedSteps_ == null) {
      return null;
    }
    // StepMiner Algorithm
    double thr = meta_.getPvalueThr();
    double p1 = fittedSteps_[1].getPvalue();
    double p2 = fittedSteps_[2].getPvalue();
    double p12 = getP12();

    if (p1 < thr && p12 > thr) {
      return fittedSteps_[1];
    }
    if (p2 < thr) {
      return fittedSteps_[2];
    }
    return fittedSteps_[0];
  }

  public Step findSelectTwoStep() throws StepException {
    if (fittedSteps_ == null) {
      return null;
    }
    // StepMiner Algorithm
    double thr = meta_.getPvalueThr();
    double p1 = fittedSteps_[1].getPvalue();
    double p2 = fittedSteps_[2].getPvalue();
    double p12 = getP12();

    if (p1 < thr && p12 > thr) {
      return fittedSteps_[0];
    }
    if (p2 < thr) {
      return fittedSteps_[2];
    }
    return fittedSteps_[0];
  }

  public Step findBestTwoStep() throws StepException {
    if (fittedSteps_ == null) {
      return null;
    }
    double thr = meta_.getPvalueThr();
    double p2 = fittedSteps_[2].getPvalue();

    if (p2 < thr) {
      return fittedSteps_[2];
    }
    return fittedSteps_[0];
  }

  public Step findBestStep(int type) throws StepException {
    if (type == 0) {
      return findBestSingleStep();
    }
    else if (type == 1) {
      return findBestBothStep();
    }
    else if (type == 2) {
      return findBestTwoStep();
    }
    else if (type == 3) {
      return findSelectTwoStep();
    }
    throw new StepException("Undefined analysis type : " + type);
  }

  public static Step fitStep(GeneData gene, AnalysisMetaData meta) 
    throws StepException {
      BestStep res = new BestStep();
      res.fittedSteps_ = new Step[3];
      res.fittedSteps_[0] = ZeroStep.fitStep(gene, meta);
      res.fittedSteps_[1] = OneStep.fitStep(gene, meta);
      res.fittedSteps_[2] = TwoStep.fitStep(gene, meta);
      res.meta_ = meta;

      return res;
    }

};

