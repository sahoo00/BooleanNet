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

import java.util.HashSet;

public class AnalysisMetaData {

  int start_;
  int end_;
  int[] perm_;
  int[] stepSearch_; // Allowed steps for the first step(Column indices of PCL)
                     // GeneData may contain id and name as first columns
                     // In this case stepSearch would be the actual indices
                     // of the GeneData.
  int[] stepSearch1_;  // Allowed steps for the second step
                       // Column indices - start_
                       // Index to the data part of the GeneData
  double pvalueThr_;

  public AnalysisMetaData() {
    start_ = -1;
    end_ = -1;
    perm_ = null;
    stepSearch_ = null;
    stepSearch1_ = null;
    pvalueThr_ = 0.05;
  }

  public AnalysisMetaData(int s, int e, int[] p, int[] st, double thr) {
    start_ = s;
    end_ = e;
    perm_ = p;
    stepSearch_ = st;
    stepSearch1_ = null;
    pvalueThr_ = thr;
  }

  public int getStart() { return start_; }
  public int getEnd() { return end_; }
  public int getNum() { return end_ - start_ + 1; }
  public double getPvalueThr() { return pvalueThr_;}
  public int[] getPermutation() { return perm_;}
  public void setStepSearch(int[] st) { stepSearch_ = st; }
  public void setStepSearch1(int[] st) { stepSearch1_ = st; }


  // Search from index s to index e for steps
  public void stepSearch(int s, int e) {
    stepSearch_ = new int[e - s + 1];
    for (int i=0; i < stepSearch_.length; i++) {
      stepSearch_[i] = i + s;
    }
  }

  // Search from index s to index e for steps
  public void stepSearch1(int s, int e) {
    stepSearch1_ = new int[e - s + 1];
    for (int i=0; i < stepSearch1_.length; i++) {
      stepSearch1_[i] = i + s;
    }
  }

  public void intersectStepSearch(int s, int e) throws StepException {
    stepSearch_ = getStepSearch();
    HashSet<Integer> a = new HashSet<Integer>();
    for (int i =0; i < stepSearch_.length; i++) {
      a.add(new Integer(stepSearch_[i]));
    }
    int count = 0;
    for (int i=s; i <= e ; i++) {
      if (a.contains(new Integer(i))) {
        count ++;
      }
    }
    int[] st = new int[count];
    int index =0;
    for (int i=s; i <= e ; i++) {
      if (a.contains(new Integer(i))) {
        st[index++] = i;
      }
    }
    stepSearch_ = st;
  }

  public void intersectStepSearch1(int s, int e) throws StepException {
    stepSearch1_ = getStepSearch1();
    HashSet<Integer> a = new HashSet<Integer>();
    for (int i =0; i < stepSearch1_.length; i++) {
      a.add(new Integer(stepSearch1_[i]));
    }
    int count = 0;
    for (int i=s; i <= e ; i++) {
      if (a.contains(new Integer(i))) {
        count ++;
      }
    }
    int[] st = new int[count];
    int index =0;
    for (int i=s; i <= e ; i++) {
      if (a.contains(new Integer(i))) {
        st[index++] = i;
      }
    }
    stepSearch1_ = st;
  }

  public void sanitize(Object[] data) throws StepException {
    if (data == null) {
      throw new StepException("GeneData is null");
    }
    if (start_ < 0) {
      throw new StepException("start_ is negative");
    }
    if (start_ > end_) {
      throw new StepException("start_ > end_");
    }
    if (end_ >= data.length) {
      throw new StepException("end_ >= data.length");
    }
    if (start_ >= data.length) {
      throw new StepException("start_ >= data.length");
    }
  }

  public int[] getStepSearch() {
    if (stepSearch_ == null) {
      stepSearch(start_, end_);
    }
    return stepSearch_;
  }

  public int[] getStepSearch1() {
    if (stepSearch1_ == null) {
      stepSearch1(0, end_ - start_);
    }
    return stepSearch1_;
  }

  public Object[] permute(Object[] data) {
    if (perm_ == null) {
      return data;
    }
    else {
      Object[] res = new Object[data.length];
      for (int i = 0; i < data.length; i++) {
        int entry = i;
        if (i >= start_ && (i-start_) < perm_.length) {
          entry = perm_[i-start_] + start_;
        }
        res[i] = data[entry];
      }
      return res;
    }
  }

  public Object[] clipData(Object[] data) {
    Object[] res = new Object[end_-start_+1];
    for (int i = start_; i <= end_; i++) {
        res[i-start_] = data[i];
    }
    return res;
  }

  public Object[] concatData(Object[] data1, Object[] data2) {
    Object[] res = new Object[data1.length+data2.length];
    for (int i = 0; i < data1.length; i++) {
        res[i] = data1[i];
    }
    for (int i = 0; i < data2.length; i++) {
        res[i+data1.length] = data2[i];
    }
    return res;
  }

  public String toString() {
    String res = "Start:" + start_ + " End:" + end_ + " Pvalue:" + pvalueThr_ ;
    res += "\n";
    if (perm_ != null) {
      res += "Perm: ";
      for (int i =0; i < perm_.length; i++) {
        res += perm_[i] + " ";
      }
      res += "\n";
    }
    getStepSearch();
    if (stepSearch_ != null) {
      res += "stepSearch: ";
      for (int i =0; i < stepSearch_.length; i++) {
        res += stepSearch_[i] + " ";
      }
      res += "\n";
    }
    return res;
  }

  public void print() {
    System.out.print(toString());
  }

};

