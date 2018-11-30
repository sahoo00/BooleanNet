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

package tools.microarray;

import java.util.*;

/*
 *  Inspired by ItrLLSimpute method 
 *      http://www.cs.ualberta.ca/~ghlin/src/WebTools/Imputation/ItrLLSimpute.m
 * References:
 *
 * 1. Z. Cai, M. Heydari, and G. Lin.
 *    Iterated Local Least Squares Imputation for Microarray Missing Values.
 *    Journal of Bioinformatics and Computational Biology.
 *    Accepted for Publication on June 1, 2006.
 * 2. Z. Cai, M. Heydari, and G. Lin.
 *    Microarray Missing Value Imputation by Iterated Local Least Squares.
 *    In Proceedings of the Fourth Asia-Pacific Bioinformatics Conference 
 *    (APBC 2006).
 *    Taipei, Taiwan. February 13-16, 2006.
 *    Pages 159-168, 2006. 
 *
 */

public class Impute {

  Data data_;
  String algo_;
  int start_;
  int end_;
  Vector<Integer> geneIndex_; // Index into data_
  Vector<Integer> missingGeneIndex_;

  Random random_;
  double sumSqErr_;
  double sumSqAns_;
  double sumAns_;
  int num_;
  byte[] bytes_;
  HashMap<Integer, Double> answer_;

  public Impute(Data data, String algo) {
    data_ = data;
    algo_ = algo;
    start_ = data_.getNumArrayHeader();
    end_ = data_.getNumColumns() - 1;
    geneIndex_ = new Vector<Integer>();
    missingGeneIndex_ = new Vector<Integer>();
    random_ = null;
    sumSqErr_ = 0;
    sumSqAns_ = 0;
    sumAns_ = 0;
    num_ = 0;
    bytes_ = null;
  }

  public Data getData() { return data_; }

  void populateGenes() {
    int numMissing = 0;
    for (int i = data_.getNumGeneHeader(); i < data_.getNumRows(); i++) {
      GeneData gene = data_.getGeneData(i);
      if (gene.getMissingPoints(start_, end_) <= numMissing) {
        geneIndex_.add(new Integer(i));
      }
      else {
        missingGeneIndex_.add(new Integer(i));
      }
    }
    System.out.println("Number of genes with missing points: " +
        missingGeneIndex_.size() );
    System.out.println("Number of genes with no missing points: " +
        geneIndex_.size() );
  }

  GeneData impute_rowavg_gene(GeneData gene) {
    Object[] data = gene.getData();
    double center;
    try {
      center = GeneData.getMean(data, start_, end_);
    }
    catch(Exception e) {
      center = 0.0;
    }
    for (int j =start_; j <= end_; j++) {
      if (data[j] == null) {
        data[j] = new Double(center);
      }
    }
    return new GeneData(data);
  }

  void impute_rowavg() {
    for (int i = data_.getNumGeneHeader(); i < data_.getNumRows(); i++) {
      GeneData gene = data_.getGeneData(i);
      GeneData res = impute_rowavg_gene(gene);
      data_.setGeneData(res, i);
    }
  }

  /**
   * Finds Distance after replacing all the null value in b as row average
   * @param a - genes of interest
   * @param b - Other genes in the array that might be closer
   * @returns The normalized euclidean distance between two genes
   */
  static double findDistance(GeneData a, GeneData b, int start, int end) {
    GeneData cl = (GeneData) b.clone();
    Object[] data = cl.getData();
    double center;
    try {
      center = GeneData.getMean(data, start, end);
    }
    catch(Exception e) {
      center = 0.0;
    }
    for (int j =start; j <= end; j++) {
      if (data[j] == null) {
        data[j] = new Double(center);
      }
    }
    return GeneData.findDistance(a, cl, start, end);
  }

  Vector<GeneData> findKneighbors(GeneData gene, int k) {
    Vector<GeneData> res = new Vector<GeneData>();
    Vector<Integer> index = new Vector<Integer>();
    Vector<Double> distance = new Vector<Double>();
    for (int i = data_.getNumGeneHeader(); i < data_.getNumRows(); i++) {
      GeneData b = data_.getGeneData(i);
      if (gene != b) {
        index.add(new Integer(i));
        double dist = findDistance(gene, b, start_, end_);
        distance.add(new Double(-dist));
      }
    }
    index = ArrayOrder.sortCorrelation(index, distance);
    int mink = k;
    if (mink > index.size()) {
      mink = index.size();
    }
    //System.out.print("GeneIndex : ");
    for (int i = 0; i < mink ; i++) {
      Integer ind = index.get(i);
      GeneData a = data_.getGeneData(ind.intValue());
      double dist = findDistance(gene, a, start_, end_);
      // System.out.print(ind.intValue() + ":" +
      //     Matrix.formatString("0.##", dist) + " ");
      res.add(a);
    }
    //System.out.println();
    return res;
  }

  Matrix getW(GeneData gene) {
    Object[] data = gene.getData();
    int count = GeneData.getCount(data, start_, end_);
    Matrix m = new Matrix(count, 1);
    int index = 0;
    for (int i = start_; i <= end_ ; i++) {
      if (data[i] != null) {
        Double val = (Double) data[i];
        //System.out.println("W:" + count + " : " + i + " : "+ index);
        m.setCell(index, 0, val.doubleValue());
        index++;
      }
    }
    return m;
  }

  Matrix getA(GeneData gene, Vector<GeneData> list) {
    Object[] data = gene.getData();
    int count = GeneData.getCount(data, start_, end_);
    Matrix m = new Matrix(list.size(), count);
    for (int i =0; i < list.size(); i++) {
      int index = 0;
      GeneData g = list.get(i);
      g = (GeneData) g.clone();
      g = impute_rowavg_gene(g);
      for (int j = start_; j <= end_ ; j++) {
        if (data[j] != null) {
          Double val = (Double) g.getDataAt(j);
          //System.out.println("A:" + count + " : " + j + " : "+ index);
          m.setCell(i, index, val.doubleValue());
          index++;
        }
      }
    }
    return m;
  }

  Matrix getB(GeneData gene, Vector<GeneData> list) {
    Object[] data = gene.getData();
    int count = GeneData.getCount(data, start_, end_);
    count = end_ - start_ + 1 - count;
    Matrix m = new Matrix(list.size(), count);
    for (int i =0; i < list.size(); i++) {
      int index = 0;
      GeneData g = list.get(i);
      g = (GeneData) g.clone();
      g = impute_rowavg_gene(g);
      for (int j = start_; j <= end_ ; j++) {
        if (data[j] == null) {
          Double val = (Double) g.getDataAt(j);
          //System.out.println("B:" + count + " : " + j + " : "+ index);
          m.setCell(i, index, val.doubleValue());
          index++;
        }
      }
    }
    return m;
  }

  Matrix getDist(GeneData gene, Vector<GeneData> list) {
    Matrix m = new Matrix(list.size(), 1);
    for (int i =0; i < list.size(); i++) {
      GeneData g = list.get(i);
      double dist = findDistance(gene, g, start_, end_);
      m.setCell(i, 0, dist);
    }
    return m;
  }

  void generateArtificialMissingPoints(GeneData gene) {
    // Generating Artificial Missing Points
    bytes_ = new byte[end_ - start_ + 1];
    random_.nextBytes(bytes_);
    answer_ = new HashMap<Integer, Double>();
    Object[] data = gene.getData();
    //System.out.print("Index : ");
    for (int j =0; j < bytes_.length; j++) {
      if ( (bytes_[j]&0xff) < 3) {
        Double val = (Double) data[j + start_];
        if (val != null) {
          //System.out.print(j + " ");
          answer_.put(new Integer(j), val);
          data[j + start_] = null;
        }
        else {
          bytes_[j] = (byte) 4;
        }
      }
    }
    //System.out.println();
  }

  void evaluatePerformance(GeneData gene) {
    Object[] data = gene.getData();
    // Evaluate Missing points answer
    for (int j =0; j < bytes_.length; j++) {
      if ( (bytes_[j]&0xff) < 3) {
        Double val = (Double) data[j + start_];
        Double ans = answer_.get(new Integer(j));
        //System.out.println("Actual = " + ans + " Computed = " + val);
        double err = val.doubleValue() - ans.doubleValue();
        sumSqErr_ += err * err;
        sumSqAns_ += ans.doubleValue() * ans.doubleValue();
        sumAns_ += ans.doubleValue();
        num_++;
        data[j + start_] = ans;
      }
    }
  }

  void impute_gene(int index, GeneData gene, int k) {
    Object[] data = gene.getData();
    int miss = gene.getMissingPoints(start_, end_);
    int maxMiss = end_ - start_ + 1;
    //System.out.println("M ("+index+"):" + miss);
    if (miss != 0) {
      if (miss == maxMiss) {
        for (int j = start_; j <= end_ ; j++) {
          if (data[j] == null) {
            gene.setDataAt(j, new Double(0));
          }
        }
      }
      else {
        Matrix U = null;
        Vector<GeneData> neighbors = findKneighbors(gene, k);
        if (algo_.equals("LLSimpute")) {
          // LLSimpute Algorithm
          Matrix W = getW(gene);
          //System.out.println("W:\n" + W.transpose().toString());
          Matrix A = getA(gene, neighbors);
          //System.out.println("A:\n" + A.toString());
          Matrix B = getB(gene, neighbors);
          //System.out.println("B:\n" + B.toString());
          Matrix AT = A.transpose();
          Matrix X = AT.pseudoInverse().multiply(W);
          //System.out.println("X:\n" + X.transpose().toString());
          U = B.transpose().multiply(X);
          //System.out.println("U:\n" + U.transpose().toString());
          //System.out.println("K:\n" + AT.multiply(X).transpose().toString());
        }
        else if (algo_.equals("KNNimpute")) {
          //System.out.println("gene:\n" + data[1]);
          Matrix B = getB(gene, neighbors);
          //System.out.println("B:\n" + B.toString());
          Matrix dist = getDist(gene, neighbors);
          //System.out.println("dist:\n" + dist.transpose().toString());
          U = B.weightedAverage(dist);
          //System.out.println("U:\n" + U.transpose().toString());
        }
        else {
          System.out.println("Please provide an algorithm from LLSimpute,KNNimpute");
          System.exit(0);
        }
        int index1 = 0;
        for (int j = start_; j <= end_ ; j++) {
          if (data[j] == null) {
            double val = U.getCell(index1, 0);
            gene.setDataAt(j, new Double(val));
            index1++;
          }
        }
      }
    }
  }

  void printEvaluation() {
    // Compute Performance
    double stdAns = Math.sqrt( sumSqAns_/num_ - sumAns_ * sumAns_/num_/num_);
    double rms = Math.sqrt(sumSqErr_ / num_);
    double nrmse = rms / stdAns;
    //System.out.println(" Err = " + sumSqErr_ + " SqAns = " + sumSqAns_ + "  Ans = " + sumAns_);
    //System.out.println(" Rms = " + rms + " stdAns = " + stdAns);
    System.out.println(" Num = " + num_ + " Nrmse = " + nrmse);
  }

  public void impute(int k) {
    populateGenes();
    random_ = new Random(100);
    sumSqErr_ = 0;
    sumSqAns_ = 0;
    sumAns_ = 0;
    num_ = 0;
    for (int i =0; i < geneIndex_.size(); i++) {
      Integer ind = geneIndex_.get(i);
      GeneData gene = data_.getGeneData(ind.intValue());
      generateArtificialMissingPoints(gene);
      impute_gene(ind.intValue(), gene, k);
      evaluatePerformance(gene);
    }
    for (int i =0; i < missingGeneIndex_.size(); i++) {
      Integer ind = missingGeneIndex_.get(i);
      GeneData gene = data_.getGeneData(ind.intValue());
      impute_gene(ind.intValue(), gene, k);
    }
    printEvaluation();
  }

  public static void main(String[] args) throws Exception {
    Data data = tools.microarray.FileReader.PCLFileReader.readFile(args[0]);
    data.convertDoubles();
    //Impute im = new Impute(data, "LLSimpute");
    Impute im = new Impute(data, "KNNimpute");
    im.impute(10);
    tools.microarray.FileWriter.PCLFileWriter.writeFile(data, "label.pcl", null);
  }

}
