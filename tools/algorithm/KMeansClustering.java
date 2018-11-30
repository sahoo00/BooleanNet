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

package tools.algorithm;

import tools.microarray.Matrix;

import java.util.*;

public class KMeansClustering {
  double learningRate = 0.5;

  Matrix data_;
  double centroids_[][];    // Centroids of the clusters
  int numElements_[];       // Number of elements in each clusters
  int k_;

  public KMeansClustering(int k, Matrix d) {
    k_ = k;
    data_ = d;
    centroids_ = new double[k_][data_.getNumCols()];
    for (int i =0; i < k_; i++) {
      for (int j=0; j < data_.getNumCols(); j++) {
        centroids_[i][j] = data_.getCell(i,j);
      }
    }
    numElements_ = new int[k_];
    resetElements();
  }

  public double[][] getCentroids() { return centroids_; }

  void resetElements() {
    for (int i =0; i < k_; i++) {
      numElements_[i] = 0;
    }
  }

  public Matrix[] getClusters() {
    Matrix[] res = new Matrix[k_];
    // Compute NumElements
    resetElements();
    for (int i=0; i < data_.getNumRows(); i++) {
      double[] c = getClosestCentroid(data_.getRow(i));
    }
    // create clusters
    int[] indices = new int[k_];
    for (int i =0; i < res.length; i++) {
      res[i] = new Matrix(numElements_[i], data_.getNumCols());
      // System.out.println(" Num = " + numElements_[i] + "x" + data_.getNumCols());
      indices[i] = 0;
    }
    for (int i=0; i < data_.getNumRows(); i++) {
      int k = getClosestCentroidIndex(data_.getRow(i));
      for (int j =0; j < data_.getNumCols(); j++) {
        // System.out.println(" k = " + k + " Index = " + indices[k] + "->" + j);
        res[k].setCell(indices[k], j, data_.getCell(i,j));
      }
      indices[k]++;
    }
    return res;
  }

  void moveTowards(double[] centroid, double[] datum) {
    for (int i = 0; i < datum.length; i++) {
      centroid[i] += ((datum[i] - centroid[i]) * learningRate); 
    }
  }

  double getDistance(double[] datum, double[] centroid) {
    double d = 0.0;
    for (int i = 0; i < datum.length; i++) {
      d += Math.pow(datum[i] - centroid[i], 2); 
    }
    return(Math.sqrt(d));
  }

  int getClosestCentroidIndex(double[] datum) {
    double min = Double.MAX_VALUE;
    int k = -1;
    for (int i = 0; i < centroids_.length; i++) {
      double d = getDistance(datum, centroids_[i]);
      if (d < min) { k = i; min = d; } 
    }
    return k;
  }

  // Call this once for each row of the matrix
  double[] getClosestCentroid(double[] datum) {
    int k = getClosestCentroidIndex(datum);
    numElements_[k] ++;
    return centroids_[k];
  }

  void printDatum(double[] datum) {
    Vector<Double> v = new Vector<Double>();
    for (int j = 0; j < datum.length; j++) {
      v.add(new Double(datum[j]));
    }
    System.out.println(v);
  }

  void printCentroids() {
    for (int i = 0; i < centroids_.length; i++) {
        printDatum(centroids_[i]);
    }
    System.out.println("-------------------");
  }

  double[][] cloneCentroids() {
    double[][] res = new double[k_][data_.getNumCols()];
    for (int i =0; i < k_; i++) {
      for (int j=0; j < data_.getNumCols(); j++) {
        res[i][j] = centroids_[i][j];
      }
    }
    return res;
  }

  boolean checkCentoids(double[][] last, double[][] current) {
    double max = 0;
    for (int i =0; i < k_; i++) {
      for (int j=0; j < data_.getNumCols(); j++) {
        double res = Math.abs(last[i][j] - current[i][j]);
        if (res > max) {
          max = res;
        }
      }
    }
    return (max < 0.0001);
  }

  public void computeClusters() {
    printCentroids();
    while (true) {
      double[][] last = cloneCentroids();
      resetElements();
      for (int i = 0; i < data_.getNumRows(); i++) {
        moveTowards(getClosestCentroid(data_.getRow(i)), data_.getRow(i)); 
      }
      printCentroids();
      if (checkCentoids(last, centroids_)) {
        break;
      }
    }
  }

  public void run(int n) {
    for (int epoch = 0; epoch < n; epoch++) {
      printCentroids();
      for (int i = 0; i < data_.getNumRows(); i++) {
        moveTowards(getClosestCentroid(data_.getRow(i)), data_.getRow(i)); 
      }
    }
  }

  public static void main(String args[]) {
    double data[][] = {{2}, {1}, {0}, {8}, {9}, {6}};
    Matrix m = new Matrix(6, 1);
    for (int i = 0; i < 6; i++) {
      m.setCell(i,0,data[i][0]);
    }
    // new KMeansClustering(2, m).run(Integer.parseInt(args[0]));
    new KMeansClustering(2, m).computeClusters();
  }
}
