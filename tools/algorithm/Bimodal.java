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

import tools.microarray.*;
import tools.microarray.StepMiner.*;
import java.util.*;

public class Bimodal {

    public static double GAP_LENGTH = 0.5;

    int length_;
    double threshold_;
    double statistic_;
    double entropy_;
    double lowthreshold_;
    double highthreshold_;

    public Bimodal(Double[] vector) {
      Double[] vector_ = vector;
      int count = getCount(vector_);
      length_ = count;
      statistic_ = 0.0;
      entropy_ = 0.0;
      threshold_ = 0.0;
      lowthreshold_ = 0.0;
      highthreshold_ = 0.0;
      if (count <= 0) {
        return;
      }
      try {
        class dComparator implements Comparator<Double> {
          public int compare(Double s1, Double s2) {
            if (s1 == null) return 1;
            if (s2 == null) return -1;
            if (s1.doubleValue() > s2.doubleValue()) return 1;
            return -1;
          }
        };
        Arrays.sort(vector_, new dComparator());
        AnalysisMetaData meta = new AnalysisMetaData(0, vector_.length-1, 
            null, null, 0.05);
        GeneData gene = new GeneData(vector_);
        // gene.print();
        Step step = OneStep.fitStep(gene, meta);
        int index = step.getStep(0);
        threshold_ = vector_[index].doubleValue();
        int next = index + 1;
        //System.out.println(index + " " + next + " " + threshold_);
        while ( next < count && vector_[index] == null) {
          next ++;
        }
        if (next < count && vector_[index] != null) {
          threshold_ = (threshold_ + vector_[next].doubleValue()) / 2;
        }
        //System.out.println(index + " " + next + " " + threshold_);
        int lowindex = index - index * 10/100;
        int highindex = index + (count-index) * 10/100;
        if (lowindex < 0) {
            lowindex = 0;
        }
        if (highindex >= count) {
            highindex = count-1;
        }
        //lowthreshold_ = vector_[lowindex].doubleValue();
        //highthreshold_ = vector_[highindex].doubleValue();
        lowthreshold_ = threshold_ - GAP_LENGTH;
        highthreshold_ = threshold_ + GAP_LENGTH;
        int count0 = getCount0(vector_, threshold_);
        int count1 = getCount1(vector_, threshold_);
        entropy_ = computeEntropy(count0, count1);
        // System.out.println(" Middle = " + threshold_ + " Index = " + index);
        statistic_ = step.getStatistic();
        // System.out.println(" T = " + statistic_ );
      }
      catch(Exception e) {
        // e.printStackTrace();
        System.out.println("Count = " + count);
      }
    }

    public double getThreshold() { return threshold_; }
    public double getLowThreshold() { return lowthreshold_; }
    public double getHighThreshold() { return highthreshold_; }
    public double getStatistic() { return statistic_; }
    public double getEntropy() { return entropy_; }

    public boolean isBimodal() {
      double t = statistic_;
      return (t > length_);
    }

    public boolean isBimodal(double thr) {
      double t = statistic_;
      return (t > thr);
    }

    public static double computeEntropy(int c0, int c1) {
        double p1 = c0 / (c0 + c1 + 0.0);
        double p2 = 1 - p1;
        return -p1 * Math.log(p1) - p2 * Math.log(p2);
    }

    public static int getCount(Double[] v) {
      int count = 0;
      for (int i =0; i < v.length; i++) {
        if (v[i] != null) {
          count++;
        }
      }
      return count;
    }

    public static int getCount0(Double[] v, double thr) {
      int count = 0;
      for (int i =0; i < v.length; i++) {
        if (v[i] != null && v[i].doubleValue() < thr) {
          count++;
        }
      }
      return count;
    }

    public static int getCount1(Double[] v, double thr) {
      int count = 0;
      for (int i =0; i < v.length; i++) {
        if (v[i] != null && v[i].doubleValue() >= thr) {
          count++;
        }
      }
      return count;
    }

    public static Matrix getMatrix(Double[] v) {
      int count = 0;
      for (int i =0; i < v.length; i++) {
        if (v[i] != null) {
          count++;
        }
      }
      Matrix m = new Matrix(count, 1);
      int index = 0;
      for (int i =0; i < v.length ; i++) {
        if (v[i] != null) {
          m.setCell(index++, 0, v[i].doubleValue());
        }
      }
      return m;
    }

    public static boolean check(Double[] v) {
        return new Bimodal(v).isBimodal();
    }

    public static boolean check(Double[] v, double thr) {
        return new Bimodal(v).isBimodal(thr);
    }

    public static boolean checkEntropy(Double[] v, double thr) {
        Bimodal b = new Bimodal(v);
        return b.isBimodal(thr) && b.getEntropy() > 0.5;
    }

    public static void main(String arg[]) throws Exception {
        Double[] v = new Double[300];
        Random r = new Random(100);
        byte[] d = new byte[300];
        r.nextBytes(d);
        int x = Integer.parseInt(arg[0]);
        for (int i =0; i < x; i++) {
            v[i] = new Double(d[i]/255.0);
        }
        for (int i =x; i < 300 ; i++) {
            v[i] = new Double(d[i]/255.0 + 2);
        }
        Bimodal b = new Bimodal(v);
        System.out.println(b.getThreshold());
        System.out.println(b.getStatistic());
        System.out.println(b.getEntropy());
    }

}

