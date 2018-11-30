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
import java.io.*;
import tools.goanalysis.*;

public class GeneSetAnalysis {

  GeneSet a_;
  GeneSet b_;

  static double[] logFact_;

  public GeneSetAnalysis(GeneSet a) {
    a_ = a;
    b_ = null;
  }

  public GeneSetAnalysis(GeneSet a, GeneSet b) {
    a_ = a;
    b_ = b;
  }

  public static double getPvalue(int k, int n, int M, int N) {
    double pval = 0;
    // Setup logFact cache
    if (logFact_ == null || logFact_.length < (N+2)) {
      logFact_ = new double[N+2];
      logFact_[0] = 0; logFact_[1] = 0;
      for (int i = 2; i < (N+2); i++) {
        logFact_[i] = logFact_[i-1] + Math.log(i);
      }
    }
    int min = n;
    if (min > M) {
      min = M;
    }
    double lo = logFact_[N] - logFact_[n] - logFact_[N-n];
    for (int i =k; i <=min; i++) {
      double upp1 = logFact_[M] - logFact_[i] - logFact_[M-i];
      double upp2 = logFact_[N-M] - logFact_[n-i] - logFact_[N-M-n+i];
      pval = pval + Math.exp(upp1 + upp2 - lo);
    }
    return pval;
  }

  public void performGOAnalysis(String ofile, String onnFile,
      String annFile, String org, double pvalue, 
      Vector<String> allNodes) throws Exception {
    GOAnalysis goa = new GOAnalysis(onnFile, annFile, org, pvalue);
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    out.write(GOAnalysis.getHeader());
    Iterator<String> itr1 = a_.iterator();
    out.write("<li> <font size=+1> Gene Set Analysis </font> ("+ a_.getTotalNum()  +")<ul>\n");
    while(itr1.hasNext()) {
      String str1 = (String) itr1.next();
      HashSet<String> set1 = a_.getSet(str1);
      out.write("<li> <font size=+1> " + str1 + "</font> <ul>\n");
      Vector<String> strSet = new Vector<String>(set1);
      for (int k = 0; k < strSet.size(); k++) {
        System.out.println("["+strSet.get(k)+"]");
      }
      goa.printGOTerms(out, org, strSet, allNodes);
      out.write("</ul> </li>\n");
    }
    out.write("</ul> </li>\n");
    out.write(GOAnalysis.getEnd());
    out.close();
  }

  public void performGOAnalysisOutTAB(String ofile, String onnFile,
      String annFile, String org, double pvalue, 
      Vector<String> allNodes) throws Exception {
    GOAnalysis goa = new GOAnalysis(onnFile, annFile, org, pvalue);
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    //goa.printGOTermsTAB(out, org, a_, allNodes);
    goa.printGOTermsTABAlt(out, org, a_, allNodes);
    out.close();
  }

  public void performAnalysis(String file, String org, double pvalue) throws Exception {
    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    out.write(GOAnalysis.getHeader());
    HashSet<String> all = new HashSet<String>(a_.getAllGenes());
    HashSet<String> pairs = new HashSet<String>();
    all.addAll(b_.getAllGenes());
    int N = all.size();
    out.write("<li> <font size=+1> Gene Set Analysis </font> ("+ a_.getTotalNum() + "," + b_.getTotalNum() + "," + N +")<ul>\n");
    Iterator<String> itr1 = a_.iterator();
    while(itr1.hasNext()) {
      String str1 = (String) itr1.next();
      HashSet<String> set1 = a_.getSet(str1);
      Iterator<String> itr2 = b_.iterator();
      StringBuffer outStr = new StringBuffer();
      boolean sig = false;
      outStr.append("<li> <font size=+1> " + str1 + "</font> <ul>\n");
      while(itr2.hasNext()) {
        String str2 = (String) itr2.next();
        if (str1.equals(str2)) {
          continue;
        }
        if (pairs.contains(str1+str2) || pairs.contains(str2+str1)) {
          continue;
        }
        pairs.add(str1+str2);
        HashSet<String> set2 = b_.getSet(str2);
        String desc = b_.getDescription(str2);
        HashSet<String> intersect = new HashSet<String>(set1);
        intersect.retainAll(set2);
        int k = intersect.size();
        int n = set1.size();
        int M = set2.size();
        double pval = getPvalue(k, n, M, N);
        if (pval < pvalue) {
          sig = true;
          //System.out.println(str1 + ", " + str2 + "," + desc);
          //System.out.println("\t" + pval + ":" + k + "," + n + "," + M + "," + N);
          outStr.append("<li> <font size=+1> "+ str2 +"</font> " +
              desc + " - " +
              "(" + pval + ":" + k + "," + n + "," + M + "," + N + ")" +
              "<ul> <li> <table border=0><tr>\n");
          Iterator<String> gitr = intersect.iterator();
          int count = -1;
          while (gitr.hasNext()) {
            String gene = (String) gitr.next();
            count ++;
            if ( (count % 10) == 0) {
              outStr.append("</tr><tr>\n");
            }
            String link = GOAnalysis.getLink(org, gene);
            outStr.append("<td> <a target=\"_blank\" href=\""+ link + "\"> "+ gene + "</a></td>\n");
          }
          outStr.append("</tr></table> </li></ul> </li>\n");
        }
      }
      outStr.append("</ul> </li>\n");
      if (sig) {
        out.write(outStr.toString());
      }
    }
    out.write("</ul> </li>\n");
    out.write(GOAnalysis.getEnd());

    out.close();
  }

  public void performAnalysisOutTAB(String file, String org, double pvalue) throws Exception {
    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    HashSet<String> all = new HashSet<String>(a_.getAllGenes());
    HashSet<String> pairs = new HashSet<String>();
    all.addAll(b_.getAllGenes());
    int N = all.size();
    // Header
    out.write("GeneSet\tDescription");
    Iterator<String> itr1 = a_.iterator();
    while(itr1.hasNext()) {
      String str1 = (String) itr1.next();
      out.write("\t" + str1);
    }
    out.write("\n");
    Iterator<String> itr2 = b_.iterator();
    while(itr2.hasNext()) {
      String str2 = (String) itr2.next();
      HashSet<String> set2 = b_.getSet(str2);
      StringBuffer outStr = new StringBuffer();
      boolean sig = false;
      itr1 = a_.iterator();
      while(itr1.hasNext()) {
        String str1 = (String) itr1.next();
        if (str1.equals(str2)) {
          continue;
        }
        if (pairs.contains(str1+str2) || pairs.contains(str2+str1)) {
          continue;
        }
        pairs.add(str1+str2);
        HashSet<String> set1 = a_.getSet(str1);
        HashSet<String> intersect = new HashSet<String>(set1);
        intersect.retainAll(set2);
        int k = intersect.size();
        int n = set1.size();
        int M = set2.size();
        double pval = getPvalue(k, n, M, N);
        if (pval < pvalue) {
          sig = true;
          //System.out.println(str1 + ", " + str2 + "," + desc);
          //System.out.println("\t" + pval + ":" + k + "," + n + "," + M + "," + N);
          outStr.append( "\t(" + k + "," + pval + ")");
        }
        else {
          outStr.append("\t");
        }
      }
      if (sig) {
        outStr.append("\n");
        String desc = b_.getDescription(str2);
        out.write(str2 + "\t" + desc + outStr.toString());
      }
    }
    out.close();
  }

  public static void main(String args[]) throws Exception {
    GeneSet a = GeneSet.readFile(args[0]);
    GeneSet b = GeneSet.readFile(args[1]);
    GeneSetAnalysis ana = new GeneSetAnalysis(a, b);
    ana.performAnalysis(args[2], "Hs", 0.001);
  }

}
