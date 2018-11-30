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
import java.awt.Color;

import tools.microarray.StepMiner.*;
import tools.microarray.FileReader.*;
import tools.microarray.FileWriter.*;

import tools.algorithm.Bimodal;

import tools.graphs.*;

public class GroupAnalysis {

  public static int countThreshold(int index, double thr1, double thr2, Double[] v1, Double[] v2) {
    int count = 0;
    for (int i = 0; i < v1.length; i++) {
      if (v1[i] != null && v2[i] != null) {
        if (index == 0 && v1[i].doubleValue() <= thr1 && v2[i].doubleValue() <= thr2) count++;
        if (index == 1 && v1[i].doubleValue() <= thr1 && v2[i].doubleValue() > thr2) count++;
        if (index == 2 && v1[i].doubleValue() > thr1 && v2[i].doubleValue() > thr2) count++;
        if (index == 3 && v1[i].doubleValue() > thr1 && v2[i].doubleValue() <= thr2) count++;
      }
    }
    return count;
  }

  /*
   * returns 0 if positive, 1 if negative, -1 if none.
   */
  public static int getStatus(double thr1, double thr2, Double[] v1, Double[] v2) {
    int c0 = countThreshold(0, thr1, thr2, v1, v2);
    int c1 = countThreshold(1, thr1, thr2, v1, v2);
    int c2 = countThreshold(2, thr1, thr2, v1, v2);
    int c3 = countThreshold(3, thr1, thr2, v1, v2);
    if ((c0 + c2) > (v1.length * 0.90)) {
        return 0; /* positive */
    }
    if ((c1 + c3) > (v1.length * 0.90)) {
        return 1; /* negative */
    }
    return -1;
  }

  public static void groupAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String ifile = list.removeFirst();
    Data data1 = PCLFileReader.readFile(ifile);

    Bimodal[] b = new Bimodal[data1.getNumGenes()];
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      // System.out.println("[ " + i + " ]---------------------->");
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v= data1.getGeneData(i).getVector(start, end);
      b[i - data1.getNumGeneHeader()] = new Bimodal(v);
    }

    HashSet<Group> groups = new HashSet<Group>();
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      if ( (i%100) == 0) {
        System.out.println("[ " + i + " ]");
      }
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v1 = data1.getGeneData(i).getVector(start, end);
      int i1 = i - data1.getNumGeneHeader();
      double thr1 = b[i1].getThreshold();
      boolean found = false;
      Iterator<Group> itr = groups.iterator();
      while (itr.hasNext()) {
        Group next = (Group) itr.next();
        if (next.hasPositive()) {
          int i2 = next.getPositive();
          int j = i2 + data1.getNumGeneHeader();
          double thr2 = b[i2].getThreshold();
          Double[] v2 = data1.getGeneData(j).getVector(start, end);
          int status = getStatus(thr1, thr2, v1, v2);
          if (status == 0) {
            next.addPositive(i1); found = true;
          }
          if (status == 1) {
            next.addNegative(i1); found = true;
          }
        }
        if (next.hasNegative()) {
          int i2 = next.getNegative();
          int j = i2 + data1.getNumGeneHeader();
          double thr2 = b[i2].getThreshold();
          Double[] v2 = data1.getGeneData(j).getVector(start, end);
          int status = getStatus(thr1, thr2, v1, v2);
          if (status == 0) {
            next.addNegative(i1); found = true;
          }
          if (status == 1) {
            next.addPositive(i1); found = true;
          }
        }
        if (found == true) {
            break;
        }
      }
      if (!found) {
        // Add a new group
        Group res = new Group();
        res.addPositive(i1);
        groups.add(res);
      }
    }
    System.out.println("Number of groups = " + groups.size());
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    Iterator<Group> itr = groups.iterator();
    while (itr.hasNext()) {
      Group next = (Group) itr.next();
      next.write(out);
    }
    out.close();
  }

}

class Group {
  TreeSet<Integer> positive;
  TreeSet<Integer> negative;

  public Group() {
    positive = new TreeSet<Integer>();
    negative = new TreeSet<Integer>();
  }

  public int getPositive() throws NoSuchElementException {
    Integer res = positive.first();
    return res.intValue();
  }
  public int getNegative() throws NoSuchElementException {
    Integer res = negative.first();
    return res.intValue();
  }
  public void addPositive(int a) { positive.add(new Integer(a)); }
  public void addNegative(int a) { negative.add(new Integer(a)); }
  public boolean hasPositive() { return !positive.isEmpty(); }
  public boolean hasNegative() { return !negative.isEmpty(); }

  public void write(BufferedWriter out) throws IOException {
    String name = "none";
    int pos = -1;
    if (hasPositive()) {
      pos = 0;
      name = "" + getPositive();
    }
    else if (hasNegative()) {
      pos = 1;
      name = "" + getNegative();
    }
    if (hasPositive()) {
      // Write positive:
      String n = name;
      if (pos == 0) { n += "+"; }
      if (pos == 1) { n += "-"; }
      String size = "" + positive.size();
      out.write(n + "\t" + size);
      Iterator<Integer> itr = positive.iterator();
      while (itr.hasNext()) {
        Integer next = (Integer) itr.next();
        out.write("\t" + next);
      }
      out.write("\n");
    }
    if (hasNegative()) {
      // Write negative:
      String n = name;
      if (pos == 0) { n += "-"; }
      if (pos == 1) { n += "+"; }
      String size = "" + negative.size();
      out.write(n + "\t" + size);
      Iterator<Integer> itr = negative.iterator();
      while (itr.hasNext()) {
        Integer next = (Integer) itr.next();
        out.write("\t" + next);
      }
      out.write("\n");
    }
  }
}
