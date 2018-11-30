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
import tools.microarray.StepMiner.*;
import tools.microarray.FileWriter.*;
import tools.goanalysis.GOAnalysis;
import tools.graphs.PSPlot;

public class ArrayOrder {

  Data data_;
  int start_;
  int end_;
  Vector<Integer[]> sortedPerm_;
  Vector<Integer> geneIndex_;
  Vector<Integer> sortedGeneIndex_;
  Vector<Integer> missingGeneIndex_;
  int[] orderedPerm_;
  HashSet<Integer> unorderedIndices_;

  Vector<Integer[]> sortingPermKL_;

  public ArrayOrder(Data data) {
    data_ = data;
    start_ = data_.getNumArrayHeader();
    end_ = data_.getNumColumns() - 1;
    sortedPerm_ = new Vector<Integer[]>();
    sortingPermKL_ = new Vector<Integer[]>();
    geneIndex_ = new Vector<Integer>();
    sortedGeneIndex_ = new Vector<Integer>();
    missingGeneIndex_ = new Vector<Integer>();
    orderedPerm_ = new int[end_ - start_ + 1];
    for (int i =0; i < orderedPerm_.length; i++) {
      orderedPerm_[i] = i;
    }
    unorderedIndices_ = new HashSet<Integer>();
  }

  public int[] getOrderedPermutation() { return orderedPerm_; }
  public Data getData() { return data_; }

  public void populateOrder() {
    int numMissing = data_.getNumMissingPoints();
    for (int i = data_.getNumGeneHeader(); i < data_.getNumRows(); i++) {
      GeneData gene = data_.getGeneData(i);
      gene.removeZero(start_, end_);
      if (gene.getMissingPoints(start_, end_) <= numMissing) {
        Integer[] order1 = gene.getSortedOrderAsc(start_, orderedPerm_);
        Integer[] order2 = gene.getSortedOrderDes(start_, orderedPerm_);
        sortedPerm_.add(order1);
        sortedPerm_.add(order2);
        geneIndex_.add(new Integer(i));
      }
      else {
        missingGeneIndex_.add(new Integer(i));
      }
    }
  }

  public int[] getCounts(int index) {
    int[] res = new int[end_ - start_ +1];
    for (int i =0; i < res.length; i++) {
      res[i] = 0;
    }
    Enumeration<Integer[]> e = sortedPerm_.elements();
    while (e.hasMoreElements()) {
      Integer[] order = (Integer[]) e.nextElement();
      res[order[index].intValue()] ++;
    }
    return res;
  }

  public int[] getCountsConstraint(int cindex, int val, int index) {
    int[] res = new int[end_ - start_ +1];
    for (int i =0; i < res.length; i++) {
      res[i] = 0;
    }
    Enumeration<Integer[]> e = sortedPerm_.elements();
    while (e.hasMoreElements()) {
      Integer[] order = (Integer[]) e.nextElement();
      if (order[cindex].intValue() == val) {
        res[order[index].intValue()] ++;
      }
    }
    return res;
  }

  public static int getMax(int[] arr, HashSet<Integer> list) {
    int max = 0;
    int maxIndex = 0;
    for (int i =0; i < arr.length; i++) {
      if (!list.contains(new Integer(i)) && max < arr[i]) {
        max = arr[i];
        maxIndex = i;
      }
    }
    return maxIndex;
  }

  public int getNumMatch(int[] perm1, Integer[] perm2) {
    int count = 0;
    for (int i =0; i < perm1.length; i++) {
      if (perm2[i] != null && perm2[i].intValue() == perm1[i]) {
        count ++;
      }
    }
    return count;
  }

  public int getSortedCount() {
    int count = 0;
    Enumeration<Integer[]> e = sortedPerm_.elements();
    while (e.hasMoreElements()) {
      Integer[] order = (Integer[]) e.nextElement();
      if (getNumMatch(orderedPerm_, order) == order.length ) {
        count ++;
      }
    }
    return count;
  }

  public double getCorrelation(GeneData gene, int[] perm1, Integer[] perm2) {
    double sum_xy = 0, sum_x = 0, sum_y = 0, sum_sqx = 0, sum_sqy = 0;
    int count = 0;
    double res =0;
    for (int i =0; i < perm1.length; i++) {
      if (perm2[i] != null) {
        Double x = (Double) gene.getDataAt(perm1[i]+start_);
        Double y = (Double) gene.getDataAt(perm2[i].intValue()+start_);
        if (x != null && y != null) {
            count ++;
            sum_xy += x.doubleValue() * y.doubleValue();
            sum_x += x.doubleValue();
            sum_y += y.doubleValue();
            sum_sqx += x.doubleValue() * x.doubleValue();
            sum_sqy += y.doubleValue() * y.doubleValue();
        }
      }
    }
    if (count != 0) {
        res = (sum_xy - 1.0/count * sum_x * sum_y)/
            Math.sqrt(sum_sqx - 1.0/count * sum_x * sum_x)/
            Math.sqrt(sum_sqy - 1.0/count * sum_y * sum_y);
    }
    if (Double.isNaN(res)) {
        res = 0.0;
    }
    return res;
  }

  public static Vector<Integer> sortCorrelation(Vector<Integer> a, 
    Vector<Double> corr) 
  {
    Integer[] order = new Integer[a.size()];
    for (int j =0; j < order.length ; j++) {
      order[j] = new Integer(j);
    }
    class CorrComparator implements Comparator<Integer> {
      Vector<Double> corr_;
      public CorrComparator(Vector<Double> corr) {
        corr_ = corr;
      }
      public int compare(Integer s1, Integer s2) {
        Double c1 = corr_.get(s1.intValue());
        Double c2 = corr_.get(s2.intValue());
        if (c1.doubleValue() < c2.doubleValue()) {
          return 1;
        }
        return -1;
      }
    };
    Arrays.sort(order, new CorrComparator(corr));
    Vector<Integer> res = new Vector<Integer>();
    for (int j =0; j < order.length ; j++) {
      res.add(a.get(order[j].intValue()));
    }
    return res;
  }

  public double getStatistic() {
    double res = 0;
    Vector<Integer> up = new Vector<Integer>();
    Vector<Integer> down = new Vector<Integer>();
    Vector<Double> upCorr = new Vector<Double>();
    Vector<Double> downCorr = new Vector<Double>();
    sortedGeneIndex_ = new Vector<Integer>();
    int num = orderedPerm_.length;
    int count = 0;
    for (int i =0; i < geneIndex_.size(); i++) {
      Integer g = geneIndex_.get(i);
      GeneData gene = data_.getGeneData(g.intValue());
      Integer[] orderAsc = gene.getSortedOrderAsc(start_, orderedPerm_);
      Integer[] orderDes = gene.getSortedOrderDes(start_, orderedPerm_);
      double cor1 = getCorrelation(gene, orderedPerm_, orderAsc);
      double cor2 = getCorrelation(gene, orderedPerm_, orderDes);
      double pvalue, f;
      int invCount = 0;
      if (cor1 > cor2) {
        up.add(g);
        upCorr.add(new Double(cor1));
        res += cor1 * cor1 * cor1 * cor1;
        f = (1+cor1)/(1-cor1);
        invCount = gene.countInversion(start_, orderedPerm_, true);
      }
      else {
        down.add(g);
        downCorr.add(new Double(-cor2));
        res += cor2 * cor2 * cor2 * cor2;
        f = (1+cor2)/(1-cor2);
        invCount = gene.countInversion(start_, orderedPerm_, false);
      }
      pvalue = 1 - tools.microarray.StepMiner.Utils.Fisher(f, num-2, num-2);
      if (pvalue < 0.001) {
        //count++;
      }
      if (invCount < 10) {
        count++;
      }
    }
    System.out.println("Count : " + count);
    res = res / geneIndex_.size();
    up = sortCorrelation(up, upCorr);
    down = sortCorrelation(down, downCorr);
    sortedGeneIndex_.addAll(up);
    sortedGeneIndex_.addAll(down);
    return res;
  }

  public void calculateOrderOld() {
    populateOrder();
    System.out.println("Statistic : " + getStatistic());
    HashSet<Integer> list = new HashSet<Integer>();
    orderedPerm_[0] = getMax(getCounts(0), list);
    list.add(new Integer(orderedPerm_[0]));
    System.out.print("Sequence: " + orderedPerm_[0]);
    for (int i =1; i < orderedPerm_.length; i++) {
      orderedPerm_[i] = getMax(getCountsConstraint(i-1, orderedPerm_[i-1], i), list);
      list.add(new Integer(orderedPerm_[i]));
      System.out.print(" " + orderedPerm_[i]);
    }
    System.out.println();
    // prints header in order
    GeneData header = data_.getGeneData(0);
    for (int i = 0; i < orderedPerm_.length; i++) {
      String n = (String) header.getDataAt(orderedPerm_[i]+start_);
      System.out.println(n);
    }
    System.out.println("Monotonic count : " + getSortedCount());
    System.out.println("Statistic : " + getStatistic());
  }

  boolean evaluateArrayPlacing(int i, int j , int k) {
    int countMonotonic = 0;
    int countNonMonotonic = 0;
    int count = 0;
    for (int ii = data_.getNumGeneHeader(); ii < data_.getNumRows(); ii++) {
      GeneData gene = data_.getGeneData(ii);
      Double vi = (Double) gene.getDataAt(start_ + i);
      Double vj = (Double) gene.getDataAt(start_ + j);
      Double vk = (Double) gene.getDataAt(start_ + k);
      if (vi != null && vj != null && vk != null) {
        // Swap if vi > vk
        if (vi.doubleValue() > vk.doubleValue()) {
          Double tmp = vi;
          vi = vk;
          vk = tmp;
        }
        count++;
        if (vj.doubleValue() > vk.doubleValue() || 
            vj.doubleValue() < vi.doubleValue()) {
          countNonMonotonic++;
        }
        else {
          countMonotonic++;
        }
      }
    }
    // System.out.println(count+ " M: " + countMonotonic + " N: " + countNonMonotonic);
    return ( (2 * countMonotonic) < countNonMonotonic );
  }

  /*
   * Start traversing from the both ends  (Aug 17, 2006)
   */
  public void calculateOrderNew() {
    populateOrder();
    System.out.println("Statistic (Before-"+ orderedPerm_.length+ "): " + getStatistic());
    HashSet<Integer> list = new HashSet<Integer>();
    orderedPerm_[0] = getMax(getCounts(0), list);
    list.add(new Integer(orderedPerm_[0]));
    int lastindex = orderedPerm_.length - 1;
    orderedPerm_[lastindex] = getMax(getCountsConstraint(0, orderedPerm_[0], lastindex), list);
    list.add(new Integer(orderedPerm_[lastindex]));
    int num = 2; // Number of arrays considered
    if (lastindex == 0) {
        num = 1;
    }
    int mid = orderedPerm_.length / 2;
    for (int i =1; i < mid; i++) {
      orderedPerm_[i] = getMax(getCountsConstraint(i-1, orderedPerm_[i-1], i), list);
      boolean fail = evaluateArrayPlacing(orderedPerm_[i-1], orderedPerm_[i], orderedPerm_[lastindex]);
      if (fail) {
        orderedPerm_[i] = orderedPerm_[i-1];
      }
      else {
        list.add(new Integer(orderedPerm_[i]));
        num++;
      }
    }
    for (int i = lastindex - 1; i >= mid; i--) {
      orderedPerm_[i] = getMax(getCountsConstraint(i+1, orderedPerm_[i+1], i), list);
      boolean fail = evaluateArrayPlacing(orderedPerm_[mid-1], orderedPerm_[i], orderedPerm_[i+1]);
      if (fail) {
        orderedPerm_[i] = orderedPerm_[i+1];
      }
      else {
        list.add(new Integer(orderedPerm_[i]));
        num++;
      }
    }
    int[] perm = new int[num];
    int last = -1;
    int index = 0;
    System.out.print("Sequence: ");
    for (int i = 0; i < orderedPerm_.length; i++) {
      if (orderedPerm_[i] != last) {
        perm[index++] = orderedPerm_[i];
        System.out.print(" " + orderedPerm_[i]);
      }
      last = orderedPerm_[i];
    }
    System.out.println();

    orderedPerm_ = perm; // change the order permutation

    // Computing other indices that are not used.
    HashSet<Integer> setAll = new HashSet<Integer>();
    int length = end_ - start_ + 1;
    for (int i = 0; i < length; i++) {
      setAll.add(new Integer(i));
    }
    unorderedIndices_ = new HashSet<Integer>();
    for (int i = 0; i < orderedPerm_.length; i++) {
      unorderedIndices_.add(new Integer(orderedPerm_[i]));
    }
    setAll.removeAll(unorderedIndices_);
    unorderedIndices_ = setAll;

    // prints header in order
    GeneData header = data_.getGeneData(0);
    for (int i = 0; i < orderedPerm_.length; i++) {
      String n = (String) header.getDataAt(orderedPerm_[i]+start_);
      System.out.println(n);
    }

    System.out.println("Monotonic count : " + getSortedCount());
    System.out.println("Statistic (After-"+ orderedPerm_.length + "): " + getStatistic());
  }

  public void calculateOrderIter() {
    calculateOrderNew();
    Data d = getOrderedDataSimple();
    d = getUnOrderedDataSimple();
    while (unorderedIndices_.size() > 0) {
        ArrayOrder arr = new ArrayOrder(d);
        arr.calculateOrderNew();
        unorderedIndices_ = arr.unorderedIndices_;
        d = arr.getUnOrderedDataSimple();
    }
  }

  int computeNumEdges(int a, int b) {
    int count = 0;
    Enumeration<Integer[]> e = sortingPermKL_.elements();
    while (e.hasMoreElements()) {
      Integer[] order = (Integer[]) e.nextElement();
      if (order[a].intValue() == b || order[b].intValue() == a) {
        count ++;
      }
    }
    return count;
  }

  public void populateOrderKL() {
    int numMissing = data_.getNumMissingPoints();
    for (int j = data_.getNumGeneHeader(); j < data_.getNumRows(); j++) {
      GeneData gene = data_.getGeneData(j);
      gene.removeZero(start_, end_);
      if (gene.getMissingPoints(start_, end_) <= numMissing) {
        Integer[] order1 = gene.getSortedOrderAsc(start_, orderedPerm_);
        Integer[] order_f = new Integer[order1.length];
        Integer[] order_b = new Integer[order1.length];
        for (int i =0; i < order1.length; i++) {
            int f = -1, b = -1;
            if ( (i+1) < order1.length ) {
                f = order1[(i+1)].intValue();
            }
            if ( (i-1) >= 0 ) {
                b = order1[(i-1)].intValue();
            }
            order_f[order1[i].intValue()] = new Integer(f);
            order_b[order1[i].intValue()] = new Integer(b);
        }
        sortingPermKL_.add(order_f);
        sortingPermKL_.add(order_b);
        geneIndex_.add(new Integer(j));
      }
      else {
        missingGeneIndex_.add(new Integer(j));
      }
    }
  }

  public void calculateOrderKL() {
    populateOrderKL();
    System.out.println("Statistic (Before-"+ orderedPerm_.length+ "): " + getStatistic());
    HashSet<Integer> g1 = new HashSet<Integer>();
    HashSet<Integer> g2 = new HashSet<Integer>();
    int mid = orderedPerm_.length / 2;
    for (int i =0; i < mid; i++) {
      g1.add(new Integer(i));
    }
    for (int i =mid; i < orderedPerm_.length; i++) {
      g2.add(new Integer(i));
    }
    int bestCost = Integer.MAX_VALUE;
    boolean gain = true;
    int itrCount = 0;
    while (gain) {
      itrCount++;
      System.out.println("Iteration #" + itrCount);
      int[] Dn = new int[orderedPerm_.length];
      int[] mark = new int[orderedPerm_.length];
      int[] cost = new int[orderedPerm_.length];
      cost[0] = 0;
      Iterator<Integer> itr1 = g1.iterator();
      while (itr1.hasNext()) {
        Integer a = (Integer) itr1.next();
        int c = 0;
        Iterator<Integer> itr2 = g2.iterator();
        while (itr2.hasNext()) {
          Integer b = (Integer) itr2.next();
          c += computeNumEdges(a.intValue(), b.intValue());
        }
        cost[0] += c;
        Dn[a.intValue()] = c;
      }
      itr1 = g2.iterator();
      while (itr1.hasNext()) {
        Integer a = (Integer) itr1.next();
        int c = 0;
        Iterator<Integer> itr2 = g1.iterator();
        while (itr2.hasNext()) {
          Integer b = (Integer) itr2.next();
          c += computeNumEdges(a.intValue(), b.intValue());
        }
        Dn[a.intValue()] = c;
      }
      itr1 = g1.iterator();
      while (itr1.hasNext()) {
        Integer a = (Integer) itr1.next();
        int c = 0;
        Iterator<Integer> itr2 = g1.iterator();
        while (itr2.hasNext()) {
          Integer b = (Integer) itr2.next();
          c += computeNumEdges(a.intValue(), b.intValue());
        }
        Dn[a.intValue()] -= c;
      }
      itr1 = g2.iterator();
      while (itr1.hasNext()) {
        Integer a = (Integer) itr1.next();
        int c = 0;
        Iterator<Integer> itr2 = g2.iterator();
        while (itr2.hasNext()) {
          Integer b = (Integer) itr2.next();
          c += computeNumEdges(a.intValue(), b.intValue());
        }
        Dn[a.intValue()] -= c;
      }
      for (int i =0; i < orderedPerm_.length; i++) {
        mark[i] = 0;
      }
      boolean existUnmarked = false;
      for (int i =0; i < orderedPerm_.length; i++) {
        if (mark[i] == 0) {
          existUnmarked = true;
        }
      }
      int bestChange = 0;
      int index = 0;
      if (cost[index] < bestCost) {
        gain = true;
        bestChange = index;
        bestCost = cost[index];
      }
      class Pair {
        public int i;
        public int j;
        public Pair(int a, int b) {
          i = a; j = b;
        }
      };
      gain = false;
      Pair[] pairs = new Pair[mark.length];
      while (existUnmarked) {
        int max = 0, max_i=0, max_j=0;
        for (int i =0; i < mark.length; i++) {
          for (int j = i + 1 ; j < mark.length; j++) {
            if (mark[i] == 0 && mark[j] == 0) {
              int e = computeNumEdges(i, j);
              int g = Dn[i] + Dn[j] - 2 * e;
              if (max < g) {
                max = g; max_i = i; max_j = j;
              }
            }
          }
        }
        // System.out.println(max_i + " " + max_j + " " + mark.length);
        mark[max_i] = 1; mark[max_j] = 1;
        index++;
        pairs[index] = new Pair(max_i, max_j);
        cost[index] = cost[index-1] - max;
        for (int i =0; i < mark.length; i++) {
          if (mark[i] == 0) {
            int e1 = computeNumEdges(i, max_j);
            int e2 = computeNumEdges(i, max_i);
            if (g1.contains(new Integer(i))) {
              Dn[i] = Dn[i] - 2 * e1 + 2 * e2;
            }
            else {
              Dn[i] = Dn[i] - 2 * e2 + 2 * e1;
            }
          }
        }
        if (cost[index] < bestCost) {
          gain = true;
          bestChange = index;
          bestCost = cost[index];
        }
        existUnmarked = false;
        for (int i =0; i < mark.length; i++) {
          if (mark[i] == 0) {
            existUnmarked = true;
          }
        }
        if (index == mark.length/2) {
          break;
        }
      }
      for (int s = 1; s < bestChange; s++) {
        g1.remove(new Integer(pairs[s].i));
        g2.remove(new Integer(pairs[s].j));
        g1.add(new Integer(pairs[s].j));
        g2.add(new Integer(pairs[s].i));
      }
      index = 0;
      itr1 = g1.iterator();
      while (itr1.hasNext()) {
        Integer a = (Integer) itr1.next();
        orderedPerm_[index++] = a.intValue();
      }
      itr1 = g2.iterator();
      while (itr1.hasNext()) {
        Integer a = (Integer) itr1.next();
        orderedPerm_[index++] = a.intValue();
      }
    }
    
    // prints header in order
    GeneData header = data_.getGeneData(0);
    for (int i = 0; i < orderedPerm_.length; i++) {
      String n = (String) header.getDataAt(orderedPerm_[i]+start_);
      System.out.println(n);
    }

    System.out.println("Statistic (After-"+ orderedPerm_.length + "): " + getStatistic());
  }

  public void calculateOrder() {
    calculateOrderIter();
  }

  public Data getUnOrderedDataSimple() {
    GeneData[] genes = new GeneData[data_.getNumRows()];
    int[] permHead = new int[data_.getNumArrayHeader()];
    for (int i = 0; i < permHead.length; i++) {
      permHead[i] = i;
    }
    for (int i = 0; i < data_.getNumRows(); i++) {
      GeneData gene = data_.getGeneData(i);
      GeneData head = gene.permute(0, permHead);
      GeneData data = gene.subset(start_, unorderedIndices_);
      genes[i] = GeneData.merge(head, data);
    }
    Data res = new Data(unorderedIndices_.size(), data_.getNumGenes(), 
        data_.getNumGeneHeader(), permHead.length, genes);
    res.setGeneNameScheme(data_.getGeneNameScheme());
    return res;
  }

  public Data getOrderedDataSimple() {
    GeneData[] genes = new GeneData[data_.getNumRows()];
    int[] permHead = new int[data_.getNumArrayHeader()];
    for (int i = 0; i < permHead.length; i++) {
      permHead[i] = i;
    }
    for (int i = 0; i < data_.getNumRows(); i++) {
      GeneData gene = data_.getGeneData(i);
      GeneData head = gene.permute(0, permHead);
      GeneData data = gene.permute(start_, orderedPerm_);
      genes[i] = GeneData.merge(head, data);
    }
    Data res = new Data(orderedPerm_.length, data_.getNumGenes(), 
        data_.getNumGeneHeader(), permHead.length, genes);
    res.setGeneNameScheme(data_.getGeneNameScheme());
    return res;
  }

  public Data getOrderedDataSorted() {
    GeneData[] genes = new GeneData[data_.getNumRows()];
    int[] permHead = new int[data_.getNumArrayHeader()];
    for (int i = 0; i < permHead.length; i++) {
      permHead[i] = i;
    }
    for (int i = 0; i < data_.getNumGeneHeader(); i++) {
      GeneData gene = data_.getGeneData(i);
      GeneData head = gene.permute(0, permHead);
      GeneData data = gene.permute(start_, orderedPerm_);
      genes[i] = GeneData.merge(head, data);
    }
    int index = data_.getNumGeneHeader();
    Enumeration<Integer> e = sortedGeneIndex_.elements();
    while (e.hasMoreElements()) {
      Integer g = (Integer) e.nextElement();
      GeneData gene = data_.getGeneData(g.intValue());
      GeneData head = gene.permute(0, permHead);
      GeneData data = gene.permute(start_, orderedPerm_);
      genes[index++] = GeneData.merge(head, data);
    }
    e = missingGeneIndex_.elements();
    while (e.hasMoreElements()) {
      Integer g = (Integer) e.nextElement();
      GeneData gene = data_.getGeneData(g.intValue());
      GeneData head = gene.permute(0, permHead);
      GeneData data = gene.permute(start_, orderedPerm_);
      genes[index++] = GeneData.merge(head, data);
    }
    if (index != genes.length) {
      throw new AssertionError("index != genes.length");
    }
    Data res = new Data(orderedPerm_.length, data_.getNumGenes(), 
        data_.getNumGeneHeader(), permHead.length, genes);
    res.setGeneNameScheme(data_.getGeneNameScheme());
    return res;
  }

  public Data getOrderedData() {
    return getOrderedDataSimple();
  }

  public void normalize() {
    data_.normalize();
  }

  public void meanCenter() {
    for (int i = data_.getNumGeneHeader(); i < data_.getNumRows(); i++) {
      try {
        GeneData gene = data_.getGeneData(i);
        double center = GeneData.getMean(gene.getData(), start_, end_);
        gene.performCentering(data_.getNumArrayHeader(), center);
      }
      catch(Exception e) {
        // No Action
      }
    } // end for
  }

  public void zeroCenter() {
    data_.zeroCenter();
  }

 public static void performGOAnalysis(Data data, String file, Double pvalThr) 
   throws Exception {
    GeneNameScheme ns = data.getGeneNameScheme();
    GOAnalysis goa = new GOAnalysis(ns.getOntologyFile(),
        ns.getAnnotationFile(), ns.getOrg(), pvalThr.doubleValue());
    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    out.write(GOAnalysis.getHeader());
    out.write("<li> <font size=+2> Group </font> <ul>\n");
    int numGenes = data.getNumGenes();
    int block = 500;
    int head = data.getNumGeneHeader();
    int end = data.getNumRows();
    Vector<String> geneSetAll = new Vector<String>();
    for (int i = head ; i < end; i++) {
      String gene = data.getGenesAt(i);
      geneSetAll.add(gene);
    } // end for
    Vector<String> geneSet = new Vector<String>();
    for (int i = head ; i < end; i++) {
      String gene = data.getGenesAt(i);
      geneSet.add(gene);
      if ( ( i % block ) == 0 || i == (end - 1)) {
        out.write("<li> <font size=+2> Block = "+(i/block)+" </font> <ul>\n");
        System.out.println("Block : " + (i/block));
        for (int j = 0; j < geneSet.size(); j++) {
          String g = geneSet.get(j);
          System.out.println("[" + g + "]");
        }
        goa.printGOTerms(out, ns.getOrg(), geneSet, geneSetAll);
        out.write("</ul> </li>\n");
        geneSet = new Vector<String>();
      }
    } // end for
    out.write("</ul> </li>\n");
    out.write(GOAnalysis.getEnd());
    out.close();
 }

 public static HashSet<Integer> getSelectedLabels(String tag) 
    throws NumberFormatException {
    HashSet<Integer> set = new HashSet<Integer>();
    if (tag.startsWith("B")) {
        tag = tag.replaceAll("B", "");
        String[] list = tag.split("-");
        if (list.length == 1) {
            int num = Integer.parseInt(list[0]);
            set.add(new Integer(num-1));
        }
        if (list.length == 2) {
            int start = Integer.parseInt(list[0]);
            int end = Integer.parseInt(list[1]);
            start--; end--;
            if (end < start) {
                end = start;
            }
            for (int i =start; i <= end; i++) {
              set.add(new Integer(i));
            }
        }
        if (list.length == 3) {
            int start = Integer.parseInt(list[0]);
            int end = Integer.parseInt(list[1]);
            int incr = Integer.parseInt(list[2]);
            start--; end--;
            if (end < start) {
                end = start;
            }
            if (incr <= 0) {
                incr = 1;
            }
            for (int i =start; i <= end; i+=incr) {
              set.add(new Integer(i));
            }
        }
    }
    System.out.print("Selected tags(" + tag + "): ");
    Iterator<Integer> itr = set.iterator();
    while (itr.hasNext()) {
      Integer n = (Integer) itr.next();
      System.out.print("" + (n.intValue() + 1) + " ");
    }
    System.out.println();
    return set;
 }

  /* 
   * Plotting Genes
   */
 public static void plotGenes(Data data, String file) throws Exception {
   String filename = Data.getFileName(file);
   String filetag = Data.getFileTag(file);
   HashSet<Integer> selectedLabels = ArrayOrder.getSelectedLabels(filetag);
   PSPlot plotter = new PSPlot(filename);
   plotter.open();
   plotter.array(5 /* Rows */, 2 /* columns */);
   int numGenes = data.getNumGenes();
   int block = 500;
   int startA = data.getNumArrayHeader(); // Array Index
   int endA = data.getNumColumns() - 1; // Array Index
   int head = data.getNumGeneHeader(); // Column index
   int end = data.getNumRows(); // Column index
   for (int i = head ; i < end; i++) {
     String geneName = data.getGenesAt(i);
     GeneData gene = data.getGeneData(i);
     int blkindex = (i/block);
     if (!selectedLabels.contains(new Integer(blkindex))) {
       if (!filetag.equals("All")) {
         continue;
       }
     }
     gene.convertDouble(startA, endA);
     Object[] geneData = gene.getData();
     Vector<Double> dat = new Vector<Double>();
     Vector<Double> time = new Vector<Double>();
     for (int j = startA; j <= endA; j++) {
       dat.add((Double)geneData[j]);
       time.add( new Double(j-startA));
     }
     plotter.plot(time, dat);
     plotter.ylabel(" Gene expression ");
     plotter.xlabel(geneName);
     plotter.title(geneName);
   } // end for
   plotter.close();
 }

 public void performRandomExpt() {
    Data data = (Data)data_.clone();
    Random random = new Random(100);
    data.scramble(random);
    ArrayOrder ord = new ArrayOrder(data);
    ord.calculateOrderOld();
    try {
      PCLFileWriter.writeFile(ord.getOrderedDataSorted(), "label1.pcl", null);
      StepMiner sm = new StepMiner(data_);
      sm.setOneStepAnalysis();
      sm.setPvalueThr(0.05);
      sm.performAnalysis();
      sm = new StepMiner(getOrderedData());
      sm.setOneStepAnalysis();
      sm.setPvalueThr(0.05);
      sm.performAnalysis();
      sm = new StepMiner(data);
      sm.setOneStepAnalysis();
      sm.setPvalueThr(0.05);
      sm.performAnalysis();
      sm = new StepMiner(ord.getOrderedData());
      sm.setOneStepAnalysis();
      sm.setPvalueThr(0.05);
      sm.performAnalysis();
    }
    catch(Exception e) {
      e.printStackTrace();
    }
 }

 public Data listGenes() {
   int numArrays = 0;
   int numGenes = data_.getNumGenes();
   int numGeneHeader = 2;
   int numArrayHeader = 1;
   GeneNameScheme nameScheme = data_.getGeneNameScheme();
   nameScheme.print();

   GeneData[] data = new GeneData[numGenes+numGeneHeader];
   int index = 0;
   Object[] obj = new Object[1];
   obj[0] = "Genes";
   data[index++] = new GeneData(obj);
   obj = new Object[1];
   obj[0] = "EWEIGHT";
   data[index++] = new GeneData(obj);
   for (int i = data_.getNumGeneHeader(); i < data_.getNumRows(); i++) {
     String geneName = data_.getGenesAt(i);
     obj = new Object[1];
     obj[0] = geneName;
     data[index++] = new GeneData(obj);
   }
   Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
   res.setGeneNameScheme(nameScheme);
   return res;
 }

 public static void main(String[] args) throws Exception {
   Random random = new Random(100);
   Data data = Data.getArtificialData(random, 10000, Integer.parseInt(args[0]));
   ArrayOrder ord = new ArrayOrder(data);
   ord.calculateOrderOld();
   ord.performRandomExpt();
   tools.microarray.FileWriter.PCLFileWriter.writeFile(data, "label.pcl", null);
 }

}
