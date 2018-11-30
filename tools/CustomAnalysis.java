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

package tools;

import java.util.*;
import java.io.*;
import java.awt.Color;
import java.text.MessageFormat;

import tools.microarray.*;
import tools.microarray.StepMiner.*;
import tools.microarray.FileReader.*;
import tools.microarray.FileWriter.*;

import tools.algorithm.Bimodal;

import tools.graphs.*;
import tools.io.*;

public class CustomAnalysis {

  public static String formatString(String format, double x) {
    MessageFormat mf = new MessageFormat(
        "{0,number," + format + "}");
    Object[] objs = {new Double(x)};
    return mf.format(objs);
  }

  public static double getCorrelation(Double[] v1, Double[] v2) {
    double sum_xy = 0, sum_x = 0, sum_y = 0, sum_sqx = 0, sum_sqy = 0;
    int count = 0;
    double res =0;
    int length = v1.length;
    if (length > v2.length) {
      length = v2.length;
    }
    for (int i =0; i <length; i++) {
      Double x = v1[i];
      Double y = v2[i];
      if (x != null && y != null) {
        count ++;
        sum_xy += x.doubleValue() * y.doubleValue();
        sum_x += x.doubleValue();
        sum_y += y.doubleValue();
        sum_sqx += x.doubleValue() * x.doubleValue();
        sum_sqy += y.doubleValue() * y.doubleValue();
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

  public static double max(Double[] val) {
    double max = Double.MIN_VALUE;
    for (int i = 0; i < val.length; i++) {
      if (val[i] != null && val[i].doubleValue() > max) {
        max = val[i].doubleValue();
      }
    }
    return max;
  }

  public static double min(Double[] val) {
    double min = Double.MAX_VALUE;
    for (int i = 0; i < val.length; i++) {
      if (val[i] != null && val[i].doubleValue() < min) {
        min = val[i].doubleValue();
      }
    }
    return min;
  }

  public static int getCount(Data data, double thr, int[] perm) {
    int count = 0;
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      GeneData gene = data.getGeneData(i);
      int numArrays = data.getNumArrays()/2;
      Double[] x1 = gene.getVector(data.getNumArrayHeader(),
          data.getNumArrayHeader()+numArrays-1);
      Double[] x2 = gene.getVector(data.getNumArrayHeader()+numArrays,
          data.getNumArrayHeader()+2 * numArrays-1);
      if ( (max(x1) - min(x1)) >= 1 || (max(x2) - min(x2)) >= 1 ) {
        x2 = Permutation.permute(x2, perm);
        double corr = getCorrelation(x1, x2);
        if (corr > thr) {
          count ++;
        }
      }
    }
    return count;
  }

  public static Data significantData(Data data, double thr) {
    int count = getCount(data, thr, Permutation.getNullPermutation(8));
    int[] order = new int[count];
    int index = 0;
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      GeneData gene = data.getGeneData(i);
      int numArrays = data.getNumArrays()/2;
      Double[] x1 = gene.getVector(data.getNumArrayHeader(),
          data.getNumArrayHeader()+numArrays-1);
      Double[] x2 = gene.getVector(data.getNumArrayHeader()+numArrays,
          data.getNumArrayHeader()+2 * numArrays-1);
      if ( (max(x1) - min(x1)) >= 1 || (max(x2) - min(x2)) >= 1 ) {
        double corr = getCorrelation(x1, x2);
        if (corr > thr) {
          order[index++] = i;
        }
      }
    }
    Data res = (Data) data.clone();
    res.restrictGenes(order);
    return res;
  }

  public static int getTotalCount(Data data) {
    int count = 0;
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      GeneData gene = data.getGeneData(i);
      int numArrays = data.getNumArrays()/2;
      Double[] x1 = gene.getVector(data.getNumArrayHeader(),
          data.getNumArrayHeader()+numArrays-1);
      Double[] x2 = gene.getVector(data.getNumArrayHeader()+numArrays,
          data.getNumArrayHeader()+2 * numArrays-1);
      if ( (max(x1) - min(x1)) >= 1 || (max(x2) - min(x2)) >= 1 ) {
        count ++;
      }
    }
    return count;
  }

  public static void countTotalFoldChange(LinkedList<String> list) throws Exception {
    String file1 = list.removeFirst();
    Data data1 = PCLFileReader.readFile(file1);
    Data data = Data.mergeData(data1, data1, true, false);
    System.out.println("Total = " + getTotalCount(data));
  }

  public static void corrAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    double thr = Double.parseDouble(list.removeFirst());
    String file1 = list.removeFirst();
    String file2 = list.removeFirst();
    Data data1 = PCLFileReader.readFile(file1);
    Data data2 = PCLFileReader.readFile(file2);
    Data data = Data.mergeData(data1, data2, true, false);
    Permutation p = new Permutation(100);
    int[] perm = p.getNullPermutation(data2.getNumArrays());
    int count = getCount(data, thr, perm);
    int sum_count = 0;
    for (int i =0; i < 100; i++) {
      perm = p.getRandomPermutation(data2.getNumArrays());
      sum_count += getCount(data, thr, perm);
      System.out.print(".");
      System.out.flush();
    }
    double mean_count = sum_count/100;
    System.out.println("\nActual = " + count);
    System.out.println("Mean = " + mean_count);
    System.out.println("Fdr = " + mean_count/count);
    System.out.println("Total = " + getTotalCount(data));
    Data sig = significantData(data, thr);
    PCLFileWriter.writeFile(sig, ofile);
  }

  public static void intersectionAnalysis(LinkedList<String> list) throws Exception {
    Data[] arr = new Data[list.size()];
    Iterator<String> itr = list.iterator();
    int index = 0;
    while (itr.hasNext()) {
      String file = (String) itr.next();
      arr[index++] = PCLFileReader.readFile(file);
    }
    for (int i =0; i < arr.length; i++) {
      Data current = arr[i];
      for (int j =0; j < arr.length; j++) {
        Data next = arr[j];
        Data intersect = Data.mergeData(current, next, true, false);
        System.out.print(intersect.getNumGenes() + "\t");
      }
      System.out.println("Total = "+ current.getNumGenes());
    }
  }

  public static void intersectionGenesAnalysis(LinkedList<String> list) throws Exception {
    Data[] arr = new Data[list.size()];
    Iterator<String> itr = list.iterator();
    int index = 0;
    while (itr.hasNext()) {
      String file = (String) itr.next();
      arr[index] = PCLFileReader.readFile(file);
      for (int i =0; i < arr[index].getNumRows(); i++) {
        GeneData gene = arr[index].getGeneData(i);
        gene.toUpperCase(0, 0);
      }
      index++;
    }
    for (int i =0; i < arr.length; i++) {
      Data current = arr[i];
      for (int j =0; j < arr.length; j++) {
        Data next = arr[j];
        Data intersect = Data.mergeData(current, next, true, false);
        System.out.print(intersect.getNumGenes() + "\t");
      }
      System.out.println("Total = "+ current.getNumGenes());
    }
  }

  public static void intersectStepAnalysis(LinkedList<String> list) throws Exception {
    double pvalue = Double.parseDouble(list.removeFirst());
    Data[] arr = new Data[list.size() * 2];
    Iterator<String> itr = list.iterator();
    int index = 0;
    while (itr.hasNext()) {
      String file = (String) itr.next();
      Data data = PCLFileReader.readFile(file);
      data.convertDoubles();
      StepMiner sm = new StepMiner(data);
      sm.setOneStepAnalysis();
      sm.setStepCentering(false);
      sm.setPvalueThr(pvalue);
      sm.performAnalysis();
      arr[index++] = sm.getStepOrderedData("Up");
      arr[index++] = sm.getStepOrderedData("Down");
    }
    for (int i =0; i < arr.length; i++) {
      Data current = arr[i];
      for (int j =0; j < arr.length; j++) {
        Data next = arr[j];
        Data intersect = Data.mergeData(current, next, true, false);
        System.out.print(intersect.getNumGenes() + "\t");
      }
      System.out.println("Total = "+ current.getNumGenes());
    }
  }

  public static void intersectStepAllAnalysis(LinkedList<String> list) throws Exception {
    double pvalue = Double.parseDouble(list.removeFirst());
    String ofile = list.removeFirst();
    Data[] arr = new Data[list.size() * 2];
    Iterator<String> itr = list.iterator();
    int index = 0;
    while (itr.hasNext()) {
      String file = (String) itr.next();
      Data data = PCLFileReader.readFile(file);
      data.convertDoubles();
      StepMiner sm = new StepMiner(data);
      sm.setOneStepAnalysis();
      sm.setStepCentering(false);
      sm.setPvalueThr(pvalue);
      sm.performAnalysis();
      arr[index++] = sm.getStepOrderedData("Up");
      arr[index++] = sm.getStepOrderedData("Down");
    }
    Data up = null;
    for (int i =0; i < arr.length; i+=2) {
      Data current = arr[i];
      if (up == null) {
        up = current;
      }
      else {
        Data intersect = Data.mergeData(up, current, true, false);
        up = intersect;
      }
      System.out.print(current.getNumGenes() + "\t");
    }
    System.out.println();
    Data down = null;
    for (int i =1; i < arr.length; i+=2) {
      Data current = arr[i];
      if (down == null) {
        down = current;
      }
      else {
        Data intersect = Data.mergeData(down, current, true, false);
        down = intersect;
      }
      System.out.print(current.getNumGenes() + "\t");
    }
    System.out.println();
    System.out.println("Total Up = "+ up.getNumGenes());
    System.out.println("Total Down = "+ down.getNumGenes());
    Data combined = Data.concatData(up, down);
    PCLFileWriter.writeFile(combined, ofile);
  }

  public static void intersectStepGOAnalysis(LinkedList<String> list) throws Exception {
    double spvalue = Double.parseDouble(list.removeFirst());
    String outFile = list.removeFirst();
    String onnFile = list.removeFirst();
    String annFile = list.removeFirst();
    String org = list.removeFirst();
    double pvalue = Double.parseDouble(list.removeFirst());

    HashSet<String> allNodes = new HashSet<String>();
    GeneSet set = new GeneSet();

    Data[] arr = new Data[list.size() * 2];
    Iterator<String> itr = list.iterator();
    int index = 0;
    while (itr.hasNext()) {
      String file = (String) itr.next();
      Data data = PCLFileReader.readFile(file);
      data.convertDoubles();

      GeneNameScheme ns = new GeneNameScheme(1, org, "\\|\\|", 1);
      ns.setAnnotationFile(annFile);
      ns.setOntologyFile(onnFile);
      ns.setNumMissingPoints(0);
      data.setGeneNameScheme(ns);

      StepMiner sm = new StepMiner(data);
      sm.setOneStepAnalysis();
      sm.setStepCentering(false);
      sm.setPvalueThr(spvalue);
      sm.performAnalysis();

      arr[index] = sm.getStepOrderedData("Up");
      String[] genes = arr[index].getGenes(arr[index].getNullOrder());
      String tag = file.split("\\.")[0] + "-Up";
      arr[index].setName(tag);
      for (int i=0; i < genes.length; i++) {
        set.add(tag, genes[i]);
        allNodes.add(genes[i]);
      }
      index ++;
      arr[index] = sm.getStepOrderedData("Down");
      genes = arr[index].getGenes(arr[index].getNullOrder());
      tag = file.split("\\.")[0] + "-Down";
      arr[index].setName(tag);
      for (int i=0; i < genes.length; i++) {
        set.add(tag, genes[i]);
        allNodes.add(genes[i]);
      }
      index ++;
    }
    for (int i =0; i < arr.length; i++) {
      Data current = arr[i];
      for (int j =0; j < arr.length; j++) {
        Data next = arr[j];
        Data intersect = Data.mergeData(current, next, true, false);
        String tag = current.getName() + "-" + next.getName();
        String[] genes = intersect.getGenes(intersect.getNullOrder());
        for (int k=0; k < genes.length; k++) {
          set.add(tag, genes[k]);
        }
        System.out.print(intersect.getNumGenes() + "\t");
      }
      System.out.println("Total = "+ current.getNumGenes());
    }

    Vector<String> allNodesV = new Vector<String>(allNodes);
    GeneSetAnalysis ana = new GeneSetAnalysis(set);
    ana.performGOAnalysis(outFile, onnFile, annFile, org, pvalue, allNodesV);
  }

  public static void intersectStepGSAnalysis(LinkedList<String> list) throws Exception {
    double spvalue = Double.parseDouble(list.removeFirst());
    String outFile = list.removeFirst();
    String setFile = list.removeFirst();
    String org = list.removeFirst();
    double pvalue = Double.parseDouble(list.removeFirst());

    HashSet<String> allNodes = new HashSet<String>();
    GeneSet set = new GeneSet();

    Data[] arr = new Data[list.size() * 2];
    Iterator<String> itr = list.iterator();
    int index = 0;
    while (itr.hasNext()) {
      String file = (String) itr.next();
      Data data = PCLFileReader.readFile(file);
      data.convertDoubles();

      GeneNameScheme ns = new GeneNameScheme(1, org, "\\|\\|", 1);
      ns.setNumMissingPoints(0);
      data.setGeneNameScheme(ns);

      StepMiner sm = new StepMiner(data);
      sm.setOneStepAnalysis();
      sm.setStepCentering(false);
      sm.setPvalueThr(spvalue);
      sm.performAnalysis();

      arr[index] = sm.getStepOrderedData("Up");
      String[] genes = arr[index].getGenes(arr[index].getNullOrder());
      String tag = file.split("\\.")[0] + "-Up";
      arr[index].setName(tag);
      for (int i=0; i < genes.length; i++) {
        set.add(tag, genes[i]);
        allNodes.add(genes[i]);
      }
      index ++;
      arr[index] = sm.getStepOrderedData("Down");
      genes = arr[index].getGenes(arr[index].getNullOrder());
      tag = file.split("\\.")[0] + "-Down";
      arr[index].setName(tag);
      for (int i=0; i < genes.length; i++) {
        set.add(tag, genes[i]);
        allNodes.add(genes[i]);
      }
      index ++;
    }
    for (int i =0; i < arr.length; i++) {
      Data current = arr[i];
      for (int j =0; j < arr.length; j++) {
        Data next = arr[j];
        Data intersect = Data.mergeData(current, next, true, false);
        String tag = current.getName() + "-" + next.getName();
        String[] genes = intersect.getGenes(intersect.getNullOrder());
        for (int k=0; k < genes.length; k++) {
          set.add(tag, genes[k]);
        }
        System.out.print(intersect.getNumGenes() + "\t");
      }
      System.out.println("Total = "+ current.getNumGenes());
    }

    Vector<String> allNodesV = new Vector<String>(allNodes);
    GeneSet b = GeneSet.readFile(setFile);
    GeneSetAnalysis ana = new GeneSetAnalysis(set, b);
    ana.performAnalysis(outFile, org, pvalue);
  }

  public static int monotonicCount(Data data, double thr, int[] perm,
      Comparator<Double> comp) {
    double[] pvalues = new double[data.getNumGenes()];
    int res = 0;
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      GeneData gene = data.getGeneData(i);
      int numArrays = data.getNumArrays();
      Double[] x1 = gene.getVector(data.getNumArrayHeader(),
          data.getNumArrayHeader()+numArrays-1);
      Double[] x2 = gene.getVector(data.getNumArrayHeader(),
          data.getNumArrayHeader()+numArrays-1);
      Arrays.sort(x2, comp);
      x1 = Permutation.permute(x1, perm);
      int count = 0;
      for (int j = 0; j < x2.length; j++) {
        if (x1[j].equals(x2[j])) {
          count++;
        }
      }
      double corr = getCorrelation(x1, x2);
      double statistic = (numArrays-count-2) * corr*corr/(1- corr*corr);
      double pvalue = 1 - Utils.Fisher(statistic, 1, (numArrays-count-2));
      pvalues[i-data.getNumGeneHeader()] = pvalue;
      if (corr >= 0 && pvalue < thr) {
        res++;
      }
    }
    return res;
  }

  public static Data monotonicData(Data data, double thr,
      Comparator<Double> comp) {
    double[] pvalues = new double[data.getNumGenes()];
    Vector<Integer> ord = new Vector<Integer>();
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      GeneData gene = data.getGeneData(i);
      int numArrays = data.getNumArrays();
      Double[] x1 = gene.getVector(data.getNumArrayHeader(),
          data.getNumArrayHeader()+numArrays-1);
      Double[] x2 = gene.getVector(data.getNumArrayHeader(),
          data.getNumArrayHeader()+numArrays-1);
      Arrays.sort(x2, comp);
      int count = 0;
      for (int j = 0; j < x2.length; j++) {
        if (x1[j].equals(x2[j])) {
          count++;
        }
      }
      double corr = getCorrelation(x1, x2);
      double statistic = (numArrays-count-2) * corr*corr/(1- corr*corr);
      double pvalue = 1 - Utils.Fisher(statistic, 1, (numArrays-count-2));
      pvalues[i-data.getNumGeneHeader()] = pvalue;
      if (corr >= 0 && pvalue < thr) {
        ord.add(new Integer(i));
      }
    }
    int[] order = new int[ord.size()];
    for (int i =0; i < ord.size(); i++) {
      Integer val = ord.get(i);
      order[i] = val.intValue();
    }
    Data res = (Data) data.clone();
    res.restrictGenes(order);
    return res;
  }

  public static void monotonicAnalysis(LinkedList<String> list) throws Exception {
    String choice = list.removeFirst();
    String ofile = list.removeFirst();
    double thr = Double.parseDouble(list.removeFirst());
    String file = list.removeFirst();
    Data data = PCLFileReader.readFile(file);
    class DesComparator implements Comparator<Double> {
      boolean asc_;
      public DesComparator(boolean ascending) {
        asc_ = ascending;
      }
      public int compare(Double c1, Double c2) {
        int res = 1;
        if (asc_) {
          res = -1;
        }
        if (c1.doubleValue() < c2.doubleValue()) {
          return res;
        }
        return -res;
      }
    };
    DesComparator comp = new DesComparator(true);
    if (choice.equals("down")) {
      comp = new DesComparator(false);
    }
    Permutation p = new Permutation(100);
    int[] perm = p.getNullPermutation(data.getNumArrays());
    int count = monotonicCount(data, thr, perm, comp);
    int sum_count = 0;
/*
    for (int i =0; i < 100; i++) {
      perm = p.getRandomPermutation(data.getNumArrays());
      sum_count += monotonicCount(data, thr, perm, comp);
      System.out.print(".");
      System.out.flush();
    }
*/
    double mean_count = sum_count/100;
    System.out.println("\nActual = " + count);
    System.out.println("Mean = " + mean_count);
    System.out.println("Fdr = " + mean_count/count);
    System.out.println("Total = " + data.getNumGenes());
    Data sig = monotonicData(data, thr, comp);
    PCLFileWriter.writeFile(sig, ofile);
  }

  public static int countCorrStep(Data data1, Data data2, Data data3, 
      double spvalue, double thr, int[] perm1, int[] perm2, 
      int[] perm3) throws Exception {
    data1 = Data.selectArraysFromData(data1, perm1);
    data2 = Data.selectArraysFromData(data2, perm2);
    data3 = Data.selectArraysFromData(data3, perm3);
    Data data = Data.mergeData(data1, data2, true, false);
    data = Data.mergeData(data, data3, true, false);
    int h = data.getNumArrayHeader();
    int n1 = data1.getNumArrays();
    int n2 = data2.getNumArrays();
    int n3 = data3.getNumArrays();
    data.convertDoubles();
    data.setRange(h + ":" + (h+n1-1));
    StepMiner sm = new StepMiner(data);
    sm.setOneStepAnalysis();
    sm.setStepCentering(false);
    sm.setPvalueThr(spvalue);
    sm.performAnalysis();
    Data smdata1 = sm.getStepOrderedData("Step");
    data.setRange((h+n1) + ":" + (h+n1+n2-1));
    sm.performAnalysis();
    Data smdata2 = sm.getStepOrderedData("Step");
    data.setRange((h+n1+n2) + ":" + (h+n1+n2+n3-1));
    sm.performAnalysis();
    Data smdata3 = sm.getStepOrderedData("Step");
    Data smdata = Data.mergeData(smdata1, smdata2, false, false);
    smdata = Data.mergeData(smdata, smdata3, false, false);
    data = Data.selectOrder(data, smdata);
    Vector<Integer> ord = new Vector<Integer>();
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      GeneData gene = data.getGeneData(i);
      int numArrays = data.getNumArrays();
      Double[] x1 = gene.getVector(h, h+n1-1);
      Double[] x2 = gene.getVector(h+n1, h+n1+n2-1);
      Double[] x3 = gene.getVector(h+n1+n2, h+n1+n2+n3-1);
      double corr1 = getCorrelation(x1, x2);
      double corr2 = getCorrelation(x2, x3);
      double corr3 = getCorrelation(x1, x3);
      double statistic = (corr1+corr2+corr3)/3;
      if (statistic >= thr) {
        ord.add(new Integer(i));
      }
    }
    return ord.size();
  }

  public static Data corrStepData(Data data1, Data data2, Data data3, 
      double spvalue, double thr) throws Exception {
    Data data = Data.mergeData(data1, data2, true, false);
    data = Data.mergeData(data, data3, true, false);
    int h = data.getNumArrayHeader();
    int n1 = data1.getNumArrays();
    int n2 = data2.getNumArrays();
    int n3 = data3.getNumArrays();
    data.convertDoubles();
    data.setRange(h + ":" + (h+n1-1));
    StepMiner sm = new StepMiner(data);
    sm.setOneStepAnalysis();
    sm.setStepCentering(false);
    sm.setPvalueThr(spvalue);
    sm.performAnalysis();
    Data smdata1 = sm.getStepOrderedData("Step");
    data.setRange((h+n1) + ":" + (h+n1+n2-1));
    sm.performAnalysis();
    Data smdata2 = sm.getStepOrderedData("Step");
    data.setRange((h+n1+n2) + ":" + (h+n1+n2+n3-1));
    sm.performAnalysis();
    Data smdata3 = sm.getStepOrderedData("Step");
    Data smdata = Data.mergeData(smdata1, smdata2, false, false);
    smdata = Data.mergeData(smdata, smdata3, false, false);
    data = Data.selectOrder(data, smdata);
    System.out.println("Total Union : " + data.getNumGenes());
    Vector<Integer> ord = new Vector<Integer>();
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      GeneData gene = data.getGeneData(i);
      int numArrays = data.getNumArrays();
      Double[] x1 = gene.getVector(h, h+n1-1);
      Double[] x2 = gene.getVector(h+n1, h+n1+n2-1);
      Double[] x3 = gene.getVector(h+n1+n2, h+n1+n2+n3-1);
      double corr1 = getCorrelation(x1, x2);
      double corr2 = getCorrelation(x2, x3);
      double corr3 = getCorrelation(x1, x3);
      double statistic = (corr1+corr2+corr3)/3;
      if (statistic >= thr) {
        ord.add(new Integer(i));
      }
    }
    int[] order = new int[ord.size()];
    for (int i =0; i < ord.size(); i++) {
      Integer val = ord.get(i);
      order[i] = val.intValue();
    }
    Data res = (Data) data.clone();
    res.restrictGenes(order);
    Data rest = Data.diffData(data, res);
    return Data.concatData(res, rest);
  }

  public static void corrStepAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    double spvalue = Double.parseDouble(list.removeFirst());
    double thr = Double.parseDouble(list.removeFirst());
    String file1 = list.removeFirst();
    String file2 = list.removeFirst();
    String file3 = list.removeFirst();
    Data data1 = PCLFileReader.readFile(file1);
    Data data2 = PCLFileReader.readFile(file2);
    Data data3 = PCLFileReader.readFile(file3);
    Permutation p = new Permutation(100);
    int[] perm1 = p.getNullPermutation(data1.getNumArrays());
    int[] perm2 = p.getNullPermutation(data2.getNumArrays());
    int[] perm3 = p.getNullPermutation(data3.getNumArrays());
    int count = countCorrStep(data1, data2, data3, spvalue, thr, perm1, perm2, perm3);
    int sum_count = 0;
    for (int i =0; i < 100; i++) {
      perm1 = p.getRandomPermutation(data1.getNumArrays());
      perm2 = p.getRandomPermutation(data2.getNumArrays());
      perm3 = p.getRandomPermutation(data3.getNumArrays());
      sum_count += countCorrStep(data1, data2, data3, spvalue, thr, perm1, perm2, perm3);
      System.out.print(".");
      System.out.flush();
    }
    double mean_count = sum_count/100;
    System.out.println("\nActual = " + count);
    System.out.println("Mean = " + mean_count);
    System.out.println("Fdr = " + mean_count/count);
    Data sig = corrStepData(data1, data2, data3, spvalue, thr);
    PCLFileWriter.writeFile(sig, ofile);
  }

  public static void aracneAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    int num = Integer.parseInt(list.removeFirst());
    String file1 = list.removeFirst();
    String file2 = list.removeFirst();
    Data data1 = PCLFileReader.readFile(file1);
    Data data2 = PCLFileReader.readFile(file2);
    int[] order = new int[2*num];
    for (int i=0; i < num; i++) {
      order[i] = i;
    }
    for (int i=0; i < num; i++) {
      order[num+i] = data2.getNumGenes() - num + i;
    }
    data2 = Data.selectGenesFromData(data2, order);
    Data col1 = Data.selectColumnFromData(data2, 0);
    Data col2 = Data.selectColumnFromData(data2, 2);
    Data col3 = Data.selectColumnFromData(data2, 4);
    col1 = Data.selectOrder(data1, col1);
    col2 = Data.selectOrder(data1, col2);
    col1 = (Data) col1.clone();
    col2 = (Data) col2.clone();
    col1.reduceLog();
    col2.reduceLog();
    PCLFileWriter.writeFile(col1, "label1.pcl");
    PCLFileWriter.writeFile(col2, "label2.pcl");
    GeneNameScheme ns = new GeneNameScheme(1, "Hs", ":", 0);
    ns.setNumMissingPoints(0);
    col1.setGeneNameScheme(ns);
    col2.setGeneNameScheme(ns);
    PSPlot plotter = new PSPlot(ofile);
    plotter.open();
    plotter.array(4 /* Rows */, 2 /* columns */);
    for (int i = col1.getNumGeneHeader(); i < col1.getNumRows(); i++) {
      GeneData gene1 = col1.getGeneData(i);
      GeneData gene2 = col2.getGeneData(i);
      int start = col1.getNumArrayHeader();
      int end = col1.getNumColumns()-1;
      Double[] data = gene1.getVector(start, end);
      Double[] time = gene2.getVector(start, end);
      String name1 = col1.getGenesAt(i);
      String name2 = col2.getGenesAt(i);
      String name3 = (String) col3.getGeneData(i).getDataAt(0);
      String id1 = (String) col1.getGeneData(i).getDataAt(0);
      String id2 = (String) col2.getGeneData(i).getDataAt(0);
      if (!name1.equals(name2)) {
        plotter.plot(time, data);
        plotter.xlabel(name2);
        plotter.ylabel(name1);
        plotter.title(id1 + "," + id2 + " MI=" + name3);
      }
    }
    plotter.close();
  }

  public static void writeBimodalGenes(boolean reducelog, LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    double thr = Double.parseDouble(list.removeFirst());
    Data data1 = PCLFileReader.readFile(file1);

    if (reducelog) {
      data1.reduceLog();
    }

    Vector<Integer> order = new Vector<Integer>();
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      // System.out.println("[ " + i + " ]---------------------->");
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v= data1.getGeneData(i).getVector(start, end);
      if (Bimodal.checkEntropy(v, thr)) {
        order.add(new Integer(i-data1.getNumGeneHeader()));
      }
    }
    int[] ord = new int[order.size()];
    for (int i=0; i < order.size(); i++) {
        Integer val = order.get(i);
        ord[i] = val.intValue();
    }
    System.out.println("Number of Bimodal genes = " + order.size());
    Data res = Data.selectGenesFromData(data1, ord);
    PCLFileWriter.writeFile(res, ofile);

  }

  public static void writeThreshold(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    if (list.size() > 0) {
        double gap = Double.parseDouble(list.removeFirst());
        Bimodal.GAP_LENGTH = gap;
    }
    Data data1 = PCLFileReader.readFile(file1);

    Bimodal[] b = new Bimodal[data1.getNumGenes()];
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      // System.out.println("[ " + i + " ]---------------------->");
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v= data1.getGeneData(i).getVector(start, end);
      b[i - data1.getNumGeneHeader()] = new Bimodal(v);
    }

    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      int i1 = i - data1.getNumGeneHeader();
      double thr1 = b[i1].getThreshold();
      double lowthr1 = b[i1].getLowThreshold();
      double highthr1 = b[i1].getHighThreshold();
      double st = b[i1].getStatistic();
      out.write(i + "\t" + thr1 + "\t" + st + "\t" + lowthr1 + "\t" + highthr1 + "\n");
    }
    out.close();

  }

  /*
   *  Streamline the PCL file reading.
   *        Memory efficient version of the above function
   */
  public static void writeThreshold1(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    if (list.size() > 0) {
        double gap = Double.parseDouble(list.removeFirst());
        Bimodal.GAP_LENGTH = gap;
    }
    PCLFileReader data1 = new PCLFileReader(file1);
    data1.begin();
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    while (data1.hasNext()) {
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      long lineno = data1.getLineNumber();
      GeneData gene = data1.getData();
      // gene.print();
      if (gene == null) {
        break;
      }
      Double[] v= gene.getVector(start, end);
      Bimodal b = new Bimodal(v);
      double thr1 = b.getThreshold();
      double lowthr1 = b.getLowThreshold();
      double highthr1 = b.getHighThreshold();
      double st = b.getStatistic();
      out.write(lineno + "\t" + thr1 + "\t" + st + "\t" + lowthr1 + "\t" + highthr1 + "\n");
    }
    out.close();

  }

  /*
   *  Streamline the PCL file reading.
   *        Write the Bit vector
   */
  public static void writeThresholdBv(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    if (list.size() > 0) {
        double gap = Double.parseDouble(list.removeFirst());
        Bimodal.GAP_LENGTH = gap;
    }
    PCLFileReader data1 = new PCLFileReader(file1);
    data1.begin();
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    GeneData gene = data1.getHeader();
    out.write(gene.getDataAt(0) + "\t");
    out.write(gene.getDataAt(1) + "\t");
    out.write("Bit Vector\n");
    while (data1.hasNext()) {
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      long lineno = data1.getLineNumber();
      gene = data1.getData();
      if (gene == null) {
        break;
      }
      Double[] v= gene.getVector(start, end);
      Bimodal b = new Bimodal(v);
      double thr1 = b.getThreshold();
      double lowthr1 = b.getLowThreshold();
      double highthr1 = b.getHighThreshold();
      out.write(gene.getDataAt(0) + "\t");
      out.write(gene.getDataAt(1) + "\t");
      v= gene.getVector(start, end);
      for (int i =0; i < v.length; i++) {
        if (v[i] == null) {
          out.write(" ");
        }
        else if (v[i].doubleValue() <= lowthr1) {
          out.write("0");
        }
        else if (v[i].doubleValue() >= highthr1) {
          out.write("2");
        }
        else {
          out.write("1");
        }
      }
      out.write("\n");
    }
    out.close();

  }

  public static void writeThresholdBv1(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    String thr_file = list.removeFirst();
    NetworkInfo info = new NetworkInfo();
    info.setThrFile(thr_file);
    info.readThresholds();
    PCLFileReader data1 = new PCLFileReader(file1);
    data1.begin();
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    GeneData gene = data1.getHeader();
    out.write(gene.getDataAt(0) + "\t");
    out.write(gene.getDataAt(1) + "\t");
    out.write("Bit Vector\n");
    while (data1.hasNext()) {
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      long lineno = data1.getLineNumber();
      gene = data1.getData();
      if (gene == null) {
        break;
      }
      Double[] v= gene.getVector(start, end);
      double thr1 = info.getThresholdByIndex(lineno).doubleValue();
      double lowthr1 = info.getLowerThresholdByIndex(lineno).doubleValue();
      double highthr1 = info.getUpperThresholdByIndex(lineno).doubleValue();
      out.write(gene.getDataAt(0) + "\t");
      out.write(gene.getDataAt(1) + "\t");
      v= gene.getVector(start, end);
      for (int i =0; i < v.length; i++) {
        if (v[i] == null) {
          out.write(" ");
        }
        else if (v[i].doubleValue() <= lowthr1) {
          out.write("0");
        }
        else if (v[i].doubleValue() >= highthr1) {
          out.write("2");
        }
        else {
          out.write("1");
        }
      }
      out.write("\n");
    }
    out.close();

  }

  public static double STAT_THRESHOLD = 5;
  public static double THRESHOLD = 0.10;
  public static int SIZE_CUTOFF = 30;

  public static void writePairs(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    String geneid = null;
    if (list.size() > 0) {
       geneid = list.removeFirst();
    }
    Data data1 = PCLFileReader.readFile(file1);

    Bimodal[] b = new Bimodal[data1.getNumGenes()];
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      // System.out.println("[ " + i + " ]---------------------->");
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v= data1.getGeneData(i).getVector(start, end);
      b[i - data1.getNumGeneHeader()] = new Bimodal(v);
    }

    int count1 = 0;
    int count2 = 0;
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      GeneData gene = data1.getGeneData(i);
      String id = (String) gene.getDataAt(0);
      if (geneid != null && !id.equals(geneid)) {
        continue;
      }
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v1 = gene.getVector(start, end);
      int i1 = i - data1.getNumGeneHeader();
      double thr1 = b[i1].getThreshold();
      System.out.println(i);
      boolean found = false;
      for (int j = data1.getNumGeneHeader(); j < data1.getNumRows(); j++) {
        if ( i == j ) continue;
        int j1 = j - data1.getNumGeneHeader();
        Double[] v2 = data1.getGeneData(j).getVector(start, end);
        double thr2 = b[j1].getThreshold();
        int c0 = countThreshold(0, thr1, thr2, v1, v2);
        int c1 = countThreshold(1, thr1, thr2, v1, v2);
        int c2 = countThreshold(2, thr1, thr2, v1, v2);
        int c3 = countThreshold(3, thr1, thr2, v1, v2);
        // System.out.println(thr1 + "," + thr2 + ":" + c0 + "-" +c1 + "-" +c2 + "-" +c3);
        boolean best = false;
        double ratio = 0;
        double len = (double) v1.length;
        /*
        if ((c0 + c2) > (v1.length * 0.90)) {
            best = true;
            found = true;
            count1 ++;
            ratio = (c0 + c2) / len;
            String sr = formatString("0.##", ratio);
            out.write(i + "\t" + j + "\t+" + "\t" + sr + "\n");
        }
        if ((c1 + c3) > (v1.length * 0.90)) {
            best = true;
            found = true;
            count2 ++;
            ratio = (c1 + c3) / len;
            String sr = formatString("0.##", ratio);
            out.write(i + "\t" + j + "\t-" + "\t" + sr + "\n");
        }
        int count = 0;
        int type = 0;
        if (!best && c0 < (v1.length * 0.10)){count++;type = 0;ratio=c0/len;}
        if (!best && c1 < (v1.length * 0.10)){count++;type = 1;ratio=c1/len;}
        if (!best && c2 < (v1.length * 0.10)){count++;type = 2;ratio=c2/len;}
        if (!best && c3 < (v1.length * 0.10)){count++;type = 3;ratio=c3/len;}
        if (!best && count == 1) {
          ratio = 1 - ratio;
          String sr = formatString("0.##", ratio);
          out.write(i + "\t" + j + "\t" + type + "\t" + sr + "\n");
        }
        // break;
        */
        double p0 = (c0/(c0+c1+1.0) + c0/(c0+c3+1.0))/2;
        double p1 = (c1/(c1+c2+1.0) + c1/(c1+c0+1.0))/2;
        double p2 = (c2/(c2+c3+1.0) + c2/(c2+c1+1.0))/2;
        double p3 = (c3/(c3+c0+1.0) + c3/(c3+c2+1.0))/2;
        double thr = THRESHOLD;
        if (p0 < thr) {
          String sr = formatString("0.#####", p0);
          out.write(i + "\t" + j + "\t" + 0 + "\t" + sr + "\n");
        }
        if (p1 < thr) {
          String sr = formatString("0.#####", p1);
          out.write(i + "\t" + j + "\t" + 1 + "\t" + sr + "\n");
        }
        if (p2 < thr) {
          String sr = formatString("0.#####", p2);
          out.write(i + "\t" + j + "\t" + 2 + "\t" + sr + "\n");
        }
        if (p3 < thr) {
          String sr = formatString("0.#####", p3);
          out.write(i + "\t" + j + "\t" + 3 + "\t" + sr + "\n");
        }
      }
/*
      if (!found) {
        out.write(i + "\t" + i + "\t*\t1.00\n");
      }
*/
      if (geneid != null && id.equals(geneid)) {
        break;
      }
      // break;
    }
    out.close();
    System.out.println(" (0x2) = " + count1 + " (1x3) = " + count2 );

  }

  /*
   * With ignoring the intermediate values
   */
  public static void writePairsNew1(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    String geneid = null;
    if (list.size() > 0) {
       geneid = list.removeFirst();
    }
    Data data1 = PCLFileReader.readFile(file1);

    Bimodal[] b = new Bimodal[data1.getNumGenes()];
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      // System.out.println("[ " + i + " ]---------------------->");
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      GeneData gene = data1.getGeneData(i);
      Double[] v= gene.getVector(start, end);
      int i1 = i - data1.getNumGeneHeader();
      b[i1] = new Bimodal(v);
      double lowthr1 = b[i1].getLowThreshold();
      double highthr1 = b[i1].getHighThreshold();
      v= gene.getVector(start, end);
      for (int j = 0; j < v.length; j++) {
        if (v[j] != null && v[j].doubleValue() > lowthr1 && 
            v[j].doubleValue() < highthr1) {
            gene.setDataAt(start + j, null);
        }
      }
    }

    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      GeneData gene = data1.getGeneData(i);
      String id = (String) gene.getDataAt(0);
      if (geneid != null && !id.equals(geneid)) {
        continue;
      }
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v1 = gene.getVector(start, end);
      int i1 = i - data1.getNumGeneHeader();
      double thr1 = b[i1].getThreshold();
      double lowthr1 = b[i1].getLowThreshold();
      double highthr1 = b[i1].getHighThreshold();
      System.out.println(i);
      boolean found = false;
      for (int j = data1.getNumGeneHeader(); j < data1.getNumRows(); j++) {
        if ( i == j ) continue;
        int j1 = j - data1.getNumGeneHeader();
        Double[] v2 = data1.getGeneData(j).getVector(start, end);
        double thr2 = b[j1].getThreshold();
        double lowthr2 = b[j1].getLowThreshold();
        double highthr2 = b[j1].getHighThreshold();
        double p0 = getImpStat(0, lowthr1, lowthr2, v1, v2);
        double p1 = getImpStat(1, lowthr1, highthr2, v1, v2);
        double p2 = getImpStat(2, highthr1, highthr2, v1, v2);
        double p3 = getImpStat(3, highthr1, lowthr2, v1, v2);
        double thr = THRESHOLD;
        if (i == 4 && j == 6) {
          String sr = formatString("0.#####", p0);
          out.write(i + "\t" + j + "\t" + 0 + "\t" + sr + "\n");
          sr = formatString("0.#####", p1);
          out.write(i + "\t" + j + "\t" + 0 + "\t" + sr + "\n");
          sr = formatString("0.#####", p2);
          out.write(i + "\t" + j + "\t" + 0 + "\t" + sr + "\n");
          sr = formatString("0.#####", p3);
          out.write(i + "\t" + j + "\t" + 0 + "\t" + sr + "\n");
        }
        if (p0 <= thr) {
          String sr = formatString("0.#####", p0);
          out.write(i + "\t" + j + "\t" + 0 + "\t" + sr + "\n");
        }
        if (p1 <= thr) {
          String sr = formatString("0.#####", p1);
          out.write(i + "\t" + j + "\t" + 1 + "\t" + sr + "\n");
        }
        if (p2 <= thr) {
          String sr = formatString("0.#####", p2);
          out.write(i + "\t" + j + "\t" + 2 + "\t" + sr + "\n");
        }
        if (p3 <= thr) {
          String sr = formatString("0.#####", p3);
          out.write(i + "\t" + j + "\t" + 3 + "\t" + sr + "\n");
        }
      }
      if (geneid != null && id.equals(geneid)) {
        break;
      }
      // break;
    }
    out.close();
  }

  public static void printBitSet(BitSet b) {
    for (int i =0; i < b.size(); i++) {
      if (b.get(i)) {
        System.out.print("1");
      }
      else {
        System.out.print("0");
      }
    }
    System.out.println();
  }

  public static double[] getErrorProbability(
      BitSet a, BitSet a_thr, BitSet b, BitSet b_thr, int correction) { 
    double[] res = new double[4];
    res[0] = res[1] = res[2] = res[3] = 1.0;
    if (a.length() == 0 || b.length() == 0) {
      return res;
    }
    int num = a.size();
    BitSet thrBits = (BitSet) a_thr.clone();
    thrBits.and(b_thr);
    BitSet tmp = (BitSet) thrBits.clone();
    //BitSet v1 = a.get(0,  a.length()-1);
    BitSet v1 = (BitSet) a.clone();
    v1.or(b);
    tmp.andNot(v1);
    int c0 = tmp.cardinality() + correction;
    tmp = (BitSet) thrBits.clone();
    //v1 = b.get(0,  b.length()-1);
    v1 = (BitSet) b.clone();
    v1.andNot(a);
    tmp.and(v1);
    int c1 = tmp.cardinality() + correction;
    tmp = (BitSet) thrBits.clone();
    //v1 = a.get(0,  a.length()-1);
    v1 = (BitSet) a.clone();
    v1.andNot(b);
    tmp.and(v1);
    int c2 = tmp.cardinality() + correction;
    tmp = (BitSet) thrBits.clone();
    //v1 = a.get(0,  a.length()-1);
    v1 = (BitSet) a.clone();
    v1.and(b);
    tmp.and(v1);
    int c3 = tmp.cardinality() + correction;
    /*
    int total = c0 + c1 + c2 + c3;
    if (num > (2 * total)) {
      return res;
    }
    */
    //System.out.println(c0 + "\t" + c1 + "\t" + c2 + "\t" + c3);
    //System.out.println(correction);

    res[0] = 0.5 * c0 / (c0 + c1) + 0.5 * c0/(c0 + c2);
    res[1] = 0.5 * c1 / (c1 + c0) + 0.5 * c1/(c1 + c3);
    res[2] = 0.5 * c2 / (c2 + c0) + 0.5 * c2/(c2 + c3);
    res[3] = 0.5 * c3 / (c3 + c1) + 0.5 * c3/(c3 + c2);
    return res;
  }

  public static double[] getErrorStats(
      BitSet a, BitSet a_thr, BitSet b, BitSet b_thr, int correction, int debug) { 
    double[] res = new double[4];
    res[0] = res[1] = res[2] = res[3] = 1.0;
    if (a.length() == 0 || b.length() == 0) {
      return res;
    }
    int num = a.size();
    BitSet thrBits = (BitSet) a_thr.clone();
    thrBits.and(b_thr);
    BitSet tmp = (BitSet) thrBits.clone();
    //BitSet v1 = a.get(0,  a.length()-1);
    BitSet v1 = (BitSet) a.clone();
    v1.or(b);
    tmp.andNot(v1);
    int c0 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    //v1 = b.get(0,  b.length()-1);
    v1 = (BitSet) b.clone();
    v1.andNot(a);
    tmp.and(v1);
    int c1 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    //v1 = a.get(0,  a.length()-1);
    v1 = (BitSet) a.clone();
    v1.andNot(b);
    tmp.and(v1);
    int c2 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    //v1 = a.get(0,  a.length()-1);
    v1 = (BitSet) a.clone();
    v1.and(b);
    tmp.and(v1);
    int c3 = tmp.cardinality();
    int total = c0 + c1 + c2 + c3;
    /*
    if (num > (2 * total)) {
      return res;
    }
    */

    res[0] = ((c0 + c1) * (c0 + c2)/total - c0 + 1)/Math.sqrt((c0 + c1) * (c0 + c2)/total + 1);
    res[1] = ((c1 + c0) * (c1 + c3)/total - c1 + 1)/Math.sqrt((c1 + c0) * (c1 + c3)/total + 1);
    res[2] = ((c2 + c0) * (c2 + c3)/total - c2 + 1)/Math.sqrt((c2 + c0) * (c2 + c3)/total + 1);
    res[3] = ((c3 + c1) * (c3 + c2)/total - c3 + 1)/Math.sqrt((c3 + c1) * (c3 + c2)/total + 1);
    res[0] = (1/res[0]) < 0 ? 1:(1/res[0]);
    res[1] = (1/res[1]) < 0 ? 1:(1/res[1]);
    res[2] = (1/res[2]) < 0 ? 1:(1/res[2]);
    res[3] = (1/res[3]) < 0 ? 1:(1/res[3]);
    if (debug > 0) {
      System.out.println(c0 + "\t" + c1 + "\t" + c2 + "\t" + c3 + "\t" + total);
      System.out.println(correction);
      System.out.println(res[0] + "\t" + res[1] + "\t" + res[2] + "\t" + res[3]);
    }
    return res;
  }

  public static double[] getErrorProbStats(
      BitSet a, BitSet a_thr, BitSet b, BitSet b_thr, int debug) { 
    double[] res = new double[4];
    res[0] = res[1] = res[2] = res[3] = 1.0;
    if (a.length() == 0 || b.length() == 0) {
      return res;
    }
    int num = a.size();
    BitSet thrBits = (BitSet) a_thr.clone();
    thrBits.and(b_thr);
    BitSet tmp = (BitSet) thrBits.clone();
    //BitSet v1 = a.get(0,  a.length()-1);
    BitSet v1 = (BitSet) a.clone();
    v1.or(b);
    tmp.andNot(v1);
    int c0 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    //v1 = b.get(0,  b.length()-1);
    v1 = (BitSet) b.clone();
    v1.andNot(a);
    tmp.and(v1);
    int c1 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    //v1 = a.get(0,  a.length()-1);
    v1 = (BitSet) a.clone();
    v1.andNot(b);
    tmp.and(v1);
    int c2 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    //v1 = a.get(0,  a.length()-1);
    v1 = (BitSet) a.clone();
    v1.and(b);
    tmp.and(v1);
    int c3 = tmp.cardinality();
    int total = c0 + c1 + c2 + c3;

    res[0] = ((c0 + c1) * (c0 + c2)/total - c0 + 1)/Math.sqrt((c0 + c1) * (c0 + c2)/total + 1);
    res[1] = ((c1 + c0) * (c1 + c3)/total - c1 + 1)/Math.sqrt((c1 + c0) * (c1 + c3)/total + 1);
    res[2] = ((c2 + c0) * (c2 + c3)/total - c2 + 1)/Math.sqrt((c2 + c0) * (c2 + c3)/total + 1);
    res[3] = ((c3 + c1) * (c3 + c2)/total - c3 + 1)/Math.sqrt((c3 + c1) * (c3 + c2)/total + 1);
    if (debug > 0) {
      System.out.println(c0 + "\t" + c1 + "\t" + c2 + "\t" + c3 + "\t" + total);
      System.out.println(res[0] + "\t" + res[1] + "\t" + res[2] + "\t" + res[3]);
    }
    if (res[0] > STAT_THRESHOLD) {
      res[0] = 0.5 * c0 / (c0 + c1 + 1) + 0.5 * c0/(c0 + c2 + 1);
    }
    else {
      res[0] = 1;
    }
    if (res[1] > STAT_THRESHOLD) {
      res[1] = 0.5 * c1 / (c1 + c0 + 1) + 0.5 * c1/(c1 + c3 + 1);
    }
    else {
      res[1] = 1;
    }
    if (res[2] > STAT_THRESHOLD) {
      res[2] = 0.5 * c2 / (c2 + c0 + 1) + 0.5 * c2/(c2 + c3 + 1);
    }
    else {
      res[2] = 1;
    }
    if (res[3] > STAT_THRESHOLD) {
      res[3] = 0.5 * c3 / (c3 + c1 + 1) + 0.5 * c3/(c3 + c2 + 1);
    }
    else {
      res[3] = 1;
    }
    if (debug > 0) {
      System.out.println(res[0] + "\t" + res[1] + "\t" + res[2] + "\t" + res[3]);
    }
    return res;
  }

  /*
   * With Bit Vectors and ignoring the intermediate values
   */
  public static void writePairsNew(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    String geneid = null;
    int geneid_index = -1;
    if (list.size() > 0) {
       geneid = list.removeFirst();
    }
    Data data1 = PCLFileReader.readFile(file1);

    Bimodal[] b = new Bimodal[data1.getNumGenes()];
    BitSet[] val = new BitSet[data1.getNumGenes()];
    BitSet[] thr = new BitSet[data1.getNumGenes()];
    System.out.println("Computing thresholds...");
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      if ( (i % 1000) == 0 ) {
        System.out.println("[" + i + "]");
      }
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      GeneData gene = data1.getGeneData(i);
      String id = (String) gene.getDataAt(0);
      Double[] v= gene.getVector(start, end);
      int i1 = i - data1.getNumGeneHeader();
      if (geneid != null && id.equals(geneid)) {
        System.out.println("Gene id " + geneid + " Found at " + i );
        geneid_index = i;
      }
      b[i1] = new Bimodal(v);
      double lowthr1 = b[i1].getLowThreshold();
      double highthr1 = b[i1].getHighThreshold();
      BitSet v1 = new BitSet(data1.getNumArrays());
      BitSet t = new BitSet(data1.getNumArrays());
      v= gene.getVector(start, end);
      for (int j = 0; j < v.length; j++) {
        t.set(j);
        v1.clear(j);
        if (v[j] == null) {
            t.clear(j);
        }
        if (v[j] != null && v[j].doubleValue() > lowthr1 && 
            v[j].doubleValue() < highthr1) {
            gene.setDataAt(start + j, null);
            t.clear(j);
        }
        if (v[j] != null && v[j].doubleValue() >= highthr1) {
            v1.set(j);
        }
      }
      val[i1] = v1;
      thr[i1] = t;
    }
    double thr1 = THRESHOLD;
    int samplesize_cutoff = SIZE_CUTOFF;
    int correction = (int) (samplesize_cutoff * thr1 * 2);

    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    for (int i = 0; i < val.length; i++) {
      int i1 = i + data1.getNumGeneHeader();
      if (geneid_index != -1 && geneid_index != i1) {
        continue;
      }
      System.out.println(i);
      boolean found = false;
      for (int j = i+1; j < val.length; j++) {
        if ( i == j ) continue;
        double[] p = getErrorProbability(val[i], thr[i], val[j], thr[j], correction);
        int j1 = j + data1.getNumGeneHeader();
        if (p[0] <= thr1) {
          String sr = formatString("0.#####", p[0]);
          out.write(0 + "\t" + i1 + "\t" + j1 + "\t" +  sr + "\n");
        }
        if (p[1] <= thr1) {
          String sr = formatString("0.#####", p[1]);
          out.write(1 + "\t" + i1 + "\t" + j1 + "\t" +  sr + "\n");
        }
        if (p[2] <= thr1) {
          String sr = formatString("0.#####", p[2]);
          out.write(2 + "\t" + i1 + "\t" + j1 + "\t" +  sr + "\n");
        }
        if (p[3] <= thr1) {
          String sr = formatString("0.#####", p[3]);
          out.write(3 + "\t" + i1 + "\t" + j1 + "\t" +  sr + "\n");
        }
      }
      // break;
    }
    out.close();
  }

  public static BitSet getBitSet(String str, int type) {
    BitSet res = new BitSet(str.length());
    for (int i =0; i < str.length(); i++) {
      char c = str.charAt(i);
      res.clear(i);
      if (type == 0 && c == '2') {
        res.set(i);
      }
      if (type == 1 && !(c == '1' || c == ' ')) {
        res.set(i);
      }
    }
    return res;
  }

  /*
   * Find pairs with the given Bit vector file (Simple version)
   */
  public static void writePairsBv1(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String bvfile = list.removeFirst();
    String geneid = null;
    int geneid_index = -1;
    if (list.size() > 0) {
       geneid = list.removeFirst();
    }

    PCLFileReader data = new PCLFileReader(bvfile);
    PCLFileReader.CACHE_SIZE = 2000;
    data.beginRandomAccess();
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    double thr1 = THRESHOLD;
    int samplesize_cutoff = SIZE_CUTOFF;
    int corr = (int) (samplesize_cutoff * thr1 * 2);
    int blocksize = 8000;

    GeneData ga = data.getDataAt(0);
    for (int i = 0;  ga != null; i++) {
      String id = (String) ga.getDataAt(0);
      if (geneid != null && id.equals(geneid)) {
        System.out.println("Gene id " + geneid + " Found at " + i );
        geneid_index = i;
      }
      if (geneid != null && geneid_index != i) {
        ga = data.getDataAt(i+1);
        continue;
      }
      System.out.println(i);
      BitSet va = getBitSet((String) ga.getDataAt(2), 0);
      BitSet va_thr = getBitSet((String) ga.getDataAt(2), 1);
      int start = i+1;
      if (geneid != null && geneid_index == i) {
        start = 0;
      }
      GeneData gb = data.getDataAt(start);
      for (int j = start;  gb != null; j++) {
        BitSet vb = getBitSet((String) gb.getDataAt(2), 0);
        BitSet vb_thr = getBitSet((String) gb.getDataAt(2), 1);
        double[] p = getErrorProbability(va, va_thr, vb, vb_thr, corr);
        if (p[0] <= thr1) {
          String sr = formatString("0.#####", p[0]);
          out.write(0 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
        }
        if (p[1] <= thr1) {
          String sr = formatString("0.#####", p[1]);
          out.write(1 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
        }
        if (p[2] <= thr1) {
          String sr = formatString("0.#####", p[2]);
          out.write(2 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
        }
        if (p[3] <= thr1) {
          String sr = formatString("0.#####", p[3]);
          out.write(3 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
        }
        gb = data.getDataAt(j+1);
      } // end gb
      ga = data.getDataAt(i+1);
    } // end ga
    out.close();
  }

  /*
   * Find pairs with the given Bit vector file
   */
  public static void writePairsBv(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String bvfile = list.removeFirst();
    if (list.size() > 0) {
       SIZE_CUTOFF = Integer.parseInt(list.removeFirst());
    }
    if (list.size() > 0) {
       THRESHOLD = Double.parseDouble(list.removeFirst());
    }
    if (list.size() > 0) {
       STAT_THRESHOLD = Double.parseDouble(list.removeFirst());
    }
    String geneid = null;
    int geneid_index = -1;
    if (list.size() > 0) {
       geneid = list.removeFirst();
    }

    PCLFileReader data = new PCLFileReader(bvfile);
    PCLFileReader.CACHE_SIZE = 2000;
    data.beginRandomAccess();
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    double thr1 = THRESHOLD;
    int samplesize_cutoff = SIZE_CUTOFF;
    int corr = (int) (samplesize_cutoff * thr1 * 2);
    int blocksize = 8000;

    GeneData gb1 = data.getDataAt(0);
    for (int b1 = 0;  gb1 != null; b1+=blocksize) {
      BitSet[] ba1 = new BitSet[blocksize];
      BitSet[] ba1_thr = new BitSet[blocksize];
      geneid_index = -1;
      GeneData ga = data.getDataAt(b1);
      for (int i = b1; ga != null && i < (b1+blocksize); i++) {
        ga = data.getDataAt(i);
        if (ga == null) {
            break;
        }
        String id = (String) ga.getDataAt(0);
        if (geneid != null && id.equals(geneid)) {
          System.out.println("Gene id " + geneid + " Found at " + i );
          geneid_index = i;
        }
        BitSet va = getBitSet((String) ga.getDataAt(2), 0);
        BitSet va_thr = getBitSet((String) ga.getDataAt(2), 1);
        ba1[i-b1] = va;
        ba1_thr[i-b1] = va_thr;
      }
      if (geneid != null && geneid_index == -1) {
        gb1 = data.getDataAt(b1+blocksize);
        continue;
      }
      int start = b1;
      if (geneid != null && geneid_index != -1) {
        start = 0;
      }
      GeneData gb2 = data.getDataAt(start);
      for (int b2 = start;  gb2 != null; b2+=blocksize) {
        System.out.println("Block = (" + b1 + ", " + b2 + ")");
        BitSet[] ba2 = new BitSet[blocksize];
        BitSet[] ba2_thr = new BitSet[blocksize];
        GeneData gb = data.getDataAt(b2);
        for (int i = b2; gb != null && i < (b2+blocksize); i++) {
          gb = data.getDataAt(i);
          if (gb == null) {
            break;
          }
          BitSet vb = getBitSet((String) gb.getDataAt(2), 0);
          BitSet vb_thr = getBitSet((String) gb.getDataAt(2), 1);
          ba2[i-b2] = vb;
          ba2_thr[i-b2] = vb_thr;
        }
        BitSet va = ba1[0];
        for (int i = b1; va != null && i < (b1+blocksize); i++) {
          if (geneid_index != -1 && geneid_index != i) {
            if (i < (b1-1 +blocksize)) {
              va = ba1[i+1-b1];
            }
            continue;
          }
          System.out.println(i);
          va = ba1[i-b1];
          BitSet va_thr = ba1_thr[i-b1];
          int num = va.size();
          int outside = va_thr.cardinality();
          if (num > (3 * outside)) {
            if (i < (b1-1 +blocksize)) {
              va = ba1[i+1-b1];
            }
            continue;
          }
          BitSet vb = ba2[0];
          for (int j = b2; vb != null && j < (b2 + blocksize); j++) {
            vb = ba2[j-b2];
            BitSet vb_thr = ba2_thr[j-b2];
            num = vb.size();
            outside = vb_thr.cardinality();
            if ((num > (3 *outside)) || (geneid_index == -1 && j <= i)) {
              if (j < (b2 -1 + blocksize)) {
                vb = ba2[j+1-b2];
              }
              continue;
            }
            //double[] p = getErrorProbability(va, va_thr, vb, vb_thr, corr);
            //double[] p = getErrorStats(va, va_thr, vb, vb_thr, corr, 0);
            double[] p = getErrorProbStats(va, va_thr, vb, vb_thr, 0);
            if (p[0] <= thr1) {
              String sr = formatString("0.#####", p[0]);
              out.write(0 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
            }
            if (p[1] <= thr1) {
              String sr = formatString("0.#####", p[1]);
              out.write(1 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
            }
            if (p[2] <= thr1) {
              String sr = formatString("0.#####", p[2]);
              out.write(2 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
            }
            if (p[3] <= thr1) {
              String sr = formatString("0.#####", p[3]);
              out.write(3 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
            }
            if (j < (b2 -1 + blocksize)) {
              vb = ba2[j+1-b2];
            }
          }
          if (i < (b1-1 +blocksize)) {
            va = ba1[i+1-b1];
          }
        } // end ga
        gb2 = data.getDataAt(b2+blocksize);
      } // end gb2
      gb1 = data.getDataAt(b1+blocksize);
    } // end gb1
    out.close();
  }

  /**
   * Find pairs with the given Bit vector file 
   *  For use in Parallel Array of processors
   */
  public static void writePairsBvParallel(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String bvfile = list.removeFirst();
    int start1 = Integer.parseInt(list.removeFirst());
    int num1 = Integer.parseInt(list.removeFirst());
    int start2 = Integer.parseInt(list.removeFirst());
    int num2 = Integer.parseInt(list.removeFirst());
    if (list.size() > 0) {
       SIZE_CUTOFF = Integer.parseInt(list.removeFirst());
    }
    if (list.size() > 0) {
       THRESHOLD = Double.parseDouble(list.removeFirst());
    }
    if (list.size() > 0) {
       STAT_THRESHOLD = Double.parseDouble(list.removeFirst());
    }

    PCLFileReader data = new PCLFileReader(bvfile);
    PCLFileReader.CACHE_SIZE = 2000;
    data.beginRandomAccess();
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    double thr1 = THRESHOLD;
    int samplesize_cutoff = SIZE_CUTOFF;
    int corr = (int) (samplesize_cutoff * thr1 * 2);

    BitSet[] ba1 = new BitSet[num1];
    BitSet[] ba1_thr = new BitSet[num1];
    GeneData ga = data.getDataAt(start1);
    for (int i = start1; ga != null && i < (start1 + num1); i++) {
      ga = data.getDataAt(i);
      if (ga == null) {
        break;
      }
      //System.out.println(ga.getDataAt(0));
      BitSet va = getBitSet((String) ga.getDataAt(2), 0);
      BitSet va_thr = getBitSet((String) ga.getDataAt(2), 1);
      ba1[i-start1] = va;
      ba1_thr[i-start1] = va_thr;
    }
    BitSet[] ba2 = new BitSet[num2];
    BitSet[] ba2_thr = new BitSet[num2];
    GeneData gb = data.getDataAt(start2);
    for (int i = start2; gb != null && i < (start2 + num2); i++) {
      gb = data.getDataAt(i);
      if (gb == null) {
        break;
      }
      //System.out.println(gb.getDataAt(0));
      BitSet vb = getBitSet((String) gb.getDataAt(2), 0);
      BitSet vb_thr = getBitSet((String) gb.getDataAt(2), 1);
      ba2[i-start2] = vb;
      ba2_thr[i-start2] = vb_thr;
    }
    BitSet va = ba1[0];
    for (int i = start1; va != null && i < (start1 + num1); i++) {
      System.out.println(i);
      va = ba1[i-start1];
      BitSet va_thr = ba1_thr[i-start1];
      int num = va.size();
      int outside = va_thr.cardinality();
      //System.out.println("Num = " + num + " Outside = " + outside);
      if (num > (3 * outside)) {
        if (i < (start1 -1 + num1)) {
          va = ba1[i+1-start1];
        }
        continue;
      }
      BitSet vb = ba2[0];
      for (int j = start2; vb != null && j < (start2 + num2); j++) {
        vb = ba2[j-start2];
        BitSet vb_thr = ba2_thr[j-start2];
        num = vb.size();
        outside = vb_thr.cardinality();
        //System.out.println("Num = " + num + " Outside = " + outside);
        if (num > (3 *outside)) {
          if (j < (start2 -1 + num2)) {
            vb = ba2[j+1-start2];
          }
          continue;
        }
        //double[] p = getErrorProbability(va, va_thr, vb, vb_thr, corr);
        //double[] p = getErrorStats(va, va_thr, vb, vb_thr, corr, 1);
        double[] p = getErrorProbStats(va, va_thr, vb, vb_thr, 1);
        if (p[0] <= thr1) {
          String sr = formatString("0.#####", p[0]);
          out.write(0 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
        }
        if (p[1] <= thr1) {
          String sr = formatString("0.#####", p[1]);
          out.write(1 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
        }
        if (p[2] <= thr1) {
          String sr = formatString("0.#####", p[2]);
          out.write(2 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
        }
        if (p[3] <= thr1) {
          String sr = formatString("0.#####", p[3]);
          out.write(3 + "\t" + i + "\t" + j + "\t" +  sr + "\n");
        }
        if (j < (start2 -1 + num2)) {
          vb = ba2[j+1-start2];
        }
      }
      if (i < (start1 -1 + num1)) {
        va = ba1[i+1-start1];
      }
    } // end ga
    out.close();
  }

  /*
   * Find conserved pairs with the given PCL file
   */
  public static void writeCommonPairs(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String pclfile = list.removeFirst();
    String pairfile = list.removeFirst();
    String cfile = list.removeFirst();
    String homologFile = null;
    if (list.size() > 0) {
       homologFile = list.removeFirst();
    }

    HashMap<String, String> homolog = new HashMap<String, String>();
    if (homologFile != null) {
      PCLFileReader hdata = new PCLFileReader(homologFile);
      hdata.begin();
      for (int i = 1; hdata.hasNext(); i++) {
        GeneData gene = hdata.getData();
        if (gene == null) {
          break;
        }
        String human = (String) gene.getDataAt(0);
        String mouse = (String) gene.getDataAt(1);
        String fly = (String) gene.getDataAt(2);
        human = human.toUpperCase();
        mouse = mouse.toUpperCase();
        fly = fly.toUpperCase();
        homolog.put(human+"_Mm", mouse);
        homolog.put(human+"_Dm", fly);
      }
    }

    PCLFileReader data = new PCLFileReader(cfile);
    PCLFileReader.CACHE_SIZE = 20000;
    HashMap<String, LinkedList<Long> > b_ids = new HashMap<String, LinkedList<Long> >();
    data.beginRandomAccess();
    GeneNameScheme ns = new GeneNameScheme(1, "Hs", ":", 0);
    for (int i = 0; data.hasNext(); i++) {
      GeneData gene = data.getDataAt(-1);
      if (gene == null) {
        break;
      }
      String str = (String) gene.getDataAt(1);
      String name = ns.getGene(str);
      if (name == null) {
        continue;
      }
      name = name.replaceAll("\\s", "");
      name = name.toUpperCase();
      if (!b_ids.containsKey(name)) {
        b_ids.put(name, new LinkedList<Long>());
      }
      LinkedList<Long> val = b_ids.get(name);
      val.add(new Long(i));
    }

    PCLFileReader pcldata = new PCLFileReader(pclfile);
    HashMap<Integer, String> indexhash = new HashMap<Integer, String>();
    pcldata.begin();
    for (int i = 1; pcldata.hasNext(); i++) {
      GeneData gene = pcldata.getData();
      if (gene == null) {
        break;
      }
      String str = (String) gene.getDataAt(1);
      String name = ns.getGene(str);
      if (name == null || name.equals("---") || name.equals("")) {
        continue;
      }
      name = name.replaceAll("\\s", "");
      name = name.toUpperCase();
      indexhash.put(new Integer(i), name);
    }

    TABFileReader pairdata = new TABFileReader(pairfile);
    pairdata.begin();
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    GeneData pair = pairdata.getHeader();
    HashMap<Integer, BitSet> cache = new HashMap<Integer, BitSet>();
    HashMap<Integer, BitSet> cache_thr = new HashMap<Integer, BitSet>();

    double thr1 = 0.10;
    int samplesize_cutoff = 20;
    int corr = (int) (samplesize_cutoff * thr1 * 2);

    for (int i = 0; pair != null; i++, pair = pairdata.getData()) {
      if ( (i % 10000) == 0 ) {
        System.out.println("[" + i + "]");
      }
      int q = Integer.parseInt((String)pair.getDataAt(0));
      int x = Integer.parseInt((String)pair.getDataAt(1));
      int y = Integer.parseInt((String)pair.getDataAt(2));
      double w = Double.parseDouble((String)pair.getDataAt(3));
      if ( !indexhash.containsKey(new Integer(x)) ||
          !indexhash.containsKey(new Integer(y)) ) {
        continue;
      }

      // Gene names
      String namex = indexhash.get(new Integer(x));
      String namey = indexhash.get(new Integer(y));

      if (namex.equals(namey)) {
        continue;
      }
      // Searching in homolog
      if (homolog.containsKey(namex + "_Mm")) {
        String str = homolog.get(namex + "_Mm");
        if (b_ids.containsKey(str)) {
            namex = str;
        }
      }
      if (homolog.containsKey(namex + "_Dm")) {
        String str = homolog.get(namex + "_Dm");
        if (b_ids.containsKey(str)) {
            namex = str;
        }
      }
      if (homolog.containsKey(namey + "_Mm")) {
        String str = homolog.get(namey + "_Mm");
        if (b_ids.containsKey(str)) {
            namey = str;
        }
      }
      if (homolog.containsKey(namey + "_Dm")) {
        String str = homolog.get(namey + "_Dm");
        if (b_ids.containsKey(str)) {
            namey = str;
        }
      }
      if ( !b_ids.containsKey(namex) || !b_ids.containsKey(namey) ) {
        continue;
      }
      LinkedList<Long> listx = b_ids.get(namex);
      LinkedList<Long> listy = b_ids.get(namey);
      ListIterator<Long> itrx = listx.listIterator();

      //System.out.println(listx.size() + "\t" + namex);
      //System.out.println(listy.size() + "\t" + namey);

      boolean found = false;
      int a = -1 , b = -1;
      double pvalue = 1.0;

      while (itrx.hasNext()) {
        Long lx = (Long) itrx.next();
        a = lx.intValue();
        ListIterator<Long> itry = listy.listIterator();
        while (itry.hasNext()) {
          Long ly = (Long) itry.next();
          b = ly.intValue();
          BitSet va = null;
          BitSet vb = null;
          BitSet va_thr = null;
          BitSet vb_thr = null;
          if (cache.containsKey(new Integer(a))) {
            va = cache.get(new Integer(a));
            va_thr = cache_thr.get(new Integer(a));
          }
          if (cache.containsKey(new Integer(b))) {
            vb = cache.get(new Integer(b));
            vb_thr = cache_thr.get(new Integer(b));
          }
          if (va == null) {
            GeneData ga = data.getDataAt(a);
            va = getBitSet((String) ga.getDataAt(2), 0);
            va_thr = getBitSet((String) ga.getDataAt(2), 1);
            cache.put(new Integer(a), va);
            cache_thr.put(new Integer(a), va_thr);
          }
          if (vb == null) {
            GeneData gb = data.getDataAt(b);
            vb = getBitSet((String) gb.getDataAt(2), 0);
            vb_thr = getBitSet((String) gb.getDataAt(2), 1);
            cache.put(new Integer(b), vb);
            cache_thr.put(new Integer(b), vb_thr);
          }
          if (cache.size() > 10000) {
            System.out.println("Vector cache clear");
            cache.clear();
            cache_thr.clear();
          }
          double[] p = getErrorProbability(va, va_thr, vb, vb_thr, corr);
          /*
          System.out.println(p[0] + "\t" + p[1] + "\t" + p[2] + "\t" + p[3]);
          System.out.println(q + "\t" + w );
          System.out.println(a + "\t" + b );
          System.out.println(x + "\t" + y );
          */
          if (p[q] <= thr1) {
            found = true;
            pvalue = p[q];
            break;
          }
        } // end itry
        if (found) {
          break;
        }
      } // end itrx
      if (found) {
        String sr = formatString("0.#####", pvalue);
        out.write(q + "\t" + a + "\t" + b + "\t" +  sr + "\t" + pair);
      }
    } // end for
    out.close();
  }

  public static void plotPairs(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String pclfile = list.removeFirst();
    String pairfile = list.removeFirst();
    Data data1 = PCLFileReader.readFile(pclfile);
    Data data2 = PCLFileReader.readFile(pairfile);

    Bimodal[] bi = new Bimodal[data1.getNumGenes()];
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      // System.out.println("[ " + i + " ]---------------------->");
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v= data1.getGeneData(i).getVector(start, end);
      bi[i - data1.getNumGeneHeader()] = new Bimodal(v);
    }

    Color[] colors = new Color[data1.getNumArrays()];
    for (int i =0; i < colors.length; i++) {
        colors[i] = Color.red;
    }

    GeneNameScheme ns = new GeneNameScheme(1, "Hs", ":", 0);
    ns.setNumMissingPoints(0);
    data1.setGeneNameScheme(ns);
    PSPlot plotter = new PSPlot(ofile);
    plotter.open();
    plotter.array(1 /* Rows */, 1 /* columns */);

    if (list.size() > 0) {
      Double[] cd = new Double[colors.length];
      Double[] t = new Double[colors.length];
      for (int i =0; i < colors.length; i++) {
        cd[i] = new Double(i);
        t[i] = new Double(i);
      }
      plotter.setupPlot(t, cd);
      String colorfile = list.removeFirst();
      Data cdata = TABFileReader.readFile(colorfile);
      HashMap<String, Integer> map = new HashMap<String, Integer>();
      int cindex = 0;
      for (int i = 0; i < cdata.getNumRows(); i++) {
        if (i < colors.length) {
          GeneData gene = cdata.getGeneData(i);
          String sym = (String)gene.getDataAt(1);
          if (!map.containsKey(sym)) {
            cindex += 0x19;
            Integer ci = new Integer(cindex);
            if (sym.equals("4")) {
                ci =  new Integer(0);
            }
            map.put(sym, ci);
            int r = ((ci.intValue() & 0x30) >> 4) * 85;
            int g = ((ci.intValue() & 0xc) >> 2) * 85;
            int b = (ci.intValue() & 0x3) * 85;
            plotter.setRGBcolor(new Color(r,g,b));
            plotter.text(t[t.length/2], new Double(ci.intValue() * 6), 0, sym);
          }
          Integer ci = map.get(sym);
          int r = ((ci.intValue() & 0x30) >> 4) * 85;
          int g = ((ci.intValue() & 0xc) >> 2) * 85;
          int b = (ci.intValue() & 0x3) * 85;
          // System.out.println(ci + ":" + r + "," + g + "," + b);
          colors[i] = new Color(r,g,b);
        }
        else {
          break;
        }
      }
      plotter.newPage();
    }

    int index = 0;
    plotter.array(4 /* Rows */, 2 /* columns */);
    for (int i = 0; i < data2.getNumRows(); i++) {
      // System.out.println(i);
      GeneData gene = data2.getGeneData(i);
      String sym = (String)gene.getDataAt(2);
      sym = sym.replaceAll("\\s", "");
      if ((sym.equals("+") || sym.equals("-"))) {
        index++;
        if (index > 5000) {
            break;
        }
        System.out.println(sym);
        int start = data1.getNumArrayHeader();
        int end = data1.getNumColumns() - 1;
        int a = Integer.parseInt((String)gene.getDataAt(0));
        int b = Integer.parseInt((String)gene.getDataAt(1));
        Double[] data= data1.getGeneData(a).getVector(start, end);
        Double[] time= data1.getGeneData(b).getVector(start, end);
        String name1 = data1.getGenesAt(a);
        String name2 = data1.getGenesAt(b);
        plotter.setupPlot(time, data);
        for (int j=0; j < data.length; j++) {
          plotter.setRGBcolor(colors[j]);
          plotter.point(time[j], data[j]);
        }
        plotter.setRGBcolor(Color.black);
        plotter.verticalLine(new Double(bi[b - data1.getNumGeneHeader()].getThreshold()));
        plotter.horizontalLine(new Double(bi[a - data1.getNumGeneHeader()].getThreshold()));
        plotter.xlabel(name2);
        plotter.ylabel(name1);
        plotter.title("Pairs");
      }
    }
    plotter.close();

  }

  public static double getImpStat(int index, double thr1, double thr2, Double[] v1, Double[] v2) {
    double res = 0.0;
    int count = 0;
    int count1 = 0;
    int count2 = 0;
    for (int i = 0; i < v1.length; i++) {
      if (v1[i] != null && v2[i] != null) {
        if ( (index == 0 || index == 1) && v1[i].doubleValue() <= thr1) count1++;
        if ( (index == 0 || index == 3) && v2[i].doubleValue() <= thr2) count2++;
        if ( (index == 2 || index == 3) && v1[i].doubleValue() > thr1) count1++;
        if ( (index == 1 || index == 2) && v2[i].doubleValue() > thr2) count2++;
        if (index == 0 && v1[i].doubleValue() <= thr1 && v2[i].doubleValue() <= thr2) count++;
        if (index == 1 && v1[i].doubleValue() <= thr1 && v2[i].doubleValue() > thr2) count++;
        if (index == 2 && v1[i].doubleValue() > thr1 && v2[i].doubleValue() > thr2) count++;
        if (index == 3 && v1[i].doubleValue() > thr1 && v2[i].doubleValue() <= thr2) count++;
      }
    }
    res = 0.5 * count / count1 + 0.5 * count / count2;
    return res;
  }

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

  public static void bimodalAnalysis(LinkedList<String> list) throws Exception {
    String cmd = list.removeFirst();

    if (cmd.equals("get")) {
       writeBimodalGenes(false, list);
    }
    if (cmd.equals("getr")) {
       writeBimodalGenes(true, list);
    }
    if (cmd.equals("thr")) {
       writeThreshold(list);
    }
    if (cmd.equals("thr1")) {
       writeThreshold1(list);
    }
    if (cmd.equals("thrBv")) {
       writeThresholdBv(list);
    }
    if (cmd.equals("thrBv1")) {
       writeThresholdBv1(list);
    }
    if (cmd.equals("pairs")) {
       writePairsNew(list);
    }
    if (cmd.equals("pairsBv")) {
       writePairsBv(list);
    }
    if (cmd.equals("pairsBv1")) {
       writePairsBv1(list);
    }
    if (cmd.equals("pairsBvParallel")) {
       writePairsBvParallel(list);
    }
    if (cmd.equals("commonPairs")) {
       writeCommonPairs(list);
    }
    if (cmd.equals("targets")) {
       writePairsNew(list);
    }
    if (cmd.equals("groups")) {
       GroupAnalysis.groupAnalysis(list);
    }
    if (cmd.equals("plot")) {
       plotPairs(list);
    }
  }

  public static void geoAnalysis(LinkedList<String> list) throws Exception {
    String ofile1 = list.removeFirst();
    String ofile2 = list.removeFirst();
    String ifile = list.removeFirst();
    GEOFileReader.MAX_ARRAYS = Integer.parseInt(list.removeFirst());
    GEOFileReader.RANDOM = Boolean.valueOf(list.removeFirst()).booleanValue();
    if (list.size() > 0) {
      GEOFileReader.SYM_INDEX = Integer.parseInt(list.removeFirst());
    }
    if (list.size() > 0) {
      GEOFileReader.TITLE_INDEX = Integer.parseInt(list.removeFirst());
    }
    if (GEOFileReader.RANDOM) {
      System.out.println("Random selection activated");
    }
    HashSet<String> excludeList = new HashSet<String>();
    if (list.size() > 0) {
      String infoFile = list.removeFirst();
      Data cdata = TABFileReader.readFile(infoFile);
      for (int i = 0; i < cdata.getNumRows(); i++) {
        GeneData gene = cdata.getGeneData(i);
        String str = (String)gene.getDataAt(0);
        str = str.replaceAll("\\s", "");
        excludeList.add(str);
      }
    }
    Data data = GEOFileReader.readFile(ifile, excludeList);
    PCLFileWriter.writeFile(data, ofile1);
    data = GEOFileReader.readInfo(ifile, excludeList);
    PCLFileWriter.writeFile(data, ofile2);
  }

  public static void geoAnalysis1(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String ifile = list.removeFirst();
    Data data = GEOFileReader.readSampleIds(ifile);
    PCLFileWriter.writeFile(data, ofile);
  }

  public static void GSAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String org = list.removeFirst();
    double pvalue = Double.parseDouble(list.removeFirst());
    String file1 = list.removeFirst();
    GeneSet a = GeneSet.readFile(file1);
    if (list.size() > 0) {
      String file2 = list.removeFirst();
      GeneSet b = GeneSet.readFile(file2);
      GeneSetAnalysis ana = new GeneSetAnalysis(a, b);
      if (ofile.endsWith(".html")) {
        ana.performAnalysis(ofile, org, pvalue);
      }
      else {
        ana.performAnalysisOutTAB(ofile, org, pvalue);
      }
    }
    else {
      GeneSet.writeFile(a, ofile, org);
    }
  }

  public static void GOAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String onnFile = list.removeFirst();
    String annFile = list.removeFirst();
    String org = list.removeFirst();
    double pvalue = Double.parseDouble(list.removeFirst());
    String file1 = list.removeFirst();
    GeneSet a = GeneSet.readFile(file1);
    GeneSetAnalysis ana = new GeneSetAnalysis(a);
    Vector<String> allNodesV = new Vector<String>(a.getAllGenes());
    System.out.println("Number of genes = " + allNodesV.size());
    if (ofile.endsWith(".html")) {
      ana.performGOAnalysis(ofile, onnFile, annFile, org, pvalue, allNodesV);
    }
    else {
      ana.performGOAnalysisOutTAB(ofile, onnFile, annFile, org, pvalue, allNodesV);
    }
  }

  public static void selectArrayPCLAnalysis(boolean remove, LinkedList<String> l) throws Exception {
    String ofile = l.removeFirst();
    String ifile = l.removeFirst();
    HashSet<String> list = new HashSet<String>();
    if (l.size() > 0) {
      String infoFile = l.removeFirst();
      Data cdata = TABFileReader.readFile(infoFile);
      for (int i = 0; i < cdata.getNumRows(); i++) {
        GeneData gene = cdata.getGeneData(i);
        String str = (String)gene.getDataAt(0);
        str = str.replaceAll("\\s", "");
        //System.out.println(":"+str+":");
        list.add(str);
      }
    }
    /*
    Data data = PCLFileReader.readFile(ifile);
    GeneData gene = data.getGeneData(0);
    */
    PCLFileReader data = new PCLFileReader(ifile);
    data.begin();
    GeneData gene = data.getHeader();
    Object[] header = gene.getData();
    int count = 0;
    for (int i = 0; i < data.getNumArrays(); i++) {
      String hdr = (String) header[i+data.getNumArrayHeader()];
      hdr = hdr.replaceAll("\\s", "");
      //System.out.println(":"+hdr+":");
      if (remove && !list.contains(hdr)) {
        count++;
      }
      if (!remove && list.contains(hdr)) {
        count++;
      }
    }
    System.out.println("Num Array : " + data.getNumArrays() + " -> " + count);
    if (count > 0) {
      int[] selectOrder = new int[count+data.getNumArrayHeader()];
      int index = 0;
      for (int i =0; i < data.getNumArrayHeader(); i++) {
        selectOrder[index] = i;
        index++;
      }
      for (int i = data.getNumArrayHeader(); i < data.getNumColumns(); i++) {
        String hdr = (String) header[i];
        hdr = hdr.replaceAll("\\s", "");
        if (remove && !list.contains(hdr)) {
          selectOrder[index] = i;
          index++;
        }
        if (!remove && list.contains(hdr)) {
          selectOrder[index] = i;
          index++;
        }
      }
      /*
      // Change selectOrder - permute the data only
      data = Data.selectArraysFromData(data, selectOrder);
      //data.print();
      PCLFileWriter.writeFile(data, ofile);
       */
      PCLFileReader reader = data;
      PCLFileWriter writer = new PCLFileWriter(ofile);
      GeneData h = reader.getHeader();
      writer.writeData(h.permute(0, selectOrder));
      writer.writeData(reader.getData().permute(0, selectOrder));
      while(reader.hasNext()) {
        gene = reader.getData();
        if (gene == null) {
          break;
        }
        if ( (reader.getLineNumber() % 1000 == 0) ) {
          System.out.println(reader.getLineNumber());
        }
        writer.writeData(gene.permute(0, selectOrder));
      }
      writer.close();
    }
  }

  public static void selectArrayAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String ifile = list.removeFirst();
    String geneName = list.removeFirst();
    double minthreshold = Double.parseDouble(list.removeFirst());
    double maxthreshold = Double.parseDouble(list.removeFirst());
    Data data = PCLFileReader.readFile(ifile);
    GeneData g = data.getGeneData(0);
    Object[] header = g.getData();
    geneName = geneName.replaceAll("\\s", "");
    for (int i = 0; i < data.getNumRows(); i++) {
      GeneData gene = data.getGeneData(i);
      String str = (String)gene.getDataAt(0);
      str = str.replaceAll("\\s", "");
      if (str.equals(geneName)) {
        System.out.println(geneName + " Found at #" + i);
        Double[] d = gene.getVector(data.getNumArrayHeader(), data.getNumColumns() - 1);
        Vector<GeneData> resData = new Vector<GeneData>();
        for (int j = 0; j < d.length; j++) {
          String h = (String) header[j+data.getNumArrayHeader()];
          Double d1 = d[j];
          if (d1 != null && 
              (d1.doubleValue() > minthreshold &&
               d1.doubleValue() < maxthreshold )) {
            Object[] da = new Object[2];
            da[0] = h;
            da[1] = d1;
            GeneData dag = new GeneData(da);
            resData.add(dag);
          }
        }
        Data res = new Data(2, resData.size(), 0, 0, resData);
        PCLFileWriter.writeFile(res, ofile);
        break;
      }
    }
  }

  public static void reduceLog(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String ifile = list.removeFirst();
    /*
    Data data = PCLFileReader.readFile(ifile);
    data.convertDoubles();
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      if ( (i%500) == 0) {
        System.out.println("[ " + i + " ]");
      }
      GeneData gene = data.getGeneData(i);
      Object[] geneData = gene.getData();
      for (int j = start; j <= end; j++) {
        Double v = (Double) geneData[j];
        if (v != null && v.doubleValue() <= 0) {
          geneData[j] = null;
          System.out.println("*Warning: " + v + " at ("+i+","+j+")");
        }
      }
    }
    data.reduceLog();
    PCLFileWriter.writeFile(data, ofile);
    */
    PCLFileReader reader = new PCLFileReader(ifile);
    PCLFileWriter writer = new PCLFileWriter(ofile);
    reader.begin();
    GeneData h = reader.getHeader();
    writer.writeData(h);
    writer.writeData(reader.getData());
    int start = reader.getNumArrayHeader();
    int end = reader.getNumColumns() - 1;
    while(reader.hasNext()) {
      GeneData gene = reader.getData();
      if (gene == null) {
        break;
      }
      if ( (reader.getLineNumber() % 1000 == 0) ) {
        System.out.println(reader.getLineNumber());
      }
      gene.reduceLog(start, end);
      writer.writeData(gene);
    }
    writer.close();
  }

  public static void normalizeGene(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String ifile = list.removeFirst();
    String geneid = list.removeFirst();
    geneid = geneid.replaceAll("\\s", "");
    Data data = PCLFileReader.readFile(ifile);
    int index = -1;
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      if ( (i%500) == 0) {
        System.out.println("[ " + i + " ]");
      }
      GeneData gene = data.getGeneData(i);
      String id = (String) gene.getDataAt(0);
      id = id.replaceAll("\\s", "");
      if (id.equals(geneid)) {
        System.out.println("Gene id : " + geneid + " found at " + i);
        index = i;
        break;
      }
    }
    if (index < 0) {
        System.out.println("Gene id : " + geneid + " not found");
        return;
    }
    data.convertDoubles();
    GeneData refGene = (GeneData) data.getGeneData(index).clone();
    Object[] refData = refGene.getData();
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      if ( (i%500) == 0) {
        System.out.println("[ " + i + " ]");
      }
      int start = data.getNumArrayHeader();
      int end = data.getNumColumns() - 1;
      GeneData gene = data.getGeneData(i);
      Object[] geneData = gene.getData();
      for (int j = start; j <= end; j++) {
        Double v = (Double) geneData[j];
        Double r = (Double) refData[j];
        if (v != null && r != null) {
          geneData[j] = new Double(v.doubleValue() - r.doubleValue());
        }
        else {
          geneData[j] = null;
        }
      }
    }
    PCLFileWriter.writeFile(data, ofile);
  }

  public static void selectNames(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String ifile = list.removeFirst();
    String nfile = list.removeFirst();
    Data data = PCLFileReader.readFile(ifile);
    Data datan = TABFileReader.readFile(nfile);

    GeneNameScheme ns;
    ns = new GeneNameScheme(1, "Dm", ":", 0);
    ns = new GeneNameScheme(1, "Mm", "\\|\\|", 3);
    ns = new GeneNameScheme(1, "Hs", ":", 0);
    data.setGeneNameScheme(ns);

    Data res = Data.containNames(data, datan);
    PCLFileWriter.writeFile(res, ofile);
  }

  public static void addConstant(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String ifile = list.removeFirst();
    double value = Double.parseDouble(list.removeFirst());
    Data data = PCLFileReader.readFile(ifile);
    data.convertDoubles();
    for (int i = data.getNumGeneHeader(); i < data.getNumRows(); i++) {
      if ( (i%500) == 0) {
        System.out.println("[ " + i + " ]");
      }
      int start = data.getNumArrayHeader();
      int end = data.getNumColumns() - 1;
      GeneData gene = data.getGeneData(i);
      Object[] geneData = gene.getData();
      for (int j = start; j <= end; j++) {
        Double v = (Double) geneData[j];
        if (v != null) {
          geneData[j] = new Double(v.doubleValue() + value);
        }
        else {
          geneData[j] = null;
        }
      }
    }
    PCLFileWriter.writeFile(data, ofile);
  }

  public static void filterAnalysis(LinkedList<String> list) throws Exception {
    String cmd = list.removeFirst();
    if (cmd.equals("reduceLog")) {
      reduceLog(list);
    }
    if (cmd.equals("normalize")) {
      normalizeGene(list);
    }
    if (cmd.equals("selectNames")) {
      selectNames(list);
    }
    if (cmd.equals("add")) {
      addConstant(list);
    }
  }

  public static void tabFileAnalysis(LinkedList<String> list) throws Exception {
    String cmd = list.removeFirst();
    String ofile = list.removeFirst();
    if (cmd.equals("concat")) {
      String file1 = list.removeFirst();
      String file2 = list.removeFirst();
      Data data1 = TABFileReader.readFile(file1);
      Data data2 = TABFileReader.readFile(file2);
      Data res = Data.concatDataColumns(data1, data2);
      PCLFileWriter.writeFile(res, ofile);
    }
    if (cmd.equals("select")) {
      String file1 = list.removeFirst();
      String[] range = list.removeFirst().split(":");
      int start = Integer.parseInt(range[0]);
      int end = start;
      if (range.length > 1) {
        end = Integer.parseInt(range[1]);
      }
      if (end < start) { end = start; }
      TABFileReader reader = new TABFileReader(file1);
      PCLFileWriter writer = new PCLFileWriter(ofile);
      reader.begin();
      GeneData h1 = reader.getHeader();
      int[] order = new int[end - start + 1];
      for (int i = start; i <= end; i++) { order[i-start] = i; }
      GeneData h = h1.permute(0, order);
      writer.writeData(h);
      while(reader.hasNext()) {
        GeneData gene = reader.getData();
        if (gene == null) {
            break;
        }
        if ( (reader.getLineNumber() % 1000 == 0) ) {
            System.out.println(reader.getLineNumber());
        }
        GeneData g = gene.permute(0, order);
        writer.writeData(g);
      }
      writer.close();
    }
    if (cmd.equals("delete")) {
      String file1 = list.removeFirst();
      String[] range = list.removeFirst().split(":");
      int start = Integer.parseInt(range[0]);
      int end = start;
      if (range.length > 1) {
        end = Integer.parseInt(range[1]);
      }
      if (end < start) { end = start; }
      TABFileReader reader = new TABFileReader(file1);
      PCLFileWriter writer = new PCLFileWriter(ofile);
      reader.begin();
      GeneData h1 = reader.getHeader();
      int length = end - start + 1;
      int[] order = new int[h1.size() - length];
      for (int i = 0; i < start; i++) { order[i] = i; }
      for (int i = end + 1; i < h1.size(); i++) { order[i-length] = i; }
      GeneData h = h1.permute(0, order);
      writer.writeData(h);
      while(reader.hasNext()) {
        GeneData gene = reader.getData();
        if (gene == null) {
            break;
        }
        if ( (reader.getLineNumber() % 1000 == 0) ) {
            System.out.println(reader.getLineNumber());
        }
        GeneData g = gene.permute(0, order);
        writer.writeData(g);
      }
      writer.close();
    }
    if (cmd.equals("intersect")) {
      String file2 = list.removeFirst();
      String file3 = list.removeFirst();
      Data data = TABFileReader.readFile(file3);
      HashMap<String, LinkedList<Integer> > b_ids = new HashMap<String, LinkedList<Integer> >();
      for (int i = 0; i < data.getNumGenes(); i++) {
        GeneData gene = data.getGeneData(i);
        String id = (String) gene.getDataAt(data.getID());
        id = id.replaceAll("\\s", "");
        if (!b_ids.containsKey(id)) {
          b_ids.put(id, new LinkedList<Integer>());
        }
        LinkedList<Integer> val = b_ids.get(id);
        val.add(new Integer(i));
      }
      TABFileReader reader = new TABFileReader(file2);
      PCLFileWriter writer = new PCLFileWriter(ofile);
      reader.begin();
      GeneData gene = reader.getHeader();
      GeneData g = data.getGeneData(0);
      String id = (String) gene.getDataAt(0);
      id = id.replaceAll("\\s", "");
      if (b_ids.containsKey(id)) {
        g = GeneData.merge(gene, g.subset(1,g.size()-1));
        writer.writeData(g);
      }
      while(reader.hasNext()) {
        gene = reader.getData();
        if (gene == null) {
            break;
        }
        if ( (reader.getLineNumber() % 1000 == 0) ) {
            System.out.println(reader.getLineNumber());
        }
        id = (String) gene.getDataAt(0);
        id = id.replaceAll("\\s", "");
        if (b_ids.containsKey(id)) {
          LinkedList<Integer> val = b_ids.get(id);
          Iterator<Integer> itr = val.iterator();
          while (itr.hasNext()) {
            Integer e = (Integer) itr.next();
            g = data.getGeneData(e.intValue());
            g = GeneData.merge(gene, g.subset(1,g.size()-1));
            writer.writeData(g);
          }
        }
      }
      writer.close();
    }
    if (cmd.equals("intersectR")) {
      String file1 = list.removeFirst();
      String file2 = list.removeFirst();
      PCLFileReader data = new PCLFileReader(file2);
      PCLFileReader.CACHE_SIZE = 0;
      HashMap<String, LinkedList<Long> > b_ids = new HashMap<String, LinkedList<Long> >();
      data.beginRandomAccess();
      for (int i = 0; data.hasNext(); i++) {
        GeneData gene = data.getDataAt(-1);
        if (gene == null) {
            break;
        }
        if ( (i % 1000) == 0) {
          System.gc();
          System.out.println("[" + i + "]");
        }
        String id = (String) gene.getDataAt(0);
        id = id.replaceAll("\\s", "");
        id = new String(id);
        if (!b_ids.containsKey(id)) {
          b_ids.put(id, new LinkedList<Long>());
        }
        LinkedList<Long> val = b_ids.get(id);
        val.add(new Long(i));
      }
      PCLFileReader reader = new PCLFileReader(file1);
      PCLFileWriter writer = new PCLFileWriter(ofile);
      reader.begin();
      GeneData h1 = reader.getHeader();
      GeneData h2 = data.getHeader();
      GeneData h = GeneData.merge(h1, h2.subset(1,h2.size()-1));
      writer.writeData(h);
      while(reader.hasNext()) {
        GeneData gene = reader.getData();
        if (gene == null) {
            break;
        }
        gene = (GeneData) gene.clone();
        if ( (reader.getLineNumber() % 1000 == 0) ) {
            System.out.println(reader.getLineNumber());
        }
        String id = (String) gene.getDataAt(0);
        id = id.replaceAll("\\s", "");
        if (id.equals("EWEIGHT")) {
            GeneData g = GeneData.getWeight(0, h2.size()-1);
            g = GeneData.merge(gene, g.subset(1,g.size()-1));
            writer.writeData(g);
            continue;
        }
        if (b_ids.containsKey(id)) {
          LinkedList<Long> val = b_ids.get(id);
          Iterator<Long> itr = val.iterator();
          while (itr.hasNext()) {
            Long e = (Long) itr.next();
            GeneData g = data.getDataAt(e.longValue());
            g = GeneData.merge(gene, g.subset(1,g.size()-1));
            writer.writeData(g);
          }
        }
      }
      writer.close();
    }
  }

  public static void pieAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    String file2 = list.removeFirst();
    String file3 = list.removeFirst();
    Data data1 = PCLFileReader.readFile(file1);
    Data data2 = PCLFileReader.readFile(file2);
    Data data3 = PCLFileReader.readFile(file3);
    Data res;
    Data tmp = Data.unionData(data2, data3);
    tmp = Data.diffData(data1, tmp);
    res = Data.selectColumnFromData(tmp, 0);
    System.out.println(tmp.getNumGenes());
    tmp = Data.unionData(data1, data3);
    tmp = Data.diffData(data2, tmp);
    tmp = Data.selectColumnFromData(tmp, 0);
    res = Data.concatDataColumns(res, tmp);
    System.out.println(tmp.getNumGenes());
    tmp = Data.unionData(data1, data2);
    tmp = Data.diffData(data3, tmp);
    tmp = Data.selectColumnFromData(tmp, 0);
    res = Data.concatDataColumns(res, tmp);
    System.out.println(tmp.getNumGenes());
    tmp = Data.intersectData(data1, data3);
    tmp = Data.diffData(tmp, data2);
    tmp = Data.selectColumnFromData(tmp, 0);
    res = Data.concatDataColumns(res, tmp);
    System.out.println(tmp.getNumGenes());
    tmp = Data.intersectData(data1, data2);
    tmp = Data.diffData(tmp, data3);
    tmp = Data.selectColumnFromData(tmp, 0);
    res = Data.concatDataColumns(res, tmp);
    System.out.println(tmp.getNumGenes());
    tmp = Data.intersectData(data2, data3);
    tmp = Data.diffData(tmp, data1);
    tmp = Data.selectColumnFromData(tmp, 0);
    res = Data.concatDataColumns(res, tmp);
    System.out.println(tmp.getNumGenes());
    tmp = Data.intersectData(data2, data3);
    tmp = Data.intersectData(tmp, data1);
    tmp = Data.selectColumnFromData(tmp, 0);
    res = Data.concatDataColumns(res, tmp);
    System.out.println(tmp.getNumGenes());
    PCLFileWriter.writeFile(res, ofile);
  }

  public static void statsCorrAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String file1 = list.removeFirst();
    int num = Integer.parseInt(list.removeFirst());
    Data data1 = PCLFileReader.readFile(file1);

    Bimodal[] b = new Bimodal[data1.getNumGenes()];
    for (int i = data1.getNumGeneHeader(); i < data1.getNumRows(); i++) {
      // System.out.println("[ " + i + " ]---------------------->");
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v= data1.getGeneData(i).getVector(start, end);
      b[i - data1.getNumGeneHeader()] = new Bimodal(v);
    }

    int min = num;
    if (min > data1.getNumGenes()) {
      min = data1.getNumGenes();
    }
    int[] res = new int[21];
    int[] res1 = new int[21];
    for (int i = 0; i < 21; i++) {
        res[i] = 0;
        res1[i] = 0;
    }
    Permutation p = new Permutation(100);
    int[] perm1 = p.getRandomPermutation(data1.getNumGenes());
    int[] perm2 = p.getRandomPermutation(data1.getNumGenes());
    for (int i = 0; i < min; i++) {
      int index = perm1[i] + data1.getNumGeneHeader();
      GeneData gene = data1.getGeneData(index);
      int start = data1.getNumArrayHeader();
      int end = data1.getNumColumns() - 1;
      Double[] v1 = gene.getVector(start, end);
      double thr1 = b[perm1[i]].getThreshold();
      System.out.println(i);
      for (int j = 0; j < min; j++) {
        index = perm2[j] + data1.getNumGeneHeader();
        Double[] v2 = data1.getGeneData(index).getVector(start, end);
        double corr = getCorrelation(v1, v2);
        Double tmp = new Double(corr*10 + 10);
        res[tmp.intValue()] ++;
        double thr2 = b[perm2[j]].getThreshold();
        int c0 = countThreshold(0, thr1, thr2, v1, v2);
        int c1 = countThreshold(1, thr1, thr2, v1, v2);
        int c2 = countThreshold(2, thr1, thr2, v1, v2);
        int c3 = countThreshold(3, thr1, thr2, v1, v2);
        // System.out.println(thr1 + "," + thr2 + ":" + c0 + "-" +c1 + "-" +c2 + "-" +c3);
        int len = v1.length;
        int max = (c0 + c2);
        if (max < (c1 + c3)) {
            max = (c1 + c3);
        }
        tmp = new Double(max * 10.0/v1.length + 10);
        double p0 = 1 - GeneSetAnalysis.getPvalue(c0, c0+c1, c0+c3, len);
        tmp = new Double(p0 * 10.0 + 10);
        res1[tmp.intValue()] ++;
      }
    }
    BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
    for (int i = 0; i < 21; i++) {
        out.write( ((i-10.0)/10) + "\t" + res[i] + "\t" + res1[i] + "\n");
    }
    out.close();
  }

  public static void statsAnalysis(LinkedList<String> list) throws Exception {
    String cmd = list.removeFirst();
    if (cmd.equals("corr")) {
        statsCorrAnalysis(list);
    }
    else {
      System.out.println("<cmd> <args> :");
      System.out.println("stats             corr ofile pclfile num");
    }
  }

  public static void stepAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String pclfile = list.removeFirst();
    double pvalue = 0.05;
    double goPvalue = 0.001;
    boolean fdr = false;

    String file = pclfile;
    String type = "BothStep";
    String centering = "Step";
    String range = null;
    String org = "Hs";
    String annFile = "gene_association.goa_human";
    String onnFile = "gene_ontology.obo";
    String tStr = null; // timepoints specification
    int geneIndex = 0;
    int splitIndex = -1;
    String splitString = ":";
    int numMissing = 0;

    while (list.size() > 0) {
        String cmd = list.removeFirst();
        String arg;
        if (list.size() > 0) {
            arg = list.removeFirst();
        }
        else {
            break;
        }
        if (cmd.equals("type")) { type = arg; }
        if (cmd.equals("centering")) { centering = arg; }
        if (cmd.equals("range")) { range = arg; }
        if (cmd.equals("org")) { org = arg; }
        if (cmd.equals("annFile")) { annFile = arg; }
        if (cmd.equals("onnFile")) { onnFile = arg; }
        if (cmd.equals("timepoints")) { tStr = arg; }
        if (cmd.equals("pvalue")) { pvalue = Double.parseDouble(arg); }
        if (cmd.equals("goPvalue")) { goPvalue = Double.parseDouble(arg); }
        if (cmd.equals("fdr")) { fdr = Boolean.parseBoolean(arg); }
        if (cmd.equals("geneIndex")) { geneIndex = Integer.parseInt(arg); }
        if (cmd.equals("splitIndex")) { splitIndex = Integer.parseInt(arg); }
        if (cmd.equals("numMissing")) { numMissing = Integer.parseInt(arg); }
        if (cmd.equals("splitString")) { splitString = arg; }
    }

    GeneNameScheme ns = new GeneNameScheme(geneIndex,
        org, splitString, splitIndex);
    ns.setAnnotationFile(annFile);
    ns.setOntologyFile(onnFile);
    ns.setNumMissingPoints(numMissing);
    Data data = PCLFileReader.readFile(pclfile);
    data.setGeneNameScheme(ns);
    data.setRange(range);
    data.convertDoubles();
    int numAH = data.getNumArrayHeader();
    GeneData header = data.getGeneData(0);

    Double[] timepoints = new Double[data.getNumArrays()];
    for (int i =0; i < timepoints.length; i++) {
      timepoints[i] = new Double(i);
    }
    if (tStr != null) {  // example 0-5x2,1,2,3,0,1-5,0
      String[] ilist = tStr.split(",");
      int index = 0;
      for (int i =0; i < ilist.length ; i++) {
        String[] mlist = ilist[i].split("x");
        int times = 1;
        if (mlist.length > 1) {
           times = Integer.parseInt(mlist[1]);
        }
        String[] vlist = mlist[0].split("-");
        int start = Integer.parseInt(vlist[0]);
        int end = start;
        if (vlist.length > 1) {
           end = Integer.parseInt(vlist[1]);
        }
        for (int j = start; j <= end ; j++) {
          for (int k = 0; k < times ; k++) {
             if (index >= timepoints.length) {
                break;
             }
             timepoints[index++] = new Double(j);
          }
          if (index > timepoints.length) {
            break;
          }
        }
        if (index > timepoints.length) {
          break;
        }
      }
      for (int i =0; i < timepoints.length; i++) {
        System.out.println(timepoints[i] + " " + header.getDataAt(i+numAH));
      }
    }
    data.setTimepoints(timepoints);
    StepMiner sm = new StepMiner(data);
    if (type.equals("OneStep")) {
      sm.setOneStepAnalysis();
    }
    if (type.equals("BothStep")) {
      sm.setBothStepAnalysis();
    }
    if (type.equals("TwoStep")) {
      sm.setTwoStepAnalysis();
    }
    if (type.equals("SelectTwoStep")) {
      sm.setSelectTwoStepAnalysis();
    }
    if (centering.equals("NoCentering")) {
      sm.setStepCentering(false);
    }
    if (centering.equals("Step")) {
      sm.setStepCentering(true);
    }
    sm.setFdrAnalysis(fdr);
    sm.setPvalueThr(pvalue);
    sm.performAnalysis();
    String f = ofile;
    if (f.endsWith(".pcl")) {
      sm.writePCL(f);
    }
    if (f.endsWith(".ano")) {
      PCLFileWriter.writeFile(sm.getStepAnnotationData(), f);
    }
    if (f.endsWith(".ing")) {
      PCLFileWriter.writeFile(sm.getStepAnnotationByGenes(), f);
    }
    if (f.endsWith(".cdt")) {
      PCLFileWriter.writeFile(sm.getCdtAnnotationData(), f);
      String gtr = f.replaceFirst(".cdt$", ".gtr");
      PCLFileWriter.writeFile(sm.getGtrAnnotationData(), gtr);
    }
    if (f.endsWith(".html")) {
      sm.performGOAnalysis(f, goPvalue);
    }
    if (f.endsWith(".htm")) {
      sm.performGOBestAnalysis(f, goPvalue);
    }
    if (f.endsWith(".ps") || f.endsWith(".eps")) {
      sm.plotSteps(f);
    }
    if (f.endsWith(".ann")) {
      sm.writeAnnotations(f);
    }
    if (f.endsWith(".gmt") || f.endsWith(".gxa") || f.endsWith(".tab")) {
      sm.writeGeneSets(f);
    }
  }

  public static void twoStepAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String pclfile = list.removeFirst();
    int index = Integer.parseInt(list.removeFirst());
    double pvalue = Double.parseDouble(list.removeFirst());
    Data data = PCLFileReader.readFile(pclfile);
    data.convertDoubles();
    int start = data.getNumArrayHeader();
    int end = data.getNumColumns() - 1;

    Double[] timepoints = new Double[data.getNumArrays()];
    int breakpoint = index - data.getNumArrayHeader();
    for (int i =0; i < timepoints.length; i++) {
      timepoints[i] = new Double(i/3);
    }
    for (int i =0; i < timepoints.length; i++) {
        System.out.println(timepoints[i]);
    }
    data.setTimepoints(timepoints);
    data.setBreakpoint(index);
    StepMiner sm = new StepMiner(data);
    //sm.setTwoStepAnalysis();
    sm.setBothStepAnalysis();
    //sm.setSelectTwoStepAnalysis();
    sm.setFdrAnalysis(true);
    sm.setPvalueThr(pvalue);
    sm.performAnalysis();
    sm.writePCL(ofile);

    /*
    data.setRange(start+":"+(index-1));
    StepMiner sm = new StepMiner(data);
    sm.setStepCentering(false);
    sm.setOneStepAnalysis();
    sm.setPvalueThr(pvalue);
    sm.performAnalysis();
    data.setRange(start+":"+end);
    Data d1 = sm.getStepOrderedData("Up");
    d1.setRange(index+":"+end);
    StepMiner sm1 = new StepMiner(d1);
    sm1.setStepCentering(false);
    sm1.setOneStepAnalysis();
    sm1.setPvalueThr(pvalue);
    sm1.performAnalysis();
    d1.setRange(start+":"+end);
    d1 = sm1.getStepOrderedData();
    Data d2 = sm.getStepOrderedData("Down");
    d2.setRange(index+":"+end);
    StepMiner sm2 = new StepMiner(d2);
    sm2.setStepCentering(false);
    sm2.setOneStepAnalysis();
    sm2.setPvalueThr(pvalue);
    sm2.performAnalysis();
    d2.setRange(start+":"+end);
    d2 = sm2.getStepOrderedData();
    Data res = Data.concatData(d1, d2);
    ArrayOrder arr = new ArrayOrder(res);
    arr.meanCenter();
    res = arr.getOrderedData();
    PCLFileWriter.writeFile(res, ofile);
    */
  }

  public static void twoStepAnalysis1(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String pclfile = list.removeFirst();
    int index = Integer.parseInt(list.removeFirst());
    double pvalue = Double.parseDouble(list.removeFirst());
    Data data = PCLFileReader.readFile(pclfile);
    data.convertDoubles();
    StepMiner sm = new StepMiner(data);
    sm.setTwoStepAnalysis();
    sm.setFdrAnalysis(false);
    sm.setPvalueThr(pvalue);
    sm.performAnalysis();
    Data dataAnn = sm.getStepAnnotationData();
    int start = dataAnn.getNumArrayHeader();
    int end = dataAnn.getNumColumns() - 1;
    Vector<GeneData> v = new Vector<GeneData>();
    for (int i = dataAnn.getNumGeneHeader(); i < dataAnn.getNumRows(); i++) {
      // System.out.println("[ " + i + " ]---------------------->");
      GeneData gd = dataAnn.getGeneData(i);
      int step1 = Integer.parseInt((String)gd.getDataAt(start-3));
      int step2 = step1 + Integer.parseInt((String)gd.getDataAt(start-2));
      if (step1 < index && step2 > index) {
        //System.out.println(gd.getDataAt(0) + "\t" + step1 + "\t" + step2);
        v.add(gd.subset(0,0));
      }
    }
    Data res = new Data(0, v.size(), 0, 0, v);
    res = Data.intersectDataSelect(sm.getStepOrderedData(), res);
    PCLFileWriter.writeFile(res, ofile);
  }

  public static void corrAllAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String listfile = list.removeFirst();
    String pclfile = list.removeFirst();
    PCLFileReader data = new PCLFileReader(listfile);
    PCLFileReader.CACHE_SIZE = 0;
    HashSet<String> b_ids = new HashSet<String>();
    data.begin();
    for (int i = 0; data.hasNext(); i++) {
      GeneData gene = data.getData();
      if (gene == null) {
        break;
      }
      if ( (i % 1000) == 0) {
        System.gc();
        System.out.println("[" + i + "]");
      }
      String id = (String) gene.getDataAt(0);
      id = id.replaceAll("\\s", "");
      id = new String(id);
      if (!b_ids.contains(id)) {
        b_ids.add(id);
      }
    }
    HashSet<Integer> pcl_ids = new HashSet<Integer>();
    PCLFileReader reader = new PCLFileReader(pclfile);
    reader.beginRandomAccess();
    for (int i = 0; reader.hasNext(); i++) {
      GeneData gene = reader.getDataAt(-1);
      if (gene == null) {
        break;
      }
      if ( (i % 1000) == 0) {
        System.gc();
        System.out.println("[" + i + "]");
      }
      String id = (String) gene.getDataAt(0);
      id = id.replaceAll("\\s", "");
      id = new String(id);
      if (b_ids.contains(id)) {
        pcl_ids.add(new Integer(i));
      }
    }
    PCLFileWriter writer = new PCLFileWriter(ofile);
    int start = reader.getNumArrayHeader();
    int end = reader.getNumColumns() - 1;
    Iterator<Integer> itr = pcl_ids.iterator();
    while (itr.hasNext()) {
      Integer a = (Integer) itr.next();
      GeneData ga = reader.getDataAt(a.intValue());
      Double[] v1 = ga.getVector(start, end);
      Iterator<Integer> itr1 = pcl_ids.iterator();
      while (itr1.hasNext()) {
        Integer b = (Integer) itr1.next();
        if (a.intValue() > b.intValue()) {
          GeneData gb = reader.getDataAt(b.intValue());
          Double[] v2 = gb.getVector(start, end);
          double corr = getCorrelation(v1, v2);
          Object[] d = new Object[5];
          String sr = formatString("0.####", corr);
          d[0] = sr;
          d[1] = ga.getDataAt(0);
          d[2] = gb.getDataAt(0);
          d[3] = ga.getDataAt(1);
          d[4] = gb.getDataAt(1);
          GeneData gg = new GeneData(d);
          writer.writeData(gg);
        }
      }
    }
    writer.close();
  }

  public static void corrListAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String listfile = list.removeFirst();
    String pclfile = list.removeFirst();
    PCLFileReader data = new PCLFileReader(listfile);
    PCLFileReader.CACHE_SIZE = 1000;
    data.begin();
    PCLFileReader reader = new PCLFileReader(pclfile);
    reader.beginRandomAccess();
    PCLFileWriter writer = new PCLFileWriter(ofile);
    int start = reader.getNumArrayHeader();
    int end = reader.getNumColumns() - 1;
    int num = data.getNumColumns();
    for (int i = 0; data.hasNext(); i++) {
      GeneData gene = data.getData();
      if (gene == null) {
        break;
      }
      if ( (i % 1000) == 0) {
        System.gc();
        System.out.println("[" + i + "]");
      }
      int a = Integer.parseInt((String) gene.getDataAt(1));
      int b = Integer.parseInt((String) gene.getDataAt(2));
      GeneData ga = reader.getDataAt(a);
      Double[] v1 = ga.getVector(start, end);
      GeneData gb = reader.getDataAt(b);
      Double[] v2 = gb.getVector(start, end);
      double corr = getCorrelation(v1, v2);
      Object[] d = new Object[num+1];
      String sr = formatString("0.####", corr);
      for (int j =0; j < num; j++) {
        d[j] = gene.getDataAt(j);
      }
      d[num] = sr;
      GeneData gg = new GeneData(d);
      writer.writeData(gg);
    }
    writer.close();
  }

  public static void corrOneAnalysis(LinkedList<String> list) throws Exception {
    String ofile = list.removeFirst();
    String pclfile = list.removeFirst();
    int n = Integer.parseInt(list.removeFirst());
    int k = -1;
    if (list.size() > 0) {
      k = Integer.parseInt(list.removeFirst());
    }
    PCLFileReader data = new PCLFileReader(pclfile);
    PCLFileReader.CACHE_SIZE = 2000;
    data.beginRandomAccess();
    GeneData gene = data.getDataAt(n);
    int start = data.getNumArrayHeader();
    if (list.size() > 0) {
      start = Integer.parseInt(list.removeFirst());
    }
    int end = data.getNumColumns() - 1;
    PCLFileWriter writer = new PCLFileWriter(ofile);
    Object[] d = new Object[3];
    d[0] = "1.0000"; d[1] = gene.getDataAt(0); d[2] = gene.getDataAt(1);
    GeneData gg = new GeneData(d);
    writer.writeData(gg);
    Double[] v1 = gene.getVector(start, end);
    if (k < 0) {
      GeneData ga = data.getDataAt(0);
      for (int i = 0;  ga != null; i++) {
        String id = (String) ga.getDataAt(0);
        Double[] v2 = ga.getVector(start, end);
        double corr = getCorrelation(v1, v2);
        String sr = formatString("0.####", corr);
        d[0] = sr; d[1] = ga.getDataAt(0); d[2] = ga.getDataAt(1);
        gg = new GeneData(d);
        writer.writeData(gg);
        ga = data.getDataAt(i+1);
      }
    }
    else {
      GeneData ga = data.getDataAt(k);
      String id = (String) ga.getDataAt(0);
      Double[] v2 = ga.getVector(start, end);
      double corr = getCorrelation(v1, v2);
      String sr = formatString("0.####", corr);
      d[0] = sr; d[1] = ga.getDataAt(0); d[2] = ga.getDataAt(1);
      gg = new GeneData(d);
      writer.writeData(gg);
    }
    writer.close();
  }

  public static void booleanAnalysis(LinkedList<String> list) throws Exception {
    String cmd = list.removeFirst();
    if (cmd.equals("bitMatrix") || cmd.equals("listMatrix") ||
        cmd.equals("pairs") || cmd.equals("listPairs") ||
        cmd.equals("singleListMatrix") || cmd.equals("listMatrixDebug")) {  
      String ofile = list.removeFirst();
      String bvfile = list.removeFirst();
      String phfile = list.removeFirst();
      String phid = list.removeFirst();
      BooleanAnalysis ana = new BooleanAnalysis(bvfile, ofile, phfile, phid);
      double thr = Double.parseDouble(list.removeFirst());
      ana.setThreshold(thr);
      double statThr = Double.parseDouble(list.removeFirst());
      ana.setStatThreshold(statThr);
      double singleThr = Double.parseDouble(list.removeFirst());
      ana.setSingleThreshold(singleThr);
      if (cmd.equals("listMatrix")) {
        ana.setGeneList(list.removeFirst());
        ana.performListAnalysis();
      }
      if (cmd.equals("listMatrixDebug")) {
        ana.setGeneList(list.removeFirst());
        ana.performListDebugAnalysis();
      }
      if (cmd.equals("singleListMatrix")) {
        ana.setGeneList(list.removeFirst());
        ana.performSingleListAnalysis();
      }
      if (cmd.equals("listPairs")) {
        ana.setGeneList(list.removeFirst());
        ana.writeListPairs();
      }
      if (cmd.equals("bitMatrix")) {
        ana.performAnalysis();
      }
      if (cmd.equals("pairs")) {
        ana.writePairs();
      }
    }
    if (cmd.equals("commonPairs")) {
      String ofile = list.removeFirst();
      String bvfile = list.removeFirst();
      BooleanAnalysis ana = new BooleanAnalysis(bvfile, ofile);
      double thr = Double.parseDouble(list.removeFirst());
      ana.setThreshold(thr);
      double statThr = Double.parseDouble(list.removeFirst());
      ana.setStatThreshold(statThr);
      ana.writeCommonPairs(list);
    }
    if (cmd.equals("printNumbers")) {
      String ofile = list.removeFirst();
      String bvfile = list.removeFirst();
      BooleanAnalysis ana = new BooleanAnalysis(bvfile, ofile);
      double thr = Double.parseDouble(list.removeFirst());
      ana.setThreshold(thr);
      double statThr = Double.parseDouble(list.removeFirst());
      ana.setStatThreshold(statThr);
      PCLFileReader reader = new PCLFileReader(bvfile);
      PCLFileReader.CACHE_SIZE = 0;
      reader.beginRandomAccess();
      int a = Integer.parseInt(list.removeFirst());
      int b = Integer.parseInt(list.removeFirst());
      GeneData ga = reader.getDataAt(a);
      BitSet va = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 0);
      BitSet va_thr = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 1);
      GeneData gb = reader.getDataAt(b);
      BitSet vb = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 0);
      BitSet vb_thr = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 1);
      double[] p = ana.getErrorProbStats(va, va_thr, vb, vb_thr, 1);
    }
    if (cmd.equals("bitMatrixPrint")) {
      String filename = list.removeFirst();
      BitMatrixNetworkSimple file = new BitMatrixNetworkSimple(filename,
          BitMatrixFile.READMODE);
      file.readMatrixFile();
    }
    if (cmd.equals("bitMatrixPrintStats")) {
      String filename = list.removeFirst();
      BitMatrixNetworkSimple file = new BitMatrixNetworkSimple(filename,
          BitMatrixFile.READMODE);
      file.printStats();
    }
    if (cmd.equals("bitMatrixFill")) {
      BitMatrixFile.BLOCKSIZE = 50000;
      BitMatrixFile.CACHESIZE = 50000;
      String filename = list.removeFirst();
      BitMatrixNetworkSimple file = new BitMatrixNetworkSimple(filename,
          BitMatrixFile.WRITEMODE);
      file.readMatrixHeader();
      file.fillLowerTriangle();
      file.close();
    }
    if (cmd.equals("bitMatrixFillStats")) {
      BitMatrixFile.BLOCKSIZE = 50000;
      BitMatrixFile.CACHESIZE = 50000;
      String filename = list.removeFirst();
      BitMatrixNetworkSimple file = new BitMatrixNetworkSimple(filename,
          BitMatrixFile.WRITEMODE);
      file.readMatrixHeader();
      file.fillStats();
      file.close();
    }
  }

  public static void shuffleAnalysis(LinkedList<String> list) throws Exception {
    String cmd = list.removeFirst();
    if (cmd.equals("bv")) {
      String ofile = list.removeFirst();
      String bvfile = list.removeFirst();
      int seed = Integer.parseInt(list.removeFirst());
      PCLFileReader reader = new PCLFileReader(bvfile);
      PCLFileReader.CACHE_SIZE = 0;
      reader.begin();
      BufferedWriter out = new BufferedWriter(new FileWriter(ofile));
      GeneData header = reader.getHeader();
      out.write(header.toString());
      Permutation p = new Permutation(seed);
      for (int i = 1; reader.hasNext(); i++) {
        GeneData ga = reader.getData();
        if (ga == null) {
          break;
        }
        if ( (i%100) == 0) {
            System.err.println("[" + i + "]");
        }
        String str = (String) ga.getDataAt(2);
        byte[] arr = str.getBytes();
        arr = p.permute(arr, p.getRandomPermutation(arr.length));
        str = new String(arr);
        ga.setDataAt(2, str);
        out.write(ga.toString());
      }
      out.close();
    }
  }

  public static void htAnalysis(LinkedList<String> list) throws Exception {
    // Hypergeometric test
    int N = Integer.parseInt(list.removeFirst());
    int M = Integer.parseInt(list.removeFirst());
    int n = Integer.parseInt(list.removeFirst());
    int k = Integer.parseInt(list.removeFirst());
    double pval = GeneSetAnalysis.getPvalue(k, n, M, N);
    String sr1 = formatString("0.####E0", pval);
    String sr2 = formatString("0.####", pval);
    System.out.println(sr1 + "\t" + sr2);
  }

  public static void testAnalysis(LinkedList<String> list) throws Exception {
    /*
    // Hypergeometric test
    int N = Integer.parseInt(list.removeFirst());
    int M = Integer.parseInt(list.removeFirst());
    int n = Integer.parseInt(list.removeFirst());
    int k = Integer.parseInt(list.removeFirst());
    double pval = 1 - GeneSetAnalysis.getPvalue(k, n, M, N);
    String sr = formatString("0.####", pval);
    System.out.println(sr);
    */
    /*
    // Find correlation of two given indices
    String file = list.removeFirst();
    int n = Integer.parseInt(list.removeFirst());
    int k = Integer.parseInt(list.removeFirst());
    PCLFileReader data = new PCLFileReader(file);
    PCLFileReader.CACHE_SIZE = 2000;
    data.beginRandomAccess();
    GeneData gene = data.getDataAt(n);
    int start = data.getNumArrayHeader();
    int end = data.getNumColumns() - 1;
    System.out.println(gene.getDataAt(1));
    Double[] v1 = gene.getVector(start, end);
    GeneData ga = data.getDataAt(0);
    for (int i = 0;  ga != null; i++) {
      String id = (String) ga.getDataAt(0);
      Double[] v2 = ga.getVector(start, end);
      double corr = getCorrelation(v1, v2);
      String sr = formatString("0.####", corr);
      System.out.print(sr + "\t");
      System.out.println(ga.getDataAt(0) + "\t" + ga.getDataAt(1));
      ga = data.getDataAt(i+1);
    }
    */
    String file = list.removeFirst();
    PCLFileReader data = new PCLFileReader(file);
    PCLFileReader.CACHE_SIZE = 2000;
    data.beginRandomAccess();
    GeneData ga = data.getDataAt(0);
    for (int i = 0;  ga != null; i++) {
      BitSet va_thr = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 1);
      int outside = va_thr.cardinality();
      System.out.println(ga.getDataAt(0) + "\t" + outside);
      ga = data.getDataAt(i+1);
    }
    /*
    String ofile = list.removeFirst();
    Data res = TABFileReader.readFile(list.removeFirst());
    res.setRange("1");
    res.setRowRange("1");
    res.addWeightRow(1);
    PCLFileWriter.writeFile(res, ofile);
    */
  }

  public static void pln(String str) {
    System.out.println(str);
  }

  public static void prn(String str) {
    System.out.print(str);
  }

  public static void main(String[] args) throws Exception {
    LinkedList<String> list = new LinkedList<String>(Arrays.asList(args));
    String cmd = "";
    if (list.size() > 0) {
      cmd = list.removeFirst();
    }
    if (cmd.equals("intersect")) {
      intersectionAnalysis(list);
    }
    else if(cmd.equals("intersectGenes")) {
      intersectionGenesAnalysis(list);
    }
    else if(cmd.equals("intersectStep")) {
      intersectStepAnalysis(list);
    }
    else if(cmd.equals("intersectStepAll")) {
      intersectStepAllAnalysis(list);
    }
    else if(cmd.equals("intersectStepGO")) {
      intersectStepGOAnalysis(list);
    }
    else if(cmd.equals("intersectStepGS")) {
      intersectStepGSAnalysis(list);
    }
    else if(cmd.equals("foldChange")) {
      countTotalFoldChange(list);
    }
    else if(cmd.equals("corr")) {
      corrAnalysis(list);
    }
    else if(cmd.equals("corrStep")) {
      corrStepAnalysis(list);
    }
    else if(cmd.equals("monotonic")) {
      monotonicAnalysis(list);
    }
    else if(cmd.equals("aracne")) {
      aracneAnalysis(list);
    }
    else if(cmd.equals("geo")) {
      geoAnalysis(list);
    }
    else if(cmd.equals("geoSampleIds")) {
      geoAnalysis1(list);
    }
    else if(cmd.equals("gs")) {
      GSAnalysis(list);
    }
    else if(cmd.equals("go")) {
      GOAnalysis(list);
    }
    else if(cmd.equals("bimodal")) {
      bimodalAnalysis(list);
    }
    else if(cmd.equals("removeArray")) {
      selectArrayPCLAnalysis(true, list);
    }
    else if(cmd.equals("selectArrayPCL")) {
      selectArrayPCLAnalysis(false, list);
    }
    else if(cmd.equals("selectArray")) {
      selectArrayAnalysis(list);
    }
    else if(cmd.equals("filter")) {
      filterAnalysis(list);
    }
    else if(cmd.equals("tabFile")) {
      tabFileAnalysis(list);
    }
    else if(cmd.equals("pie")) {
      pieAnalysis(list);
    }
    else if(cmd.equals("stats")) {
      statsAnalysis(list);
    }
    else if(cmd.equals("step")) {
      stepAnalysis(list);
    }
    else if(cmd.equals("twoStep")) {
      twoStepAnalysis(list);
    }
    else if(cmd.equals("twoStep1")) {
      twoStepAnalysis1(list);
    }
    else if(cmd.equals("corrAll")) {
      corrAllAnalysis(list);
    }
    else if(cmd.equals("corrList")) {
      corrListAnalysis(list);
    }
    else if(cmd.equals("corrOne")) {
      corrOneAnalysis(list);
    }
    else if(cmd.equals("shuffle")) {
      shuffleAnalysis(list);
    }
    else if(cmd.equals("boolean")) {
      booleanAnalysis(list);
    }
    else if(cmd.equals("ht")) {
      htAnalysis(list);
    }
    else if(cmd.equals("test")) {
      testAnalysis(list);
    }
    else {
      pln("tools.CustomAnalysis <cmd> <args>");
      pln("<cmd> <args> :");
      pln("intersect         file1 file2 file3 ...");
      pln("intersectGenes    file1 file2 file3 ...");
      pln("intersectStep     spvalue file1 file2 file3 ...");
      pln("intersectStepAll  spvalue outfile file1 file2 file3 ...");
      pln("intersectStepGO   spvalue outfile onnFile annFile org pvalue file1 file2 ...");
      pln("intersectStepGS   spvalue outfile setFile org pvalue file1 file2 ...");
      pln("corr              outfile threshold file1 file2");
      pln("corrAll           outfile listIds pcl");
      pln("corrList          outfile pairsFile pcl");
      pln("corrOne           outfile pclFile num1 [num2 [start]]");
      pln("corrStep          outfile spvalue threshold file1 file2 file3");
      pln("monotonic         up/down outfile threshold file");
      pln("aracne            outfile num pclfile pairfile");
      pln("geo               outfile outInfoFile <Soft GEO file> <numArr> <Random:true/false> [symindex=8 titleindex=7] [<Exclusion list>]");
      pln("geoSampleIds      outfile <Soft GEO file>");
      pln("bimodal           <get[r]> outfile pclfile threshold");
      pln("                  getr -> reduce log data file");
      pln("bimodal           <thr[1|Bv]/pairs[Bv]/groups> outfile pclfile");
      pln("bimodal           thrBv1 outfile pclfile thrfile");
      pln("bimodal           targets outfile pclfile geneid");
      pln("bimodal           plot outfile pclfile pairfile [<info file>]");
      pln("bimodal           commonPairs outfile pclfile pairfile bitVectorFile");
      pln("bimodal           pairsBvParallel outfile BvFile start1 num1 start2 num2");
      pln("gs                outfile org pvalue file1 [file2]");
      pln("go                outfile onnFile annFile org pvalue gmtfile");
      pln("removeArray       outfile pclFile listFile");
      pln("selectArrayPCL    outfile pclFile listFile");
      pln("selectArray       outfile pclFile geneid minthreshold maxthreshold");
      pln("filter            reduceLog outfile pclFile");
      pln("filter            normalize outfile pclFile geneid");
      pln("filter            selectNames outfile pclFile nameFile");
      pln("filter            add outfile pclFile constant");
      pln("tabFile           concat outfile file1 file2");
      pln("tabFile           select outfile file1 range");
      pln("tabFile           delete outfile file1 range");
      pln("tabFile           intersect[R] outfile file1 file2large");
      pln("pie               ofile file1 file2 file3");
      pln("stats             corr ofile pclfile num");
      pln("step              ofile pclfile [<cmd> <arg>]");
      pln("                  <cmd>: timepoints, type, range, org, geneIndex,");
      pln("                         splitIndex, splitString, pvalue");
      pln("twoStep           ofile pclfile index pvalue");
      pln("twoStep1          ofile pclfile index pvalue");
      pln("shuffle           bv <outfile> <bvfile> <seed>");
      pln("                  pcl <outfile> <pclfile> <seed>");
      pln("boolean           <cmd> args");
      pln("                  bitMatrix/pairs/listMatrix/singleListMatrix ofile bvfile phfile phid pvalue statThr singleThr [listFile]");
      pln("                  commonPairs ofile old_bvfile pvalue statThr old_pairfile new_bvfile [HomologFile oldOrg newOrg]");
      pln("                  bitMatrixPrint <relationFile>");
      pln("                  bitMatrixPrintStats <relationFile>");
      pln("                  bitMatrixFill <relationFile>");
      pln("                  bitMatrixFillStats <relationFile>");
      pln("ht                N M n k");
      pln("test              <args>");
    }
  }
}
