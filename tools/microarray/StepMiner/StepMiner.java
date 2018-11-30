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

import tools.microarray.*;
import tools.microarray.FileWriter.*;
import java.util.*;
import java.io.*;
import java.text.MessageFormat;
import tools.graphs.*;

public class StepMiner implements Cloneable {

  Data data_;
  boolean fdrAnalysis_;
  int numRandom_;
  double pvalueThr_;
  int numMissingPoints_;
  int analysisType_;  // Analysis type : Only one step vs both steps.
                      //  0 - only one step, 1 - both step, 2 - only two step
                      //  3 - both step - select two step only
  boolean stepCentering_; // Whether to perform step centering
  boolean stopAnalysis_;  // Don't fit steps; Can't perform centering;
  boolean deleteGenesInMultipleGroups_;  // delete genes that appear in other
                                         // groups
  int start_;    // Range for the analysis
  int end_;      // Range for the analysis
  int[] genes_;  // Filtered genes (index to the original data)
  Step[] steps_; // BestStep  : All the FittedSteps
  Step[] bestSteps_; // Only the best step

  String fdrStats;

  int[] sortedOrder_;  // Sorted order on the data_
  int[] sortedGenes_;  // Sorted genes (index into original data)
  Step[] sortedSteps_; // Only the best step

  int nodeNumber_; // tmp var for GtrAnnotation

  public StepMiner(Data data) {
    fdrAnalysis_ = false;
    numRandom_ = 100;
    pvalueThr_ = 0.05;
    analysisType_ = 1;
    stepCentering_ = true;
    stopAnalysis_ = false;
    deleteGenesInMultipleGroups_ = false;
    data_ = data;
    init_data_();
    fdrStats = "No Fdr";
  }

  void init_data_() {
    numMissingPoints_ = data_.getNumMissingPoints();
    start_ = data_.getNumArrayHeader();
    end_ = data_.getNumColumns() - 1;
  }

  public void setPvalueThr(double p) {pvalueThr_ = p;}
  public void setOneStepAnalysis() { analysisType_ = 0;}
  public void setBothStepAnalysis() { analysisType_ = 1;}
  public void setTwoStepAnalysis() { analysisType_ = 2;}
  public void setSelectTwoStepAnalysis() { analysisType_ = 3;}
  public void setFdrAnalysis(boolean f) { fdrAnalysis_ = f;}
  public void setStepCentering(boolean v) { stepCentering_ = v;}
  public void setGeneNameScheme(GeneNameScheme s){data_.setGeneNameScheme(s);}
  public void setDeleteGenesInMultipleGroups(boolean v) 
  { deleteGenesInMultipleGroups_ = v; }

  public int[] getGeneOrder() { return genes_; }
  public int[] getSortedGeneOrder() { return sortedGenes_; }
  public String getFdrStats() { return fdrStats; }

  public double[] getP1() throws StepException {
    double[] res = new double[steps_.length];
    for (int i =0; i < steps_.length; i++) {
        BestStep s = (BestStep) steps_[i];
        Step step = s.findStep(1); // get oneStep
        res[i] = step.getPvalue();
    }
    return res;
  }
  public double[] getP2() throws StepException {
    double[] res = new double[steps_.length];
    for (int i =0; i < steps_.length; i++) {
        BestStep s = (BestStep) steps_[i];
        Step step = s.findStep(2); // get twoStep
        res[i] = step.getPvalue();
    }
    return res;
  }
  public double[] getP12() throws StepException {
    double[] res = new double[steps_.length];
    for (int i =0; i < steps_.length; i++) {
        BestStep s = (BestStep) steps_[i];
        res[i] = s.getP12();
    }
    return res;
  }
  public FdrInfo getFdrInfo() throws StepException {
    return new FdrInfo(getP1(), getP2(), getP12());
  }

  protected Object clone() throws CloneNotSupportedException {
    StepMiner res = new StepMiner(data_);
    res.pvalueThr_        = pvalueThr_          ;
    res.start_            = start_              ;
    res.end_              = end_                ;
    res.genes_            = genes_              ;
    return res;
  }

  /*
   * Filter genes based on number of missing points.
   */
  public void performFiltering() throws StepException {
    LinkedList<Integer> list = new LinkedList<Integer>();
    for (int i = data_.getNumGeneHeader(); i < data_.getNumRows(); i++) {
      GeneData gene = data_.getGeneData(i);
      if (numMissingPoints_ < 0 ||
        gene.getMissingPoints(start_, end_) <= numMissingPoints_) {
        list.add(new Integer(i));
      }
    }
    genes_ = new int[list.size()];
    ListIterator iter = list.listIterator();
    int i = 0;
    while (iter.hasNext()) {
        Integer result = (Integer) iter.next();
        genes_[i] = result.intValue();
        i++;
    }
    System.out.println("Filtering Done: " + data_.getNumGenes() + " -> " + 
        genes_.length);
  }

  public int countLabel(int label) {
    int count = 0;
    for (int i =0; i < bestSteps_.length ; i++) {
      if (bestSteps_[i].getLabel() == label) {
        count++;
      }
    }
    return count;
  }

  public int countLabelSorted(int label) {
    int count = 0;
    if (sortedSteps_ != null) {
      for (int i =0; i < sortedSteps_.length ; i++) {
        if (sortedSteps_[i].getLabel() == label) {
          count++;
        }
      }
    }
    return count;
  }

  public String getStats() throws StepException {
    int deleted = 0;
    deleted = countLabel(0) + countLabel(1) + countLabel(2) + countLabel(3) + countLabel(4) 
        - countLabelSorted(0) - countLabelSorted(1) - countLabelSorted(2) - countLabelSorted(3) - countLabelSorted(4);
            
    String res = ("One Step (Up)       = " + countLabelSorted(1)) + "\n" +
    ("One Step (Down)     = " + countLabelSorted(2)) + "\n" +
    ("Two Step (Up-Down)  = " + countLabelSorted(3)) + "\n" +
    ("Two Step (Down-Up)  = " + countLabelSorted(4)) + "\n" +
    ("Rest                = " + countLabelSorted(0)) + "\n" +
    ("Deleted             = " + deleted) + "\n" +
    ("Total               = " + bestSteps_.length + "x" + 
                (end_ - start_ + 1)) + "\n";
    return res;
  }

  public void printStepStats() throws StepException {
    HashMap<String, Integer> map = new HashMap<String, Integer>();
    for (int i = 0; i < sortedSteps_.length; i++) {
      Step step = sortedSteps_[i];
      int num = step.getNumSteps();
      int label = step.getLabel(); 
      // label : 0 - rest, 1 - Up, 2 -Down, 3 - UpDown, 4 - DownUp
      String tag = label +":";
      if (num == 1) {
        tag += step.getStep(0);
      }
      if (num == 2) {
        tag += step.getStep(0) + "-" + step.getStep(1);
      }
      Integer count = new Integer(1);
      if (map.containsKey(tag)) {
        count = map.get(tag);
        count = new Integer(count.intValue()+1);
      }
      map.put(tag, count);
    } // end for
    int total = 0;
    for (int i =start_; i <= end_; i++) {
      int step = i - start_;
      String uptag = "1:" + i;
      String downtag = "2:" + i;
      Integer count = new Integer(0);
      if (map.containsKey(uptag)) {
        count = map.get(uptag);
      }
      total += count.intValue();
      System.out.println("Up\t"+ step + "\t" + count);
      count = new Integer(0);
      if (map.containsKey(downtag)) {
        count = map.get(downtag);
      }
      total += count.intValue();
      System.out.println("Down\t"+ step + "\t" + count);
    }
    System.out.println("Total\t"+ 0 + "\t" + total);
  }


  public void printStats() throws StepException {
    System.out.print(getStats());
    /*
    printStepStats();
    System.out.println("ssr ");
    for (int i =0; i < steps_.length ; i++) {
        BestStep s = (BestStep) steps_[i];
        double p1 = s.findStep(1).getSsr();
        double p2 = s.findStep(2).getSsr();
        double p12 =  s.findStep(0).getSsr();
        System.out.println(String.format("%d %8.4g %8.4g %8.4g", 
            i+1, p1, p2, p12));
    }
    System.out.println("sse ");
    for (int i =0; i < steps_.length ; i++) {
        BestStep s = (BestStep) steps_[i];
        double p1 = s.findStep(1).getSse();
        double p2 = s.findStep(2).getSse();
        double p12 =  s.findStep(0).getSse();
        System.out.println(String.format("%d %8.4g %8.4g %8.4g", 
            i+1, p1, p2, p12));
    }
    System.out.println("f   ");
    for (int i =0; i < steps_.length ; i++) {
        BestStep s = (BestStep) steps_[i];
        double p1 = s.findStep(1).getStatistic();
        double p2 = s.findStep(2).getStatistic();
        double p12 =  s.getF12();
        System.out.println(String.format("%d %8.4g %8.4g %8.4g", 
            i+1, p1, p2, p12));
    }
    System.out.println("p   ");
    for (int i =0; i < steps_.length ; i++) {
        BestStep s = (BestStep) steps_[i];
        double p1 = s.findStep(1).getPvalue();
        double p2 = s.findStep(2).getPvalue();
        double p12 =  s.getP12();
        System.out.println(String.format("%d %8.4g %8.4g %8.4g", 
            i+1, p1, p2, p12));
    }
    */
  }

  /*
   * Fit steps to every gene
   *        All the fitted steps are in steps_
   *        The best fit is in bestSteps_
   */
  public void fitStep(AnalysisMetaData meta) throws StepException {
    System.out.println("Fitting Steps .........");
    LinkedList<Step> list = new LinkedList<Step>();
    for (int i =0; i < genes_.length; i++) {
      GeneData gene = data_.getGeneData(genes_[i]);
      Step step = BestStep.fitStep(gene, meta);
      list.add(step);
      step.setGeneIndex(genes_[i]);
      // gene.print();
      // ((BestStep)step).findStep(0).print();
      // ((BestStep)step).findStep(1).print();
      // ((BestStep)step).findStep(2).print();
      // System.exit(0);
    }
    steps_ = new Step[list.size()];
    ListIterator iter = list.listIterator();
    int i = 0;
    while (iter.hasNext()) {
        steps_[i] = (Step) iter.next();
        i++;
    }
    // Copying steps_ to bestSteps_
    bestSteps_ = new Step[steps_.length];
    for (i =0; i < bestSteps_.length ; i++) {
        BestStep s = (BestStep) steps_[i];
        bestSteps_[i] = s.findBestStep(analysisType_);
        bestSteps_[i].setGeneIndex(genes_[i]);
        // System.out.println(sortedSteps_[i].getLabel() + "\t" + genes_[i]);
    }
  }

  /*
   * All genes are centered around the middle of the steps.
   */
  public void performCentering() throws StepException {
    for (int i =0; i < genes_.length; i++) {
      GeneData gene = data_.getGeneData(genes_[i]);
      Step step = bestSteps_[i];
      if (step.getNumSteps() == 0) {
        BestStep s = (BestStep) steps_[i];
        step = s.findStep(1); // get oneStep
      }
      if (step.getNumSteps() >= 0) {
        double center = step.getCenter();
        gene.performCentering(data_.getNumArrayHeader(), center);
        step.performCentering(center);
      }
    }
  }
  
  /*
   * Group genes in to various category
   *     - Sort based on direction and position of steps.
   */
  public void performSorting() throws StepException {
    Integer[] order = new Integer[bestSteps_.length];
    for (int i =0; i < order.length ; i++) {
        order[i] = new Integer(i);
    }
    class GeneComparator implements Comparator<Integer> {
      public GeneComparator() {
        super();
      }
      public int compare(Integer si1, Integer si2) {
        int s1 = si1.intValue();
        int s2 = si2.intValue();
        try {
          int res = bestSteps_[s1].getLabel() -  bestSteps_[s2].getLabel();
          if (res != 0 && bestSteps_[s1].getLabel() == 0) {
            // Label 0 at the end
            res = bestSteps_[s2].getLabel();
          }
          if (res != 0 && bestSteps_[s2].getLabel() == 0) {
            // Label 0 at the end
            res = -bestSteps_[s1].getLabel();
          }
          if (res == 0) { // Same label : Up/Dwn/Up-Down/Down-Up
            int num = bestSteps_[s1].getNumSteps();
            if (num >= 1) { // Match with one or two step
              res = bestSteps_[s1].getStep(0) - bestSteps_[s2].getStep(0);
              if (res == 0 && num > 1) {
                res = bestSteps_[s2].getStep(1) - bestSteps_[s1].getStep(1);
              }
            }
            else { // No Match - collect one step order
              BestStep bs1 = (BestStep) steps_[s1];
              BestStep bs2 = (BestStep) steps_[s2];
              Step step1 = bs1.findStep(1);
              Step step2 = bs2.findStep(1);
              res = step1.getLabel() - step2.getLabel();
              if (res == 0) { // Up/Down
                if (step1.getNumSteps() >= 0 && step2.getNumSteps() >= 0) {
                  res = step1.getStep(0) - step2.getStep(0);
                }
                if (step2.getNumSteps() < 0) {
                  res = -1;
                }
                if (step1.getNumSteps() < 0) {
                  res = 1;
                }
              }
              if (res == 0) { // Same step position
                res = 1;
                if (step1.getPvalue() < step2.getPvalue()) {
                  res = -1;
                }
              }
            }
          }
          if (res == 0) { // group with same step
            res = 1;
            if (bestSteps_[s1].getPvalue() < bestSteps_[s2].getPvalue()) {
              res = -1;
            }
          }
          return res;
        }
        catch(Exception e) {
          e.printStackTrace();
          return 1;
        }
      }
    };
    Arrays.sort(order, new GeneComparator());
    // System.out.println("-----");
    sortedOrder_ = new int[order.length];
    sortedGenes_ = new int[order.length];
    sortedSteps_ = new Step[order.length];
    for (int i=0; i < order.length; i++) {
        sortedOrder_[i] = order[i].intValue();
        sortedGenes_[i] = genes_[order[i].intValue()];
        sortedSteps_[i] = bestSteps_[order[i].intValue()];
        // System.out.println(sortedSteps_[i].getLabel());
    }
  }

  byte[] bytes_; // Temporary variable for sorting

  /*
   *  Compute FDR using a pvalue threshold
   */
  public void performFdrSimple() throws StepException {
    performFiltering();
    AnalysisMetaData meta = new AnalysisMetaData(start_, end_, 
        null, null, pvalueThr_);
    fitStep(meta);
    AnalysisMetaData[] fdrMeta = new AnalysisMetaData[numRandom_];
    FdrInfo[] results = new FdrInfo[numRandom_];
    try {
      Random random = new Random(100);
      for (int i = 0; i < fdrMeta.length; i++) {
        // Get Random permutation.
        bytes_ = new byte[end_ - start_ + 1];
        random.nextBytes(bytes_);
        Integer[] order = new Integer[bytes_.length];
        for (int j =0; j < order.length ; j++) {
          order[j] = new Integer(j);
        }
        class PermComparator implements Comparator<Integer> {
          public int compare(Integer s1, Integer s2) {
            return bytes_[s1.intValue()] - bytes_[s2.intValue()];
          }
        };
        Arrays.sort(order, new PermComparator());
        int[] perm = new int[bytes_.length];
        for (int j =0; j < order.length ; j++) {
            perm[j] = order[j].intValue();
        }

        StepMiner copy = (StepMiner) this.clone();
        fdrMeta[i] = new AnalysisMetaData(start_, end_,
            perm, null, pvalueThr_);
        copy.fitStep(fdrMeta[i]);
        results[i] = copy.getFdrInfo();
        System.out.println(i);
      }
    }
    catch(Exception e) {
      e.printStackTrace();
    }
  }

  /*
   *  Compute FDR using a pvalue threshold
   */
  public void performFdr() throws StepException {
    performFiltering();
    /*
    AnalysisMetaData meta = new AnalysisMetaData(start_, end_, 
        null, null, pvalueThr_);
        */
    int[] st = data_.getStepSearch();
    int[] st1 = data_.getStepSearch1();
    AnalysisMetaData meta = new AnalysisMetaData(start_, end_, 
        null, st, pvalueThr_);
    meta.setStepSearch1(st1);
    fitStep(meta);
    AnalysisMetaData[] fdrMeta = new AnalysisMetaData[numRandom_];
    FdrInfo[] results = new FdrInfo[numRandom_];
    try {
      Random random = new Random(100);
      for (int i = 0; i < fdrMeta.length; i++) {
        // Get Random permutation.
        bytes_ = new byte[end_ - start_ + 1];
        random.nextBytes(bytes_);
        Integer[] order = new Integer[bytes_.length];
        for (int j =0; j < order.length ; j++) {
          order[j] = new Integer(j);
        }
        class PermComparator implements Comparator<Integer> {
          public int compare(Integer s1, Integer s2) {
            return bytes_[s1.intValue()] - bytes_[s2.intValue()];
          }
        };
        Arrays.sort(order, new PermComparator());
        int[] perm = new int[bytes_.length];
        for (int j =0; j < order.length ; j++) {
          perm[j] = order[j].intValue();
        }

        fdrMeta[i] = new AnalysisMetaData(start_, end_,
            perm, null, pvalueThr_);
        results[i] = getFdrInfo();
      }
      fdrMeta[0] = meta;
      for (int i = 0; i < genes_.length; i++) {
        GeneData gene = data_.getGeneData(genes_[i]);
        // results[0] is actual analysis
        // results[i > 1] is random analysis
        for (int j = 0; j < fdrMeta.length; j++) {
          BestStep step = (BestStep) BestStep.fitStep(gene, fdrMeta[j]);
          results[j].p1[i] = step.findStep(1).getPvalue();
          results[j].p2[i] = step.findStep(2).getPvalue();
          results[j].p12[i] = step.getP12();
          results[j].perm = fdrMeta[j].getPermutation();
        }
        if ( i % 500 == 0) {
          System.out.println(i);
        }
      }
      double fdr = FdrInfo.getFdr(analysisType_, pvalueThr_, results);
      fdrStats = FdrInfo.getFdrStats(analysisType_, pvalueThr_, results);
      System.out.println(fdrStats);
    }
    catch(Exception e) {
      e.printStackTrace();
    }
  }

  public void buildGeneIndex(
    HashMap<Integer,SMHashMapUnique<Integer,String> > map,
     Vector<String> allNodes) throws StepException {
    int maxstep = end_ - start_ + 1;
    for (int i = 0; i < 5; i++) {
        SMHashMapUnique<Integer,String> set = 
            new SMHashMapUnique<Integer,String>();
        map.put(new Integer(i), set);
    }
    for (int i = 0; i < genes_.length; i++) {
      String gene = data_.getGenesAt(genes_[i]);
      if (gene != null) {
        allNodes.addElement(gene);
        Step step = bestSteps_[i];
        int label = step.getLabel();
        if (label > 0) {
          SMHashMapUnique<Integer,String> set = (SMHashMapUnique<Integer,String>)
            map.get(new Integer(label));
          int stepIndex = step.getStep(0) - data_.getNumArrayHeader();
          set.put(new Integer(stepIndex), gene);
        } // end label
        else {
          BestStep s = (BestStep) steps_[i];
          step = s.findStep(1); // get oneStep
          if (step.getNumSteps() > 0) {
            SMHashMapUnique<Integer,String> set = (SMHashMapUnique<Integer,String>)
              map.get(new Integer(label));
            int stepIndex = step.getStep(0) - data_.getNumArrayHeader();
            set.put(new Integer(stepIndex), gene);
          }
        }
      } // end gene
    } // end for
  }

  public void performGOAnalysis(String file, Double pvalThr) throws Exception {
    GeneNameScheme ns = data_.getGeneNameScheme();
    SMGOAnalysis goa = new SMGOAnalysis(data_, ns.getOntologyFile(),
        ns.getAnnotationFile(), ns.getOrg(), pvalThr);
    HashMap<Integer,SMHashMapUnique<Integer,String> > map = 
        new HashMap<Integer,SMHashMapUnique<Integer,String> >();
    Vector<String> allNodes = new Vector<String>();
    buildGeneIndex(map, allNodes);
    map.remove(new Integer(0));
    goa.writeHtml(file, allNodes, map);
  }

  public String buildTree(Stack<String> s, int size, Vector<GeneData> res) {
    if (size == 0) {
      return null;
    }
    if (size == 1) {
      String name = (String) s.pop();
      return name;
    }
    String name1 = buildTree(s, size/2, res);
    String name2 = buildTree(s, size - size/2, res);
    String node = "NODE" + nodeNumber_ + "X";
    nodeNumber_++;
    Object[] obj = {node, name1, name2, "1"};
    GeneData data = new GeneData(obj);
    // data.print();
    res.add(data);
    return node;
  }

  public Data getGtrAnnotationData() throws StepException {
    int numArrays = 1;
    int numGenes = sortedGenes_.length;
    int numGeneHeader = 0;
    int numArrayHeader = 3;
    GeneData[] data = new GeneData[numGenes+numGeneHeader];
    Stack<String> stack = new Stack<String>();
    Stack<String> tstack = new Stack<String>();
    Stack<String> ttstack = new Stack<String>();
    int lastlabel = 0;
    int lastdir = 1;
    int index = 0;
    nodeNumber_ = 0;
    HashMap<String,DAGNode> mapping = new HashMap<String,DAGNode>();
    for (int j =0; j < sortedGenes_.length; j++) {
      int dir = sortedSteps_[j].getLabel();
      int step1 = 0;
      int step2 = 0;
      if (sortedSteps_[j].getNumSteps() > 0) {
        step1 = sortedSteps_[j].getStep(0);
      }
      if (sortedSteps_[j].getNumSteps() > 1) {
        step2 = sortedSteps_[j].getStep(1);
      }
      if (dir == 0) {
        BestStep s = (BestStep) steps_[sortedOrder_[j]];
        Step step = s.findStep(1); // get oneStep
        if (step.getNumSteps() > 0) {
          step1 = step.getStep(0);
        }
      }
      int label = dir+ step1+ step2;
      if (stack.size() > 0 && label != lastlabel) {
        Vector<GeneData> obj = new Vector<GeneData>();
        String sname = buildTree(stack, stack.size(), obj);
        tstack.push(sname);
        Enumeration<GeneData> e = obj.elements();
        while (e.hasMoreElements()) {
            GeneData gdata = (GeneData) e.nextElement();
            data[index++] = gdata;
        }
      }
      if (tstack.size() > 0 && dir != lastdir) {
        Vector<GeneData> obj = new Vector<GeneData>();
        String sname = buildTree(tstack, tstack.size(), obj);
        ttstack.push(sname);
        Enumeration<GeneData> e = obj.elements();
        while (e.hasMoreElements()) {
            GeneData gdata = (GeneData) e.nextElement();
            data[index++] = gdata;
        }
      }
      String name = "GENE"+ sortedGenes_[j]+"X";
      stack.push(name);
      lastlabel = label;
      lastdir = dir;
    }
    if (stack.size() > 0) {
      Vector<GeneData> obj = new Vector<GeneData>();
      String sname = buildTree(stack, stack.size(), obj);
      tstack.push(sname);
      Enumeration<GeneData> e = obj.elements();
      while (e.hasMoreElements()) {
        GeneData gdata = (GeneData) e.nextElement();
        data[index++] = gdata;
      }
    }
    if (tstack.size() > 0) {
      Vector<GeneData> obj = new Vector<GeneData>();
      String sname = buildTree(tstack, tstack.size(), obj);
      ttstack.push(sname);
      Enumeration<GeneData> e = obj.elements();
      while (e.hasMoreElements()) {
        GeneData gdata = (GeneData) e.nextElement();
        data[index++] = gdata;
      }
    }
    Vector<GeneData> obj = new Vector<GeneData>();
    String sname = buildTree(ttstack, ttstack.size(), obj);
    Enumeration<GeneData> e = obj.elements();
    while (e.hasMoreElements()) {
      GeneData gdata = (GeneData) e.nextElement();
      data[index++] = gdata;
    }
    numGenes = index;
    // creating nodes
    for (int i = 0; i  <numGenes; i++) {
        String n1 = (String) data[i].getDataAt(0);
        String n2 = (String) data[i].getDataAt(1);
        String n3 = (String) data[i].getDataAt(2);
        if (!mapping.containsKey(n1)) {
            mapping.put(n1, new DAGNode(n1));
        }
        if (!mapping.containsKey(n2)) {
            mapping.put(n2, new DAGNode(n2));
        }
        if (!mapping.containsKey(n3)) {
            mapping.put(n3, new DAGNode(n3));
        }
    }
    // creating edges
    for (int i = 0; i  <numGenes; i++) {
        String s1 = (String) data[i].getDataAt(0);
        String s2 = (String) data[i].getDataAt(1);
        String s3 = (String) data[i].getDataAt(2);
        DAGNode n1 = mapping.get(s1);
        DAGNode n2 = mapping.get(s2);
        DAGNode n3 = mapping.get(s3);
        n1.addChild(n2);
        n1.addChild(n3);
    }
    // creating find roots.
    Iterator<String> itr = mapping.keySet().iterator();
    DAGGraph g = new DAGGraph();
    while (itr.hasNext()) {
        String key = (String) itr.next();
        DAGNode n1 = mapping.get(key);
        if (n1.getNumParents() == 0) {
            g.addRoot(n1);
        }
    }
    int depth = g.getDepth();
    System.out.println("Depth of the tree = " + depth);
    for (int i = 0; i  <numGenes; i++) {
        String s1 = (String) data[i].getDataAt(0);
        DAGNode n1 = mapping.get(s1);
        Integer ndepth = (Integer)n1.getAttribute("depth0");
        double score = (depth - ndepth.intValue() + 2) * 1.0 / (depth + 2);
        data[i].setDataAt(3, "" + score);
    }
    Data res = new Data(numArrays, numGenes, numGeneHeader, numArrayHeader, data);
    return res;
  }

  public Data getCdtAnnotationData() throws StepException {
    int numArrays = data_.getNumArrays();
    int numGenes = sortedGenes_.length;
    int numGeneHeader = 1;
    int numArrayHeader = data_.getNumArrayHeader() + 1;
    GeneData[] data = new GeneData[numGenes+numGeneHeader];
    Object[] h1Add = {"GID"};
    try {
    // Creating headers
    if ( data_.getNumGeneHeader() > 0) {
      data[0] = data_.getGeneData(0).insert(0, h1Add);
    }
    else {
      Object[] h1 = new Object[numArrayHeader+numArrays];
      for (int i =0; i < numArrayHeader; i++) {
        h1[i] = h1Add[i];
      }
      for (int i =numArrayHeader; i < h1.length; i++) {
        h1[i] = "" + (i - numArrayHeader);
      }
      data[0] = new GeneData(h1);
    }
/*
    if (data[1] == null) {
      Object[] h2 = new Object[numArrayHeader+numArrays];
      for (int i =0; i < numArrayHeader; i++) {
        h2[i] = " ";
      }
      for (int i =numArrayHeader; i < h2.length; i++) {
        h2[i] = "ARRY" + (i-numArrayHeader) + "X";
      }
      h2[0] = "AID";
      data[1] = new GeneData(h2);
    }
*/
    for (int j =0; j < sortedGenes_.length; j++) {
      GeneData gene = data_.getGeneData(sortedGenes_[j]);
      Object[] obj = { "GENE"+ sortedGenes_[j]+"X" };
      data[j+1] = gene.insert(0, obj);
    }
    }
    catch(ArrayException e) {
        throw new StepException(e.getMessage());
    }
    Data res = new Data(numArrays, numGenes, numGeneHeader, numArrayHeader, data);
    return res;
  }

  public Data getStepAnnotationData() throws StepException {
    int numArrays = data_.getNumArrays();
    int numGenes = sortedGenes_.length;
    int numGeneHeader = 2;
    int numArrayHeader = data_.getNumArrayHeader() + 5;
    GeneData[] data = new GeneData[numGenes+numGeneHeader];
    Object[] h1Add = {"label", "dir", "step1", "step2", "pvalue"};
    Object[] h2Add = {"1", "1", "1", "1", "1"};
    
    try {
    // Creating headers
    if ( data_.getNumGeneHeader() > 0) {
      data[0] = data_.getGeneData(0).insert(data_.getNumArrayHeader(), h1Add);
      if (data_.getNumGeneHeader() > 1) {
        data[1] = data_.getGeneData(1).insert(data_.getNumArrayHeader(), h2Add);
      }
    }
    else {
      Object[] h1 = new Object[numArrayHeader+numArrays];
      for (int i =0; i < numArrayHeader; i++) {
        h1[i] = h1Add[i];
      }
      for (int i =numArrayHeader; i < h1.length; i++) {
        h1[i] = "" + (i - numArrayHeader);
      }
      data[0] = new GeneData(h1);
    }
    if (data[1] == null) {
      Object[] h2 = new Object[numArrayHeader+numArrays];
      for (int i =0; i < h2.length; i++) {
        h2[i] = "1";
      }
      h2[0] = "EWEIGHT";
      data[1] = new GeneData(h2);
    }
    for (int j =0; j < sortedGenes_.length; j++) {
      GeneData gene = data_.getGeneData(sortedGenes_[j]);
      int label = 0;
      int dir = 0;
      int step1 = 0;
      int step2 = 0;
      dir = sortedSteps_[j].getLabel()-1;
      if (dir < 0) {
        dir = 0;
      }
      if (dir == 1 || dir == 0) {
        label = 2;
      }
      if (dir == 2 || dir == 3) {
        label = 3;
      }
      if (sortedSteps_[j].getNumSteps() > 0) {
        step1 = sortedSteps_[j].getStep(0) - 2;
      }
      if (sortedSteps_[j].getNumSteps() > 1) {
        step1--;
        step2 = sortedSteps_[j].getStep(1) - step1 - 2;
      }
      double pvalue = sortedSteps_[j].getPvalue();
      Object[] obj = { ""+label, ""+dir, ""+step1, ""+step2 , ""+pvalue};
      data[j+2] = gene.insert(data_.getNumArrayHeader(), obj);
    }
    }
    catch(ArrayException e) {
        throw new StepException(e.getMessage());
    }
    Data res = new Data(numArrays, numGenes, numGeneHeader, numArrayHeader, data);
    return res;
  }

  public Data getStepOrderedData() throws StepException {
    int[] order = new int[sortedGenes_.length];
    int numHeader = data_.getNumGeneHeader();
    for (int i =0; i < order.length; i++) {
      order[i] = sortedGenes_[i] - numHeader;
    }
    Data res = Data.selectGenesFromData(data_, order);
    return res;
  }

  /*
   *  Example : getStepOrderedData("Up")
   *            getStepOrderedData("Down")
   */
  public Data getStepOrderedData(String filetag) throws StepException {
    HashSet<Integer> selectedLabels = Step.getSelectedLabels(filetag);
    int count = 0;
    for (int i =0; i < sortedGenes_.length; i++) {
      Step step = sortedSteps_[i];
      if (!selectedLabels.contains(new Integer(step.getLabel()))) {
        continue;
      }
      count ++;
    }
    int[] order = new int[count];
    int numHeader = data_.getNumGeneHeader();
    int index = 0;
    for (int i =0; i < sortedGenes_.length; i++) {
      Step step = sortedSteps_[i];
      if (!selectedLabels.contains(new Integer(step.getLabel()))) {
        continue;
      }
      order[index++] = sortedGenes_[i] - numHeader;
    }
    Data res = Data.selectGenesFromData(data_, order);
    return res;
  }

  /*
   * Write annotation files.
   */
  public void writeAnnotations(String file) throws Exception {
    String filename = Data.getFileName(file);
    String filetag = Data.getFileTag(file);
    HashSet<Integer> selectedLabels = Step.getSelectedLabels(filetag);
    BufferedWriter out = new BufferedWriter(new FileWriter(filename));
    System.out.println("Writing " + file + " ...");
    out.write("Name\t" + Step.headString() + "\n");
    for (int i =0; i < sortedGenes_.length; i++) {
      Step step = sortedSteps_[i];
      if (!selectedLabels.contains(new Integer(step.getLabel()))) {
        continue;
      }
      String geneName = data_.getGenesAt(sortedGenes_[i]);
      out.write(geneName + "\t" + step.tabString() + "\n");
    }
    out.close();
    System.out.println("Done");
  }

  /* 
   * Plotting Steps
   */
  public void plotSteps(String file) throws Exception {
    String filename = Data.getFileName(file);
    String filetag = Data.getFileTag(file);
    HashSet<Integer> selectedLabels = Step.getSelectedLabels(filetag);
    PSPlot plotter = new PSPlot(filename);
    plotter.open();
    plotter.array(5 /* Rows */, 2 /* columns */);
    for (int i =0; i < sortedGenes_.length; i++) {
      Step step = sortedSteps_[i];
      if (!selectedLabels.contains(new Integer(step.getLabel()))) {
        continue;
      }
      GeneData gene = data_.getGeneData(sortedGenes_[i]);
      String geneName = data_.getGenesAt(sortedGenes_[i]);
      gene.convertDouble(start_, end_);
      Object[] geneData = gene.getData();
      Vector<Double> data = new Vector<Double>();
      Vector<Double> time = new Vector<Double>();
      for (int j = start_; j <= end_; j++) {
        data.add((Double)geneData[j]);
        time.add( new Double(j-start_));
      }
      plotter.plot(time, data);
      if (step.getNumSteps() > 0) {
        Double x1 = time.get(0);
        int step1 = step.getStep(0);
        Double x2 = new Double((time.get(step1 - start_).doubleValue()
            + time.get(step1 - start_ + 1).doubleValue() )/2);
        Double m1 = new Double(step.getMean(0));
        Double m2 = new Double(step.getMean(1));
        plotter.line(x1, m1, x2, m1);
        plotter.line(x2, m1, x2, m2);
        if (step.getNumSteps() > 1) {
          int step2 = step.getStep(1);
          Double x3;
          if ((step2 - start_ + 1) < time.size()) {
            x3 = new Double((time.get(step2 - start_).doubleValue()
                  + time.get(step2 - start_ + 1).doubleValue() )/2);
          }
          else {
            x3 = time.get(step2 - start_);
          }
          plotter.line(x2, m2, x3, m2);
          plotter.line(x3, m2, x3, m1);
          Double x4 = time.get(end_ - start_);
          plotter.line(x3, m1, x4, m1);
        }
        else {
          Double x3 = time.get(end_ - start_);
          plotter.line(x2, m2, x3, m2);
        }
      }
      String label = "None";
      if (step.getNumSteps() == 1) {
        label = "OneStep";
      }
      if (step.getNumSteps() == 2) {
        label = "TwoStep";
      }
      plotter.xlabel(label + " , p = " + step.getPvalueStr());
      plotter.ylabel(" Gene expression ");
      plotter.title(geneName);
    }
    plotter.close();
  }

  public void buildGeneSets(SMHashMapUnique<String,String> map)
    throws StepException {
    HashSet<String> allSets = new HashSet<String>();
    for (int i = 0; i < genes_.length; i++) {
      String gene = data_.getGenesAt(genes_[i]);
      if (gene != null) {
        Step step = bestSteps_[i];
        HashSet<String> tags = step.getTags(data_.getNumArrayHeader());
        allSets.addAll(tags);
        Iterator<String> itr = tags.iterator();
        while (itr.hasNext()) {
          String tag = (String) itr.next();
          map.put(tag, gene);
        }
      } // end gene
    } // end for

  }

  public GeneSet createGeneSets() throws StepException {
    GeneSet res = new GeneSet();
    for (int i = 0; i < genes_.length; i++) {
      String gene = data_.getGenesAt(genes_[i]);
      if (gene != null) {
        Step step = bestSteps_[i];
        HashSet<String> tags = step.getTags(data_.getNumArrayHeader());
        Iterator<String> itr = tags.iterator();
        while (itr.hasNext()) {
          String tag = (String) itr.next();
          res.add(tag, gene);
        }
      } // end gene
    } // end for

    Vector<String> allSets = res.getAllSets();
    // Order Gene Sets
    Integer[] order = Step.orderTags(allSets);
    res.reorder(order);
    return res;
  }

  public void writeGeneSets(String file) throws StepException {
    try {
       GeneSet set = createGeneSets();
       if (deleteGenesInMultipleGroups_) {
         set = deleteGenesInMultipleGroups();
       }
       set = GeneSet.toUpperCase(set);
       GeneSetWriter.writeFile(set, file);
    }
    catch(Exception e) {
        e.printStackTrace();
        throw new StepException(e.getMessage());
    }
  }

  void buildDeletedGenes_(HashSet<String> deletedGenes, 
      HashSet<String> primary, HashSet<String> other1, HashSet<String> other2,
      HashSet<String> other3, HashSet<String> other4) throws StepException {
    if (primary == null) {
        return;
    }
    Iterator<String> itr = primary.iterator();
    while (itr.hasNext()) {
      String name = (String) itr.next();
      // System.out.println(name);
      if (other1 != null && other1.contains(name)) {
        deletedGenes.add(name);
      }
      if (other2 != null && other2.contains(name)) {
        deletedGenes.add(name);
      }
      if (other3 != null && other3.contains(name)) {
        deletedGenes.add(name);
      }
      if (other4 != null && other4.contains(name)) {
        deletedGenes.add(name);
      }
    }
  }

  public GeneSet deleteGenesInMultipleGroups() throws StepException {
    GeneSet set = createGeneSets();
    HashSet<String> up = set.getSet("Up");
    HashSet<String> down = set.getSet("Down");
    HashSet<String> updown = set.getSet("UpDown");
    HashSet<String> downup = set.getSet("DownUp");
    HashSet<String> rest = set.getSet("Rest");
    HashSet<String> deletedGenes = new HashSet<String>();
    buildDeletedGenes_(deletedGenes, up, down, updown, downup, rest);
    buildDeletedGenes_(deletedGenes, down, updown, downup, rest, up);
    buildDeletedGenes_(deletedGenes, updown, downup, rest, up, down);
    buildDeletedGenes_(deletedGenes, downup, rest, up, down, updown);
    buildDeletedGenes_(deletedGenes, rest, up, down, updown, downup);
    System.out.println("Number of genes deleted = " + deletedGenes.size());
    System.out.println("Deleted genes:");
    Iterator<String> itr = deletedGenes.iterator();
    int index = 0;
    while (itr.hasNext()) {
      String name = (String) itr.next();
      index++;
      System.out.print(name + "\t");
      if ( (index % 10) == 0) {
        System.out.println();
      }
    }
    System.out.println();
    int count = 0;
    for (int i =0; i < sortedGenes_.length; i++) {
      String geneName = data_.getGenesAt(sortedGenes_[i]);
      if (!deletedGenes.contains(geneName)) {
        count++;
      }
    }
    
    int[] sortedOrder_n = new int[count];
    int[] sortedGenes_n = new int[count];
    Step[] sortedSteps_n = new Step[count];
    index = 0;
    for (int i =0; i < sortedGenes_.length; i++) {
      String geneName = data_.getGenesAt(sortedGenes_[i]);
      if (!deletedGenes.contains(geneName)) {
        sortedOrder_n[index] = sortedOrder_[i];
        sortedGenes_n[index] = sortedGenes_[i];
        sortedSteps_n[index] = sortedSteps_[i];
        index++;
      }
    }
    sortedOrder_ = sortedOrder_n;
    sortedGenes_ = sortedGenes_n;
    sortedSteps_ = sortedSteps_n;
    set = GeneSet.deleteGenes(set, deletedGenes);
    return set;
  }

  public void performGOBestAnalysis(String file, Double pvalThr) throws Exception {
    GeneNameScheme ns = data_.getGeneNameScheme();
    SMGOAnalysis goa = new SMGOAnalysis(data_, ns.getOntologyFile(),
        ns.getAnnotationFile(), ns.getOrg(), pvalThr);
    SMHashMapUnique<String,String> map = new SMHashMapUnique<String,String>();
    buildGeneSets(map);
    goa.writeHtml(file, map);
  }

  /*
   * Get only selected Data (e.d Step: or Up:)
   */
  public Data getTaggedData(String filetag) {
    HashSet<Integer> selectedLabels = Step.getSelectedLabels(filetag);
    Vector<Integer> order = new Vector<Integer>();
    for (int i =0; i < sortedGenes_.length; i++) {
      Step step = sortedSteps_[i];
      if (!selectedLabels.contains(new Integer(step.getLabel()))) {
        continue;
      }
      int index = sortedGenes_[i] - data_.getNumGeneHeader();
      order.add(new Integer(index));
    }
    int[] ord = new int[order.size()];
    for (int i=0; i < order.size(); i++) {
        Integer val = order.get(i);
        ord[i] = val.intValue();
    }
    Data res = Data.selectGenesFromData(data_, ord);
    return res;
  }

  /*
   * Write PCL files.
   */
  public void writePCL(String file) throws Exception {
    String filename = Data.getFileName(file);
    String filetag = Data.getFileTag(file);
    Data res = getTaggedData(filetag);
    /*
    // For writing the entire data
    int[] ord = new int[res.getNumGenes()];
    for (int i=0; i < ord.length; i++) {
        ord[i] = i;
    }
    PCLFileWriter.writeFile(res, filename, ord);
    */
    PCLFileWriter.writeFile(res, filename);
  }

  /*
   * Data for Ingenuity analysis
   */

  public Data getStepAnnotationByGenes() throws StepException {
    Vector<GeneData> v = new Vector<GeneData>();
    for (int i =0; i < sortedGenes_.length; i++) {
      Step step = sortedSteps_[i];
      GeneData gene = data_.getGeneData(sortedGenes_[i]);
      String geneName = data_.getGenesAt(sortedGenes_[i]);
      if (geneName == null) {
        geneName = (String) gene.getDataAt(0);
      }
      int label = step.getLabel(); 
      // label : 0 - rest, 1 - Up, 2 -Down, 3 - UpDown, 4 - DownUp
      double value = 0;
      try {
        //value = gene.foldChangeLog(start_, end_);
        if (label > 0) {
          value = step.getStep(0);
        }
      }
      catch(Exception e) { }
      if (label == 2 || label == 4) {
        value = -value;
      }
      if (label == 0 || label > 2) {
        // value = 0;
      }
      Object[] d = new Object[2];
      d[0] = geneName;
      d[1] = new Double(value);
      v.add(new GeneData(d));
    }
    Data res = new Data(2, v.size(), 0, 0, v);
    return res;
  }

  /*
   *  Perform StepMiner analysis
   */
  public void performAnalysis() throws StepException {
    init_data_();
    performFiltering();
    if (stopAnalysis_) {
        return;
    }
    int[] st = data_.getStepSearch();
    int[] st1 = data_.getStepSearch1();
    AnalysisMetaData meta = new AnalysisMetaData(start_, end_, 
        null, st, pvalueThr_);
    meta.setStepSearch1(st1);
    meta.print();
    fitStep(meta);
    performSorting();
    if (stepCentering_) {
      performCentering();
    }
    if (deleteGenesInMultipleGroups_) {
      deleteGenesInMultipleGroups();
    }
    printStats();
    if (fdrAnalysis_) {
      performFdr();
    }
  }

};

class FdrInfo {
    public double[] p1;
    public double[] p2;
    public double[] p12;
    public int[] perm;

    public FdrInfo(double[] a1, double[] a2, double[] a12) {
        p1 = a1;
        p2 = a2;
        p12 = a12;
    }

    public int numSignificantSingle(double thr) {
        int count = 0;
        if (p1 != null) {
            for (int i =0; i < p1.length; i++) {
                if (p1[i] < thr) {
                    count++;
                }
            }
        }
        return count;
    }

    public int numSignificantBoth(double thr) {
        int count = 0;
        if (p1 != null && p2 != null && p12 != null) {
            for (int i =0; i < p1.length; i++) {
                if (p1[i] < thr && p12[i] > thr) {
                    count++;
                }
                else if (p2[i] < thr) {
                    count++;
                }
            }
        }
        return count;
    }

    public int numSignificantSelectTwo(double thr) {
        int count = 0;
        if (p1 != null && p2 != null && p12 != null) {
            for (int i =0; i < p1.length; i++) {
                if (p1[i] < thr && p12[i] > thr) {
                    // count++;
                    // Don't count
                }
                else if (p2[i] < thr) {
                    count++;
                }
            }
        }
        return count;
    }

    public int numSignificantTwo(double thr) {
        int count = 0;
        if (p2 != null) {
            for (int i =0; i < p1.length; i++) {
                if (p2[i] < thr) {
                    count++;
                }
            }
        }
        return count;
    }

    public int numSignificant(int type, double thr) {
        if (type == 0) {
            return numSignificantSingle(thr);
        }
        else if (type == 1) {
            return numSignificantBoth(thr);
        }
        else if (type == 2) {
            return numSignificantTwo(thr);
        }
        else if (type == 3) {
            return numSignificantSelectTwo(thr);
        }
        System.err.println("Undefined analysis type : " + type);
        return 0;
    }

    public static double getFdr(int type, double thr, FdrInfo[] results) {
        int sum = 0;
        int max = 0;
        int maxindex = 0;
        for (int i =0; i < results.length; i++) {
            int count = results[i].numSignificant(type, thr);
            sum += count;
            if (count > max) {
                max = count;
                maxindex = i;
            }
        }
        double mean = sum/ (double) results.length;
        int actual = results[0].numSignificant(type, thr);
        double fdr = mean/actual;
        return fdr;
    }

    public static String getFdrStats(int type, double thr, FdrInfo[] results) {
        int sum = 0;
        int max = 0;
        int maxindex = 0;
        for (int i =0; i < results.length; i++) {
            int count = results[i].numSignificant(type, thr);
            sum += count;
            if (count > max) {
                max = count;
                maxindex = i;
            }
        }
        double mean = sum/ (double) results.length;
        int actual = results[0].numSignificant(type, thr);
        double fdr = mean/actual;
        String res = (" Fdr = " + fdr + " mean = " + mean + " actual = " +
            actual) + "\n" +
        ("Permutation that gives the max significant genes ("+
            max + "):") + "\n";
        int[] perm = results[maxindex].perm;
        if (perm == null) {
          res += ("Original order");
        }
        else {
          for (int i =0; i < perm.length; i++) {
              res += (perm[i] + " ");
          }
        }
        res += "\n";
        
        return res;
    }

    public void print(int i) {
        System.out.println("p1 : " + p1[i] + " p2 : " + 
            p2[i] + " p12 : " + p12[i]);
    }
    public void print() {
      for (int i =0; i < p1.length; i++) {
        print(i);
      }
    }
}

