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

public class BooleanSubAnalysis {

  String output_filename_;
  String bv_filenamex_;
  String bv_filename_;
  String ph_filename_;
  String phid_;

  double threshold_;
  double stat_threshold_;
  double single_threshold_;
  int single_cutoff_;
  int blocksize_;

  PCLFileReader reader_;
  PCLFileReader readerx_;
  NetworkFile out_;
  int fileType_;
  BitSet phenotype_;
  int numArrays_;
  int numPhArrays_;
  int numRows_x_;

  Vector<String> high_;
  Vector<String> low_;
  Vector<String> balancedx_;
  Vector<String> balancedy_;
  HashMap<Integer, Integer> balanced_map_;

  String geneList_;
  HashMap<String, Integer> geneMaps_;
  NetworkInfo info_;

  public BooleanSubAnalysis(String bvfilex, String bvfile, String ofile) {
    output_filename_ = ofile;
    bv_filenamex_ = bvfilex;
    bv_filename_ = bvfile;
    ph_filename_ = null;
    phid_ = null;
    init();
  }

  public BooleanSubAnalysis(String bvfilex, String bvfile, String ofile, String phfile, String id) {
    output_filename_ = ofile;
    bv_filenamex_ = bvfilex;
    bv_filename_ = bvfile;
    ph_filename_ = phfile;
    phid_ = id;
    init();
  }

  public void init() {
    threshold_ = 0.1;
    stat_threshold_ = 3;
    single_threshold_ = 0.05;
    single_cutoff_ = 20;
    blocksize_ = 8000;
    fileType_ = NetworkFile.FILE_1_2;
  }

  public void setThreshold(double t) { threshold_ = t; }
  public void setStatThreshold(double t) { stat_threshold_ = t; }
  public void setSingleThreshold(double t) { single_threshold_ = t; }
  public void setSingleCutoff(int t) { single_cutoff_ = t; }
  public void setBlockSize(int s) { blocksize_ = s; }
  public void setGeneList(String s) { geneList_ = s; }
  public void setFileType(int t) { fileType_ = t; }

  public void performSingleListAnalysis() throws IOException {
    beginAnalysis();
    readGeneList();
    performListSingleAnalysis();
    out_.startMatrix(0, 3);
    finish();
  }
  public void performListAnalysis() throws IOException {
    beginAnalysis();
    readGeneList();
    performListSingleAnalysis();
    performListPairsAnalysis();
    finish();
  }

  public void readGeneList() throws IOException {
    PCLFileReader listr = new PCLFileReader(geneList_);
    listr.begin();
    info_.readBvIndex();
    HashSet<String> geneIds = new HashSet<String>();
    geneMaps_ = new HashMap<String, Integer>();
    GeneData g = listr.getData();
    for (int i = 0;  g != null; i++) {
      String id = (String) g.getDataAt(0);
      id = id.trim();
      geneIds.add(id);
      g = listr.getData();
      Long index = info_.getIndexByName(id);
      if (index != null) {
        geneMaps_.put(id, new Integer(index.intValue()));
      }
    }
    reader_.beginLightRandomAccess();
    reader_.setLineMap(info_.getLineMap());
    System.out.println("Found " + geneMaps_.size() + " probes");
    g = readerx_.getDataAt(0);
    for (int i = 0;  g != null; i++) {
      String id = (String) g.getDataAt(0);
      id = id.trim();
      geneMaps_.put(id, new Integer(i));
      g = readerx_.getDataAt(i+1);
    }
  }

  public void performListSingleAnalysis() throws IOException {
    high_ = new Vector<String>();
    low_ = new Vector<String>();
    balancedx_ = new Vector<String>();
    balancedy_ = new Vector<String>();
    balanced_map_ = new HashMap<Integer, Integer>();
    int index = 0;
    GeneData g = readerx_.getDataAt(0);
    for (int i = 0;  g != null; i++) {
      String id = (String) g.getDataAt(0);
      id = id.trim();
      BitSet va = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 0);
      BitSet va_thr = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 1);
      va = setPhenotype(va);
      va_thr = setPhenotype(va_thr);
      if (!haveGoodDynamicRange(va_thr)) {
        g = readerx_.getDataAt(i+1);
        continue;
      }
      int group = findBias(va, va_thr);
      if (group == 0) {
        low_.add(id);
      }
      if (group == 1) {
        high_.add(id);
      }
      if (group == -1) {
        balancedy_.add(id);
        balancedx_.add(id);
        balanced_map_.put(new Integer(i), new Integer(index));
        index++;
      }
      g = readerx_.getDataAt(i+1);
    }
    numRows_x_ = (int)readerx_.getLineNumber();
    Iterator<String> keys = geneMaps_.keySet().iterator();
    while (keys.hasNext()) {
      String id = (String) keys.next();
      Integer loc = geneMaps_.get(id);
      g = reader_.getDataAt(loc.intValue());
      BitSet va = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 0);
      BitSet va_thr = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 1);
      va = setPhenotype(va);
      va_thr = setPhenotype(va_thr);
      if (!haveGoodDynamicRange(va_thr)) {
        continue;
      }
      int group = findBias(va, va_thr);
      if (group == 0) {
        low_.add(id);
      }
      if (group == 1) {
        high_.add(id);
      }
      if (group == -1) {
        loc = new Integer(loc.intValue() + numRows_x_);
        balancedy_.add(id);
        balanced_map_.put(loc, new Integer(index));
        index++;
      }
    }
    out_.writeList("low", low_);
    out_.writeList("high", high_);
    out_.writeList("balancedx", balancedx_);
    out_.writeList("balancedy", balancedy_);
  }

  public void performListPairsAnalysis() throws IOException {
    out_.startMatrix(balancedy_.size(), 3);
    Enumeration<String> alist = balancedx_.elements();
    while (alist.hasMoreElements()) {
      String aid = alist.nextElement();
      Integer aloc = geneMaps_.get(aid);
      GeneData ga = readerx_.getDataAt(aloc.intValue());
      BitSet va = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 0);
      BitSet va_thr = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 1);
      va = setPhenotype(va);
      va_thr = setPhenotype(va_thr);
      if (!haveGoodDynamicRange(va_thr)) {
        continue;
      }
      int index = -1;
      Enumeration<String> blist = balancedy_.elements();
      while (blist.hasMoreElements()) {
        index++;
        String bid = blist.nextElement();
        Integer bloc = geneMaps_.get(bid);
        int bloc_offset = bloc.intValue();
        GeneData gb;
        if (index < balancedx_.size()) {
           gb = readerx_.getDataAt(bloc.intValue());
        }
        else {
          bloc_offset = bloc.intValue() + numRows_x_;
          gb = reader_.getDataAt(bloc.intValue());
        }
        BitSet vb = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 0);
        BitSet vb_thr = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 1);
        vb = setPhenotype(vb);
        vb_thr = setPhenotype(vb_thr);
        if (!haveGoodDynamicRange(vb_thr)) {
          continue;
        }
        performSinglePairAnalysis(aloc.intValue(), bloc_offset, va, va_thr, vb, vb_thr);
      }
    }
  }

  public void performAnalysis() throws IOException {
    System.out.println("Initializing...");
    beginAnalysis();
    System.out.println("Done");
    performSingleAnalysis();
    System.out.println("Done Single");
    performBlockAnalysis();
    System.out.println("Done Block");
    finish();
  }

  public void finish() throws IOException {
    out_.close();
  }

  public void beginAnalysis() throws IOException {
    info_ = new NetworkInfo();
    String networkFile = new String(bv_filename_);
    networkFile = networkFile.trim();
    networkFile = networkFile.replaceFirst(".bv$", "");
    info_.setAllFiles(networkFile);
    reader_ = new PCLFileReader(bv_filename_);
    readerx_ = new PCLFileReader(bv_filenamex_);
    PCLFileReader.CACHE_SIZE = 0;
    reader_.begin();
    readerx_.beginRandomAccess();
    out_ = new BitMatrixSubNetworkFile(output_filename_);
    out_.writeHeader();
    GeneData g = readerx_.getDataAt(1);
    String eweight = (String) g.getDataAt(2);
    BitSet va = BitSetUtils.stringToBitSet(eweight, 0);
    numArrays_ = eweight.length();
    phenotype_ = new BitSet(numArrays_);
    for (int i =0; i < numArrays_; i++) {
      phenotype_.set(i);
    }
    if (phid_ != null && !phid_.equals("All")) {
      PCLFileReader phr = new PCLFileReader(ph_filename_);
      phr.begin();
      if (phr.getNumColumns() != (numArrays_+3)) {
        System.out.println("BitVector file has different number of arrays than the phenotype file");
        System.out.println(phr.getNumColumns() + " != " + (numArrays_+3));
        System.exit(1);
      }
      GeneData header = phr.getHeader();
      int found = 0;
      if (header != null) {
        String id = (String) header.getDataAt(0);
        if (id.equals(phid_)) {
          int start = 3;
          int end = numArrays_+2;
          System.out.println("Found Phenotype : " + phid_);
          found = 1;
          Object[] data = header.getData();
          for (int i = start; i <= end; i++) {
            String str = (String) data[i];
            if (str.equals("0")) {
              phenotype_.clear(i-start);
            }
          }
        }
      }
      while (found == 0 && phr.hasNext()) {
        int start = 3;
        int end = numArrays_+2;
        GeneData ph = phr.getData();
        if (ph == null) {
          break;
        }
        String id = (String) ph.getDataAt(0);
        if (id.equals(phid_)) {
          System.out.println("Found Phenotype : " + phid_);
          found = 1;
          Object[] data = ph.getData();
          for (int i = start; i <= end; i++) {
            String str = (String) data[i];
            if (str.equals("0")) {
              phenotype_.clear(i-start);
            }
          }
          break;
        }
      }
    }
    numPhArrays_ = phenotype_.cardinality();
  }

  public int findBias(BitSet v, BitSet thr) {
    int c1 = v.cardinality();
    BitSet tmp = (BitSet) thr.clone();
    tmp.andNot(v);
    int c0 = tmp.cardinality();
    int total = c0 + c1;
    if (total <= 0) {
      return -1;
    }
    double p = (c0/(c0+c1+0.0));
    //System.out.println(total + "\t" + c0 + "\t" + c1 + "\t" + p);
    if (c0 < single_cutoff_ && p < single_threshold_ ) {
      return 1;
    }
    if (c1 < single_cutoff_ && (1-p) < single_threshold_ ) {
      return 0;
    }
    return -1;
  }

  public boolean haveGoodDynamicRange(BitSet va_thr) {
    int num = numPhArrays_;
    int outside = va_thr.cardinality();
    if (num > (3 * outside)) {
      return false;
    }
    else {
      return true;
    }
  }

  public BitSet setPhenotype(BitSet v) {
    // v.and(phenotype_);
    BitSet res = new BitSet(numPhArrays_);
    int index = 0;
    for (int i =0; i < numArrays_; i++) {
      if (phenotype_.get(i)) {
        res.clear(index);
        if (v.get(i)) {
          res.set(index);
        }
        index++;
      }
    }
    return res;
  } 

  public void performSingleAnalysis() throws IOException {
    high_ = new Vector<String>();
    low_ = new Vector<String>();
    balancedx_ = new Vector<String>();
    balancedy_ = new Vector<String>();
    balanced_map_ = new HashMap<Integer, Integer>();
    int index = 0;
    GeneData g = readerx_.getDataAt(0);
    for (int i = 0;  g != null; i++) {
      String id = (String) g.getDataAt(0);
      BitSet va = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 0);
      BitSet va_thr = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 1);
      va = setPhenotype(va);
      va_thr = setPhenotype(va_thr);
      if (!haveGoodDynamicRange(va_thr)) {
        //System.out.println(i + "\tbd");
        g = readerx_.getDataAt(i+1);
        continue;
      }
      int group = findBias(va, va_thr);
      //System.out.println(i + "\t" + group);
      if (group == 0) {
        low_.add(id);
      }
      if (group == 1) {
        high_.add(id);
      }
      if (group == -1) {
        balancedx_.add(id);
        balancedy_.add(id);
        balanced_map_.put(new Integer(i), new Integer(index));
        index++;
      }
      g = readerx_.getDataAt(i+1);
    }
    numRows_x_ = (int)readerx_.getLineNumber();
    g = reader_.getData();
    for (int i = 1;  g != null; i++) {
      String id = (String) g.getDataAt(0);
      BitSet va = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 0);
      BitSet va_thr = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 1);
      va = setPhenotype(va);
      va_thr = setPhenotype(va_thr);
      if (!haveGoodDynamicRange(va_thr)) {
        //System.out.println(i + "\tbd" + "\t" + id);
        g = reader_.getData();
        continue;
      }
      int group = findBias(va, va_thr);
      //System.out.println(i + "\t" + group + "\t" + id);
      if (group == 0) {
        low_.add(id);
      }
      if (group == 1) {
        high_.add(id);
      }
      if (group == -1) {
        Integer loc = new Integer(i + numRows_x_);
        balancedy_.add(id);
        balanced_map_.put(loc, new Integer(index));
        index++;
      }
      g = reader_.getData();
    }
    out_.writeList("low", low_);
    out_.writeList("high", high_);
    out_.writeList("balancedx", balancedx_);
    out_.writeList("balancedy", balancedy_);
  }

  public void performBlockAnalysis() throws IOException {
    out_.startMatrix(balancedy_.size(), 3);
    GeneData gb1 = readerx_.getDataAt(0);
    for (int b1 = 0;  gb1 != null; b1+=blocksize_) {
      BitSet[] ba1 = new BitSet[blocksize_];
      BitSet[] ba1_thr = new BitSet[blocksize_];
      GeneData ga = readerx_.getDataAt(b1);
      for (int i = b1; ga != null && i < (b1+blocksize_); i++) {
        ga = readerx_.getDataAt(i);
        if (ga == null) {
            break;
        }
        BitSet va = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 0);
        BitSet va_thr = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 1);
        va = setPhenotype(va);
        va_thr = setPhenotype(va_thr);
        ba1[i-b1] = va;
        ba1_thr[i-b1] = va_thr;
      }
      int start = b1;
      BitSet va = ba1[0];
      for (int i = b1; va != null && i < (b1+blocksize_); i++) {
        System.out.println(i);
        va = ba1[i-b1];
        BitSet va_thr = ba1_thr[i-b1];
        if (!haveGoodDynamicRange(va_thr) || 
            !balanced_map_.containsKey(new Integer(i))) {
          if (i < (b1-1 +blocksize_)) {
            va = ba1[i+1-b1];
          }
          continue;
        }
        BitSet vb = ba1[i+1];
        for (int j = i+1; vb != null && j < (b1+blocksize_); j++) {
          vb = ba1[j-b1];
          BitSet vb_thr = ba1_thr[j-b1];
          if (!haveGoodDynamicRange(vb_thr) || 
              !balanced_map_.containsKey(new Integer(j))) {
            if (j < (b1-1 +blocksize_)) {
              vb = ba1[j+1-b1];
            }
            continue;
          }
          performSinglePairAnalysis(i, j, va, va_thr, vb, vb_thr);
          if (j < (b1-1 +blocksize_)) {
            vb = ba1[j+1-b1];
          }
        }
        reader_.begin();
        GeneData gb = reader_.getData();
        for (int j = 1; gb != null ; j++) {
          Integer loc = new Integer(j+numRows_x_);
          vb = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 0);
          BitSet vb_thr = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 1);
          vb = setPhenotype(vb);
          vb_thr = setPhenotype(vb_thr);
          if (!haveGoodDynamicRange(vb_thr) ||
              !balanced_map_.containsKey(loc)) {
            gb = reader_.getData();
            continue;
          }
          performSinglePairAnalysis(i, loc.intValue(), va, va_thr, vb, vb_thr);
          gb = reader_.getData();
        }
        if (i < (b1-1 +blocksize_)) {
          va = ba1[i+1-b1];
        }
      } // end va
      gb1 = readerx_.getDataAt(b1+blocksize_);
    } // end gb1
  }

  public double[] getErrorProbStats(
      BitSet a, BitSet a_thr, BitSet b, BitSet b_thr, int debug) { 
    double[] res = new double[4];
    res[0] = res[1] = res[2] = res[3] = 1.0;
    if (a.length() == 0 || b.length() == 0) {
      return res;
    }
    BitSet thrBits = (BitSet) a_thr.clone();
    thrBits.and(b_thr);
    BitSet tmp = (BitSet) thrBits.clone();
    BitSet v1 = (BitSet) a.clone();
    v1.or(b);
    tmp.andNot(v1);
    int c0 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    v1 = (BitSet) b.clone();
    v1.andNot(a);
    tmp.and(v1);
    int c1 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    v1 = (BitSet) a.clone();
    v1.andNot(b);
    tmp.and(v1);
    int c2 = tmp.cardinality();
    tmp = (BitSet) thrBits.clone();
    v1 = (BitSet) a.clone();
    v1.and(b);
    tmp.and(v1);
    int c3 = tmp.cardinality();

    int total = c0 + c1 + c2 + c3;
    if (total <= 0) {
      return res;
    }

    res[0] = ((c0 + c1) * (c0 + c2)/total - c0 + 1)/Math.sqrt((c0 + c1) * (c0 + c2)/total + 1);
    res[1] = ((c1 + c0) * (c1 + c3)/total - c1 + 1)/Math.sqrt((c1 + c0) * (c1 + c3)/total + 1);
    res[2] = ((c2 + c0) * (c2 + c3)/total - c2 + 1)/Math.sqrt((c2 + c0) * (c2 + c3)/total + 1);
    res[3] = ((c3 + c1) * (c3 + c2)/total - c3 + 1)/Math.sqrt((c3 + c1) * (c3 + c2)/total + 1);
    if (debug > 0) {
      System.out.println(c0 + "\t" + c1 + "\t" + c2 + "\t" + c3 + "\t" + total);
      System.out.println(res[0] + "\t" + res[1] + "\t" + res[2] + "\t" + res[3]);
    }
    if (res[0] > stat_threshold_) {
      res[0] = 0.5 * c0 / (c0 + c1 + 1) + 0.5 * c0/(c0 + c2 + 1);
    }
    else {
      res[0] = 1;
    }
    if (res[1] > stat_threshold_) {
      res[1] = 0.5 * c1 / (c1 + c0 + 1) + 0.5 * c1/(c1 + c3 + 1);
    }
    else {
      res[1] = 1;
    }
    if (res[2] > stat_threshold_) {
      res[2] = 0.5 * c2 / (c2 + c0 + 1) + 0.5 * c2/(c2 + c3 + 1);
    }
    else {
      res[2] = 1;
    }
    if (res[3] > stat_threshold_) {
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

  public void performSinglePairAnalysis(int i, int j, 
      BitSet va, BitSet va_thr, BitSet vb, BitSet vb_thr) throws IOException {
    double[] p = getErrorProbStats(va, va_thr, vb, vb_thr, 0);
    int code = 0;
    if (p[0] <= threshold_) { code = 1; }
    if (p[1] <= threshold_) { code = 2; }
    if (p[2] <= threshold_) { code = 3; }
    if (p[3] <= threshold_) { code = 4; }
    if (p[1] <= threshold_ && p[2] <= threshold_) {
      code = 5;
    }
    if (p[0] <= threshold_ && p[3] <= threshold_) {
      code = 6;
    }
    if (code > 0) {
      Integer a = balanced_map_.get(new Integer(i));
      Integer b = balanced_map_.get(new Integer(j));
      //System.out.println(code + "\t" + i + "\t" + j + "\t" + a + "\t" + b);
      out_.setBitMatrix(a.intValue(), b.intValue(), code);
    }
    else {
      //Integer a = balanced_map_.get(new Integer(i));
      //Integer b = balanced_map_.get(new Integer(j));
      //System.out.println(code + "\t" + i + "\t" + j + "\t" + a + "\t" + b);
    }
  }

  public static String formatString(String format, double x) {
    MessageFormat mf = new MessageFormat(
        "{0,number," + format + "}");
    Object[] objs = {new Double(x)};
    return mf.format(objs);
  }

  public static void main(String args[]) throws Exception {
    if (args.length < 1) {
      System.out.println("Arguments : <cmd:bitMatrix/list/pairs> args");
      System.exit(1);
    }
    if (args[0].equals("listMatrix")) {
      if (args.length < 8) {
        System.out.println("Arguments : <cmd> <ofile> <bvfilex> <bvfile> <phfile> <phid>  pvalue statThr listFile");
        System.exit(1);
      }
      BooleanSubAnalysis ana = new BooleanSubAnalysis(args[2], args[3], args[1], args[4], args[5]);
      ana.setThreshold(Double.parseDouble(args[6]));
      ana.setStatThreshold(Double.parseDouble(args[7]));
      //ana.performAnalysis(); 
      ana.setGeneList(args[8]);
      ana.performListAnalysis(); 
    }
    if (args[0].equals("bitMatrix")) {
      if (args.length < 7) {
        System.out.println("Arguments : <cmd> <ofile> <bvfilex> <bvfile> <phfile> <phid> pvalue statThr ");
        System.exit(1);
      }
      BooleanSubAnalysis ana = new BooleanSubAnalysis(args[2], args[3], args[1], args[4], args[5]);
      ana.setThreshold(Double.parseDouble(args[6]));
      ana.setStatThreshold(Double.parseDouble(args[7]));
      ana.performAnalysis(); 
    }
    if (args[0].equals("bitMatrix1")) {
      if (args.length < 9) {
        System.out.println("Arguments : <cmd> <ofile> <bvfilex> <bvfile> <phfile> <phid> pvalue statThr single_cutoff single_thr");
        System.exit(1);
      }
      BooleanSubAnalysis ana = new BooleanSubAnalysis(args[2], args[3], args[1], args[4], args[5]);
      ana.setThreshold(Double.parseDouble(args[6]));
      ana.setStatThreshold(Double.parseDouble(args[7]));
      ana.setSingleCutoff(Integer.parseInt(args[8]));
      ana.setSingleThreshold(Double.parseDouble(args[9]));
      ana.performAnalysis(); 
    }
  }

}

