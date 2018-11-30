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

public class BooleanAnalysis {

  String output_filename_;
  String bv_filename_;
  String ph_filename_;
  String phid_;

  double threshold_;
  double stat_threshold_;
  double single_threshold_;
  int single_cutoff_;
  int blocksize_;

  PCLFileReader reader_;
  NetworkFile out_;
  int fileType_;
  BitSet phenotype_;
  int numArrays_;
  int numPhArrays_;

  Vector<String> high_;
  Vector<String> low_;
  Vector<String> balanced_;
  HashMap<Integer, Integer> balanced_map_;

  String geneList_;
  HashMap<String, Integer> geneMaps_;

  public BooleanAnalysis(String bvfile, String ofile) {
    output_filename_ = ofile;
    bv_filename_ = bvfile;
    ph_filename_ = null;
    phid_ = null;
    init();
  }

  public BooleanAnalysis(String bvfile, String ofile, String phfile, String id) {
    output_filename_ = ofile;
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
    fileType_ = NetworkFile.FILE_1_1;
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
    performListSingleAnalysis(0);
    out_.startMatrix(0, 3);
    finish();
  }
  public void performListAnalysis() throws IOException {
    beginAnalysis();
    readGeneList();
    performListSingleAnalysis(0);
    performListPairsAnalysis(0);
    finish();
  }
  public void performListDebugAnalysis() throws IOException {
    beginAnalysis();
    readGeneList();
    performListSingleAnalysis(1);
    performListPairsAnalysis(1);
    finish();
  }

  public void writeListPairs() throws IOException {
    fileType_ = NetworkFile.PAIRS;
    beginAnalysis();
    readGeneList();
    performListSingleAnalysis(0);
    performListPairsAnalysis(0);
    finish();
  }

  public void readGeneList() throws IOException {
    PCLFileReader listr = new PCLFileReader(geneList_);
    listr.begin();
    HashSet<String> geneIds = new HashSet<String>();
    geneMaps_ = new HashMap<String, Integer>();
    GeneData g = listr.getData();
    for (int i = 0;  g != null; i++) {
      String id = (String) g.getDataAt(0);
      id = id.trim();
      geneIds.add(id);
      g = listr.getData();
    }
    g = reader_.getDataAt(0);
    for (int i = 0;  g != null; i++) {
      String id = (String) g.getDataAt(0);
      id = id.trim();
      if (geneIds.contains(id)) {
        geneMaps_.put(id, new Integer(i));
      }
      g = reader_.getDataAt(i+1);
    }
    System.out.println("Found " + geneMaps_.size() + " probes");
  }

  public void performListSingleAnalysis(int debug) throws IOException {
    high_ = new Vector<String>();
    low_ = new Vector<String>();
    balanced_ = new Vector<String>();
    balanced_map_ = new HashMap<Integer, Integer>();
    int index = 0;
    Iterator<String> keys = geneMaps_.keySet().iterator();
    while (keys.hasNext()) {
      String id = (String) keys.next();
      Integer loc = geneMaps_.get(id);
      GeneData g = reader_.getDataAt(loc.intValue());
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
        balanced_.add(id);
        balanced_map_.put(loc, new Integer(index));
        index++;
      }
    }
    out_.writeList("low", low_);
    out_.writeList("high", high_);
    out_.writeList("balanced", balanced_);
  }

  public void performListPairsAnalysis(int debug) throws IOException {
    out_.startMatrix(balanced_.size(), 3);
    Enumeration<String> alist = balanced_.elements();
    while (alist.hasMoreElements()) {
      String aid = alist.nextElement();
      Integer aloc = geneMaps_.get(aid);
      GeneData ga = reader_.getDataAt(aloc.intValue());
      BitSet va = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 0);
      BitSet va_thr = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 1);
      va = setPhenotype(va);
      va_thr = setPhenotype(va_thr);
      if (!haveGoodDynamicRange(va_thr)) {
        continue;
      }
      Enumeration<String> blist = balanced_.elements();
      while (blist.hasMoreElements()) {
        String bid = blist.nextElement();
        Integer bloc = geneMaps_.get(bid);
        GeneData gb = reader_.getDataAt(bloc.intValue());
        BitSet vb = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 0);
        BitSet vb_thr = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 1);
        vb = setPhenotype(vb);
        vb_thr = setPhenotype(vb_thr);
        if (!haveGoodDynamicRange(vb_thr)) {
          continue;
        }
        performSinglePairAnalysis(aloc.intValue(), bloc.intValue(), va, va_thr, vb, vb_thr, debug);
      }
    }
  }

  public void performAnalysis() throws IOException {
    beginAnalysis();
    performSingleAnalysis();
    performBlockAnalysis();
    finish();
  }

  public void finish() throws IOException {
    out_.close();
  }

  public void beginAnalysis() throws IOException {
    reader_ = new PCLFileReader(bv_filename_);
    PCLFileReader.CACHE_SIZE = 0;
    reader_.beginRandomAccess();
    if (fileType_ == NetworkFile.FILE_1_1) {
      out_ = new BitMatrixNetworkSimple(output_filename_);
    }
    if (fileType_ == NetworkFile.FILE_1_0) {
      out_ = new BitMatrixNetworkFile(output_filename_);
    }
    if (fileType_ == NetworkFile.PAIRS) {
      out_ = new NetworkPairsFile(output_filename_);
    }
    out_.writeHeader();
    GeneData g = reader_.getDataAt(1);
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
    balanced_ = new Vector<String>();
    balanced_map_ = new HashMap<Integer, Integer>();
    int index = 0;
    GeneData g = reader_.getDataAt(0);
    for (int i = 0;  g != null; i++) {
      String id = (String) g.getDataAt(0);
      BitSet va = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 0);
      BitSet va_thr = BitSetUtils.stringToBitSet((String) g.getDataAt(2), 1);
      va = setPhenotype(va);
      va_thr = setPhenotype(va_thr);
      if (!haveGoodDynamicRange(va_thr)) {
        //System.out.println(i + "\tbd");
        g = reader_.getDataAt(i+1);
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
        balanced_.add(id);
        balanced_map_.put(new Integer(i), new Integer(index));
        index++;
      }
      g = reader_.getDataAt(i+1);
    }
    out_.writeList("low", low_);
    out_.writeList("high", high_);
    out_.writeList("balanced", balanced_);
  }

  public void performBlockAnalysis() throws IOException {
    out_.startMatrix(balanced_.size(), 3);
    GeneData gb1 = reader_.getDataAt(0);
    for (int b1 = 0;  gb1 != null; b1+=blocksize_) {
      BitSet[] ba1 = new BitSet[blocksize_];
      BitSet[] ba1_thr = new BitSet[blocksize_];
      GeneData ga = reader_.getDataAt(b1);
      for (int i = b1; ga != null && i < (b1+blocksize_); i++) {
        ga = reader_.getDataAt(i);
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
      GeneData gb2 = reader_.getDataAt(start);
      for (int b2 = start;  gb2 != null; b2+=blocksize_) {
        System.out.println("Block = (" + b1 + ", " + b2 + ")");
        BitSet[] ba2 = new BitSet[blocksize_];
        BitSet[] ba2_thr = new BitSet[blocksize_];
        GeneData gb = reader_.getDataAt(b2);
        for (int i = b2; gb != null && i < (b2+blocksize_); i++) {
          gb = reader_.getDataAt(i);
          if (gb == null) {
            break;
          }
          BitSet vb = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 0);
          BitSet vb_thr = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 1);
          vb = setPhenotype(vb);
          vb_thr = setPhenotype(vb_thr);
          ba2[i-b2] = vb;
          ba2_thr[i-b2] = vb_thr;
        }
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
          BitSet vb = ba2[0];
          for (int j = b2; vb != null && j < (b2 + blocksize_); j++) {
            vb = ba2[j-b2];
            BitSet vb_thr = ba2_thr[j-b2];
            if (!haveGoodDynamicRange(vb_thr) || (j <= i) ||
                !balanced_map_.containsKey(new Integer(j))) {
              if (j < (b2 -1 + blocksize_)) {
                vb = ba2[j+1-b2];
              }
              continue;
            }
            performSinglePairAnalysis(i, j, va, va_thr, vb, vb_thr, 0);
            if (j < (b2 -1 + blocksize_)) {
              vb = ba2[j+1-b2];
            }
          }
          if (i < (b1-1 +blocksize_)) {
            va = ba1[i+1-b1];
          }
        } // end ga
        gb2 = reader_.getDataAt(b2+blocksize_);
      } // end gb2
      gb1 = reader_.getDataAt(b1+blocksize_);
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
      BitSet va, BitSet va_thr, BitSet vb, BitSet vb_thr, int debug) throws IOException {
    if (fileType_ == NetworkFile.PAIRS) {
      double[] p = getErrorProbStats(va, va_thr, vb, vb_thr, debug);
      double pvalue = 1.0;
      int code = 0;
      if (p[0] <= threshold_) { code = 1; pvalue = p[0]; }
      if (p[1] <= threshold_) { code = 2; pvalue = p[1]; }
      if (p[2] <= threshold_) { code = 3; pvalue = p[2]; }
      if (p[3] <= threshold_) { code = 4; pvalue = p[3]; }
      if (p[1] <= threshold_ && p[2] <= threshold_) {
        code = 5; pvalue = p[1] > p[2] ?  p[1] : p[2];
      }
      if (p[0] <= threshold_ && p[3] <= threshold_) {
        code = 6; pvalue = p[0] > p[3] ?  p[0] : p[3];
      }
      if (code > 0) {
        out_.writePair(i, j, code, pvalue, null);
      }
      return;
    }
    double[] p = getErrorProbStats(va, va_thr, vb, vb_thr, debug);
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

  public void writePairs() throws IOException {
    fileType_ = NetworkFile.PAIRS;
    beginAnalysis();
    performSingleAnalysis();
    performBlockAnalysis();
    finish();
  }

  public static HashMap<Long, Integer> getListMap(
      PCLFileReader data,
      LinkedList<Long> listx,
      HashMap<Integer, BitSet> cache,
      HashMap<Integer, BitSet> cache_thr) throws IOException {
    ListIterator<Long> itrx = listx.listIterator();
    HashMap<Long, Integer> mapx = new HashMap<Long, Integer>();
    while (itrx.hasNext()) {
      Long lx = itrx.next();
      int a = lx.intValue();
      BitSet va_thr = null;
      if (cache_thr.containsKey(new Integer(a))) {
        va_thr = cache_thr.get(new Integer(a));
      }
      if (va_thr == null) {
        GeneData ga = data.getDataAt(a);
        va_thr = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 1);
        cache_thr.put(new Integer(a), va_thr);
      }
      int outside = va_thr.cardinality();
      mapx.put(lx, new Integer(outside));
      if (cache_thr.size() > 10000) {
        System.out.println("Vector cache clear");
        cache.clear();
        cache_thr.clear();
      }
    }
    return mapx;
  }

  public void writeCommonPairs(LinkedList<String> list) throws IOException {
    String pairfile = list.removeFirst();
    String cfile = list.removeFirst();
    String homologFile = null;
    HomologData homolog_data = new HomologData();
    String srcOrg = null;
    String destOrg = null;

    if (list.size() > 0) {
       homologFile = list.removeFirst();
       homolog_data.setFilename(homologFile, HomologData.EUGENE);
       srcOrg = list.removeFirst();
       destOrg = list.removeFirst();
    }
    homolog_data.parse();

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
      if (!b_ids.containsKey(name)) {
        b_ids.put(name, new LinkedList<Long>());
      }
      LinkedList<Long> val = b_ids.get(name);
      val.add(new Long(i));
    }
    GeneData g = data.getDataAt(1);
    String eweight = (String) g.getDataAt(2);
    numArrays_ = eweight.length();
    numPhArrays_ = numArrays_;

    PCLFileReader pcldata = new PCLFileReader(bv_filename_);
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
      indexhash.put(new Integer(i), name);
    }

    TABFileReader pairdata = new TABFileReader(pairfile);
    pairdata.begin();
    BufferedWriter out = new BufferedWriter(new FileWriter(output_filename_));
    GeneData pair = pairdata.getHeader();
    HashMap<Integer, BitSet> cache = new HashMap<Integer, BitSet>();
    HashMap<Integer, BitSet> cache_thr = new HashMap<Integer, BitSet>();

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
      LinkedList<String> namexl = homolog_data.getHomolog(namex, srcOrg, destOrg);
      LinkedList<String> nameyl = homolog_data.getHomolog(namey, srcOrg, destOrg);
      int foundx = 0;
      ListIterator<String> nitrx = namexl.listIterator();
      while (nitrx.hasNext()) {
        String nx = nitrx.next();
        if ( b_ids.containsKey(nx) ) {
          foundx = 1;
        }
      }
      int foundy = 0;
      ListIterator<String> nitry = nameyl.listIterator();
      while (nitry.hasNext()) {
        String ny = nitry.next();
        if ( b_ids.containsKey(ny) ) {
          foundy = 1;
        }
      }

      if ( foundx == 0 || foundy == 0 ) {
        continue;
      }
      // Collect Indices
      LinkedList<Long> listx = new LinkedList<Long>();
      LinkedList<Long> listy = new LinkedList<Long>();
      nitrx = namexl.listIterator();
      while (nitrx.hasNext()) {
        String nx = nitrx.next();
        if ( b_ids.containsKey(nx) ) {
          listx.addAll(b_ids.get(nx));
        }
      }
      nitry = nameyl.listIterator();
      while (nitry.hasNext()) {
        String ny = nitry.next();
        if ( b_ids.containsKey(ny) ) {
          listy.addAll(b_ids.get(ny));
        }
      }
      HashMap<Long, Integer> mapx = getListMap(data, listx, cache, cache_thr);
      HashMap<Long, Integer> mapy = getListMap(data, listy, cache, cache_thr);

      class ListComparator implements Comparator<Long> {
        HashMap<Long, Integer> map_;
        public ListComparator(HashMap<Long, Integer> map) {
          map_ = map;
        }
        public int compare(Long o1, Long o2) {
          Integer a = map_.get(o1);
          Integer b = map_.get(o2);
          return b.compareTo(a);
        }
      }

      // Sort indices
      Collections.sort(listx, new ListComparator(mapx));
      Collections.sort(listy, new ListComparator(mapy));

      mapx.clear();
      mapy.clear();

      //System.out.println(listx.size() + "\t" + namex);
      //System.out.println(listy.size() + "\t" + namey);

      boolean found = false;
      int a = -1 , b = -1;
      double pvalue = 1.0;
      ListIterator<Long> itrx = listx.listIterator();

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
            va = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 0);
            va_thr = BitSetUtils.stringToBitSet((String) ga.getDataAt(2), 1);
            cache.put(new Integer(a), va);
            cache_thr.put(new Integer(a), va_thr);
          }
          if (vb == null) {
            GeneData gb = data.getDataAt(b);
            vb = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 0);
            vb_thr = BitSetUtils.stringToBitSet((String) gb.getDataAt(2), 1);
            cache.put(new Integer(b), vb);
            cache_thr.put(new Integer(b), vb_thr);
          }
          if (cache.size() > 10000) {
            System.out.println("Vector cache clear");
            cache.clear();
            cache_thr.clear();
          }
          if (!haveGoodDynamicRange(va_thr) || 
              !haveGoodDynamicRange(vb_thr)) {
            continue;
          }
          double[] p = getErrorProbStats(va, va_thr, vb, vb_thr, 0);
          int code = 0;
          if (p[0] <= threshold_) { code = 1; pvalue = p[0]; }
          if (p[1] <= threshold_) { code = 2; pvalue = p[1]; }
          if (p[2] <= threshold_) { code = 3; pvalue = p[2]; }
          if (p[3] <= threshold_) { code = 4; pvalue = p[3]; }
          if (p[1] <= threshold_ && p[2] <= threshold_) {
            code = 5; pvalue = p[1] > p[2] ?  p[1] : p[2];
          }
          if (p[0] <= threshold_ && p[3] <= threshold_) {
            code = 6; pvalue = p[0] > p[3] ?  p[0] : p[3];
          }
          if (code == q) {
            found = true;
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

  public static String formatString(String format, double x) {
    MessageFormat mf = new MessageFormat(
        "{0,number," + format + "}");
    Object[] objs = {new Double(x)};
    return mf.format(objs);
  }

  public static void main(String args[]) throws Exception {
    if (args.length < 1) {
      System.out.println("Arguments : <cmd:bitMatrix/list/pairs> args");
      System.out.println("Arguments : <cmd:commonPairs> args");
      System.exit(1);
    }
    if (args[0].equals("listMatrix")) {
      if (args.length < 8) {
        System.out.println("Arguments : <cmd> <ofile> <bvfile> <phfile> <phid>  pvalue statThr listFile");
        System.exit(1);
      }
      BooleanAnalysis ana = new BooleanAnalysis(args[2], args[1], args[3], args[4]);
      ana.setThreshold(Double.parseDouble(args[5]));
      ana.setStatThreshold(Double.parseDouble(args[6]));
      //ana.performAnalysis(); 
      ana.setGeneList(args[7]);
      ana.performListAnalysis(); 
    }
    if (args[0].equals("bitMatrix")) {
      if (args.length < 7) {
        System.out.println("Arguments : <cmd> <ofile> <bvfile> <phfile> <phid> pvalue statThr ");
        System.exit(1);
      }
      BooleanAnalysis ana = new BooleanAnalysis(args[2], args[1], args[3], args[4]);
      ana.setThreshold(Double.parseDouble(args[5]));
      ana.setStatThreshold(Double.parseDouble(args[6]));
      ana.performAnalysis(); 
    }
    if (args[0].equals("pairs")) {
      if (args.length < 7) {
        System.out.println("Arguments : <cmd> <ofile> <bvfile> <phfile> <phid> pvalue statThr ");
        System.exit(1);
      }
      BooleanAnalysis ana = new BooleanAnalysis(args[2], args[1], args[3], args[4]);
      ana.setThreshold(Double.parseDouble(args[5]));
      ana.setStatThreshold(Double.parseDouble(args[6]));
      ana.writePairs();
    }
    if (args[0].equals("commonPairs")) {
      if (args.length < 7) {
        System.out.println("Arguments : <cmd> <ofile> <bvfile> pvalue statThr <pairfile> <bvfile> <homologFile>");
        System.exit(1);
      }
      LinkedList<String> list = new LinkedList<String>(Arrays.asList(args));
      String cmd = list.removeFirst();
      String ofile = list.removeFirst();
      String bvfile = list.removeFirst();
      String pvalue = list.removeFirst();
      String statThr = list.removeFirst();
      BooleanAnalysis ana = new BooleanAnalysis(bvfile, ofile);
      ana.setThreshold(Double.parseDouble(pvalue));
      ana.setStatThreshold(Double.parseDouble(statThr));
      ana.writeCommonPairs(list);
    }
  }

}

