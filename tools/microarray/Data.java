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

public class Data implements Cloneable {

  String name_;         // Name of the data
  static int name_id_  = 0; // name id of the data
  int numArrays_;       // Number of Arrays used in the Microarray Experiment
  int numGenes_;        // Number of genes
  int numGeneHeader_;   // PCL File : CLID, EWEIGHT
  int numArrayHeader_;  // PCL File : CLID, NAME, GWEIGHT
  int id_;              // index of the primary gene ID

  GeneNameScheme nameScheme_;

  GeneData[] data_;

  Double[] timepoints_; // Adjust stepSearch
  int breakpoint_;      // For only two step analysis

  void init_(int nA, int nG, int nGh, int nAh, GeneData[] d) {
    name_ = "Data-"+name_id_;
    name_id_ ++;
    numArrays_ = nA;
    numGenes_ = nG;
    numGeneHeader_ = nGh;
    numArrayHeader_ = nAh;
    id_ = 0;

    data_ = d;
    nameScheme_ = new GeneNameScheme();
    timepoints_ = new Double[numArrays_];
    for (int i =0; i < timepoints_.length; i++) {
      timepoints_[i] = new Double(i);
    }
    breakpoint_ = -1;
  }

  public Data(int nA, int nG, int nGh, int nAh, GeneData[] d) {
    init_(nA, nG, nGh, nAh, d);
  }

  public Data(int nA, int nG, int nGh, int nAh, Vector<GeneData> dv) {
    GeneData[] d = new GeneData[dv.size()];
    for (int i =0; i < dv.size(); i++) {
      d[i] = dv.get(i);
    }
    init_(nA, nG, nGh, nAh, d);
  }

  public void setGeneNameScheme(GeneNameScheme scheme) { nameScheme_ = scheme; }
  public GeneNameScheme getGeneNameScheme() { return nameScheme_; }

  public void setName(String n) { name_ = n; }
  public String getName() { return name_; }
  public void setID(int id) { id_ = id; }
  public int getID() { return id_; }

  public int getNumArrays     () { return numArrays_;     }
  public int getNumGenes      () { return numGenes_;      }
  public int getNumGeneHeader () { return numGeneHeader_; }
  public int getNumArrayHeader() { return numArrayHeader_;}
  public int getNumColumns    () { return numArrays_+numArrayHeader_;}
  public int getNumRows       () { return numGenes_+numGeneHeader_;}

  public int getNumMissingPoints() { return nameScheme_.getNumMissingPoints();}
  public GeneData getGeneData(int index) { return data_[index];}
  public void setGeneData(GeneData s, int index) { data_[index] = s;}
  public Double[] getTimepoints() { return timepoints_;}
  public void setTimepoints(Double[] t) { timepoints_ = t; }
  public void setBreakpoint(int t) { breakpoint_ = t; }

  /**
   * Compute the step positions that need to be searched
   *  from the timepoints specification
   */
  public int[] getStepSearch() {
    int count = 1;
    Double last = timepoints_[0];
    for (int i =1; i < numArrays_; i++) {
      if (!timepoints_[i].equals(last)) {
        count++;
      }
      last = timepoints_[i];
    }
    int[] st = new int[count];
    int index = 0;
    last = timepoints_[0];
    for (int i =1; i < numArrays_; i++) {
      if (!timepoints_[i].equals(last)) {
        st[index++] = (i-1) + numArrayHeader_;
      }
      last = timepoints_[i];
    }
    st[index++] = numArrays_ - 1 + numArrayHeader_;
    if (breakpoint_ > 0) {
      count = 0;
      for (int i =0; i < st.length; i++) {
        if (st[i] < breakpoint_) {
          count++;
        }
      }
      int[] tmpst = new int[count];
      index = 0;
      for (int i =0; i < st.length; i++) {
        if (st[i] < breakpoint_) {
          tmpst[index++] = st[i];
        }
      }
      st = tmpst;
    }
    return st;
  }

  /**
   *  StepSearch for the Two Step Analysis
   */
  public int[] getStepSearch1() {
    int count = 1;
    Double last = timepoints_[0];
    for (int i =1; i < numArrays_; i++) {
      if (!timepoints_[i].equals(last)) {
        count++;
      }
      last = timepoints_[i];
    }
    int[] st = new int[count];
    int index = 0;
    last = timepoints_[0];
    for (int i =1; i < numArrays_; i++) {
      if (!timepoints_[i].equals(last)) {
        st[index++] = (i-1);
      }
      last = timepoints_[i];
    }
    st[index++] = numArrays_ - 1;
    int breakpoint = breakpoint_ - numArrayHeader_;
    if (breakpoint > 0) {
      count = 0;
      for (int i =0; i < st.length; i++) {
        if (st[i] >= breakpoint) {
          count++;
        }
      }
      int[] tmpst = new int[count];
      index = 0;
      for (int i =0; i < st.length; i++) {
        if (st[i] >= breakpoint) {
          tmpst[index++] = st[i];
        }
      }
      st = tmpst;
    }
    return st;
  }

  public Object clone() {
    GeneData[] obj = new GeneData[data_.length];
    for (int i = 0; i < obj.length; i++) {
        obj[i] = (GeneData)data_[i].clone();
    }
    Data res = new Data(numArrays_, numGenes_, numGeneHeader_, numArrayHeader_,
        obj);
    res.nameScheme_ = nameScheme_;
    return res;
  }

  public String toString() {
    String res = "Header : " + numGeneHeader_ + "x" + numArrayHeader_ + "\n" +
        "Data : " + numGenes_ + "x" + numArrays_ + "\n";
    if (nameScheme_ != null) {
        res += nameScheme_.toString();
    }
    return res;
  }

  public void print() {
    System.out.println(toString());
  }

  public void setRange(String range) {
    if (range == null) {
      return;
    }
    try {
      String[] list = range.split(":");
      int numArrayHeader = numArrayHeader_;
      int numArrays = numArrays_;
      if (list.length > 0) {
        numArrayHeader = Integer.parseInt(list[0]);
        numArrays -= (numArrayHeader - numArrayHeader_);
      }
      if (list.length > 1) {
        numArrays = Integer.parseInt(list[1]) - numArrayHeader + 1;
      }
      numArrayHeader_ = numArrayHeader;
      numArrays_ = numArrays;
      System.out.println("Range set as " + numArrayHeader_ + ":" +
        (getNumColumns() - 1) );
    }
    catch(Exception e) {
      System.out.println(" *Warning : Malformed range");
    }
  }

  public void setColRange(String range) {
    setRange(range);
  }

  public void setRowRange(String range) {
    if (range == null) {
      return;
    }
    try {
      String[] list = range.split(":");
      int numGeneHeader = numGeneHeader_;
      int numGenes = numGenes_;
      if (list.length > 0) {
        numGeneHeader = Integer.parseInt(list[0]);
        numGenes -= (numGeneHeader - numGeneHeader_);
      }
      if (list.length > 1) {
        numGenes = Integer.parseInt(list[1]) - numGeneHeader + 1;
      }
      numGeneHeader_ = numGeneHeader;
      numGenes_ = numGenes;
      System.out.println("Row range set as " + numGeneHeader_ + ":" +
        (getNumRows() - 1) );
    }
    catch(Exception e) {
      System.out.println(" *Warning : Malformed range");
    }
  }

  public String getGenesAt(int index) {
    int arrayIndex = nameScheme_.getArrayIndex();
    GeneData gene = data_[index];
    String str = (String) gene.getDataAt(arrayIndex);
    String name = nameScheme_.getGene(str);
    return name;
  }

  public int[] getNullOrder() {
    int[] order = new int[numGenes_];
    for (int i = numGeneHeader_; i < (numGenes_+numGeneHeader_); i++) {
      order[i-numGeneHeader_] = i;
    }
    return order;
  }

  public String[] getGenes(int[] order) {
    int arrayIndex = nameScheme_.getArrayIndex();
    String[] res = new String[order.length];
    for (int i =0; i < order.length; i++) {
        GeneData gene = data_[order[i]];
        String str = (String) gene.getDataAt(arrayIndex);
        res[i] = nameScheme_.getGene(str);
    }
    return res;
  }

  public void convertDoubles() throws ArrayException {
    for (int i = numGeneHeader_; i < (numGenes_+numGeneHeader_); i++) {
        data_[i].convertDouble(numArrayHeader_, getNumColumns()-1);
    }
  }

  public void restrictGenes(int[] order) {
    numGenes_ = order.length;
    GeneData[] res = new GeneData[getNumRows()];
    for (int i = 0; i < getNumGeneHeader(); i++) {
        res[i] = data_[i];
    }
    for (int i =0; i < order.length; i++) {
        res[i+getNumGeneHeader()] = data_[order[i]];
    }
    data_ = res;
  }

  public void reduceLog() throws ArrayException {
    for (int i = numGeneHeader_; i < numGenes_+numGeneHeader_; i++) {
        data_[i].reduceLog(numArrayHeader_, getNumColumns()-1);
    }
  }

  public void scramble(Random random) {
    for (int i = numGeneHeader_; i < numGenes_+numGeneHeader_; i++) {
        data_[i].scramble(random, numArrayHeader_, getNumColumns()-1);
    }
  }

  // Parse the filename and retrieve filename
  //   Example : full:label.ps
  //         Output : label.ps
  public static String getFileName(String file) {
    if (file.indexOf(':') == -1) {
        return file;
    }
    return file.split(":")[1];
  }

  // Parse the filename and retrieve filetag
  //   Example : full:label.ps
  //         Output : full
  public static String getFileTag(String file) {
    if (file.indexOf(':') == -1) {
        return "All";
    }
    return file.split(":")[0];
  }

  public static Data getArtificialData(Random random, int numgenes, int numarrays) {
    GeneData[] genes = new GeneData[numgenes+2];
    Object[] first = new Object[numarrays+3];
    Object[] second = new Object[numarrays+3];
    first[0] = "CLID"; first[1] = "Name"; first[2] = "GWEIGHT";
    second[0] = "EWEIGHT"; second[1] = ""; second[2] = "";
    for (int i = 0; i < numarrays; i++) {
        first[i+3] = "" + i;
        second[i+3] = "1";
    }
    genes[0] = new GeneData(first);
    genes[1] = new GeneData(second);

    for (int i = 0; i < numgenes ; i++) {
      Object[] d = new Object[numarrays+3];
      d[0] = "" + i ; d[1] = "" + i; d[2] = "1";
      for (int j = 0; j < numarrays; j++) {
        d[j+3] = new Double(random.nextGaussian());
      }
      genes[i+2] = new GeneData(d);
    }
    Data res = new Data(numarrays, numgenes, 2, 3, genes);
    return res;
  }

  public static Data intersectData(Data a, Data b) {
    return mergeData(a, b, true, false);
  }

  public static Data intersectDataSelect(Data a, Data b) {
    return mergeData(a, b, true, true);
  }

  public static Data unionData(Data a, Data b) {
    return mergeData(a, b, false, false);
  }

  public static Data unionDataSelect(Data a, Data b) {
    return mergeData(a, b, false, true);
  }

  /**
   *  @param a - First Data
   *  @param b - Second Data
   *  @param intersect - if (true) intersect else do union
   *  @param select - if (true) select the first data only
   *  @return - new Data with the merging
   */
  public static Data mergeData(Data a, Data b, 
      boolean intersect, boolean select) 
  {
    int numArrays = a.numArrays_ + b.numArrays_;
    int numGenes = 0;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = a.numArrayHeader_;
    GeneNameScheme nameScheme = a.nameScheme_;

    if (numGeneHeader < b.numGeneHeader_) {
      numGeneHeader = b.numGeneHeader_;
    }

    HashMap<String, LinkedList<Integer> > b_ids = new HashMap<String, LinkedList<Integer> >();
    for (int i = 0; i < b.numGenes_; i++) {
      GeneData gene = b.data_[i+b.numGeneHeader_];
      String id = (String) gene.getDataAt(b.id_);
      if (!b_ids.containsKey(id)) {
        b_ids.put(id, new LinkedList<Integer>());
      }
      LinkedList<Integer> val = b_ids.get(id);
      val.add(new Integer(i+b.numGeneHeader_));
    }

    Vector<GeneData> dataV = new Vector<GeneData>();

    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene_a = null;
      GeneData gene_b = null;
      if ( i < a.numGeneHeader_) {
        gene_a = a.data_[i];
      }
      else {
        Object[] obj = new Object[a.numArrays_+a.numArrayHeader_];
        gene_a = new GeneData(obj);
      }
      if ( i < b.numGeneHeader_) {
        int end = b.numArrays_ + b.numArrayHeader_ - 1;
        gene_b = b.data_[i].subset(b.numArrayHeader_, end);
      }
      else {
        Object[] obj = new Object[b.numArrays_];
        gene_b = new GeneData(obj);
      }
      GeneData gene = gene_a;
      if (!select) {
        gene = GeneData.merge(gene_a, gene_b);
      }
      dataV.add(gene);
    }

    // Genes in a
    String lastId = null;
    for (int i = 0; i < a.numGenes_; i++) {
      GeneData gene_a = a.data_[i+a.numGeneHeader_];
      GeneData gene_b = null;
      String id = (String) gene_a.getDataAt(a.id_);
      int start = b.numArrayHeader_;
      int end = b.numArrays_ + b.numArrayHeader_ - 1;
      if (b_ids.containsKey(id)) {
        // Check the list for lastID
        if (lastId != null && !lastId.equals(id) && b_ids.containsKey(lastId)) {
          LinkedList<Integer> val = b_ids.get(lastId);
          Iterator<Integer> itr = val.iterator();
          while (itr.hasNext()) {
            Integer e = (Integer) itr.next();
            Object[] obj = new Object[a.numArrays_];
            GeneData g_a2 = new GeneData(obj);
            GeneData g_a1 = a.data_[i-1+a.numGeneHeader_].subset(0,a.numArrayHeader_-1);
            GeneData g_a = GeneData.merge(g_a1, g_a2);
            GeneData g_b = b.data_[e.intValue()].subset(start, end);
            GeneData g = g_a;
            if (!select) {
              g = GeneData.merge(g_a, g_b);
              dataV.add(g);
            }
          }
        }
        LinkedList<Integer> val = b_ids.get(id);
        if (val.size() > 0) {
          Integer first = val.removeFirst();
          gene_b = b.data_[first.intValue()].subset(start, end);
        }
      }
      if (!intersect && gene_b == null) {
        Object[] obj = new Object[b.numArrays_];
        gene_b = new GeneData(obj);
      }
      if (gene_b != null) {
        GeneData gene = gene_a;
        if (!select) {
          gene = GeneData.merge(gene_a, gene_b);
        }
        dataV.add(gene);
      }
      lastId = id;
    }

    if (!intersect) {
      // Other Genes in b
      int start = b.numArrayHeader_;
      int end = b.numArrays_ + b.numArrayHeader_ - 1;
      Iterator<String> itr = b_ids.keySet().iterator();
      while (itr.hasNext()) {
        String id = (String) itr.next();
        LinkedList<Integer> val = b_ids.get(id);
        Iterator<Integer> itr1 = val.iterator();
        while (itr1.hasNext()) {
          Integer e = (Integer) itr1.next();
          Object[] obj = new Object[a.numArrays_];
          GeneData g_a2 = new GeneData(obj);
          GeneData g_a1 = b.data_[e.intValue()].subset(0,a.numArrayHeader_-1);
          GeneData g_a = GeneData.merge(g_a1, g_a2);
          GeneData g_b = b.data_[e.intValue()].subset(start, end);
          GeneData g = g_a;
          if (!select) {
            g = GeneData.merge(g_a, g_b);
            dataV.add(g);
          }
        }
      }
    }

    GeneData[] data = new GeneData[dataV.size()];
    numGenes = 0;
    Iterator<GeneData> itr = dataV.iterator();
    while (itr.hasNext()) {
      GeneData gene = (GeneData) itr.next();
      data[numGenes++] = gene;
    }
    numGenes = numGenes - numGeneHeader;

    if (select) {
      numArrays = a.numArrays_;
    }

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  /**
   *  @param a - First Data
   *  @param b - Second Data
   *  @return - new Data with First Data id and not with Second Data id
   */
  public static Data diffData(Data a, Data b)
  {
    int numArrays = a.numArrays_;
    int numGenes = 0;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = a.numArrayHeader_;
    GeneNameScheme nameScheme = a.nameScheme_;

    HashMap<String, Integer> b_ids = new HashMap<String,Integer>();
    for (int i = 0; i < b.numGenes_; i++) {
      GeneData gene = b.data_[i+b.numGeneHeader_];
      String id = (String) gene.getDataAt(b.id_);
      b_ids.put(id, new Integer(i+b.numGeneHeader_));
    }

    Vector<GeneData> dataV = new Vector<GeneData>();

    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene = a.data_[i];
      dataV.add(gene);
    }

    // Genes in a
    for (int i = 0; i < a.numGenes_; i++) {
      GeneData gene_a = a.data_[i+a.numGeneHeader_];
      String id = (String) gene_a.getDataAt(a.id_);
      if (!b_ids.containsKey(id)) {
        dataV.add(gene_a);
      }
    }

    GeneData[] data = new GeneData[dataV.size()];
    numGenes = 0;
    Iterator<GeneData> itr = dataV.iterator();
    while (itr.hasNext()) {
      GeneData gene = (GeneData) itr.next();
      data[numGenes++] = gene;
    }
    numGenes = numGenes - numGeneHeader;

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  /**
   *  @param a - First Data
   *  @param b - Second Data
   *  @return - new Data with First Data selected with Second Data id and
   *            Second Data order
   */
  public static Data selectOrder(Data a, Data b)
  {
    int numArrays = a.numArrays_;
    int numGenes = 0;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = a.numArrayHeader_;
    GeneNameScheme nameScheme = a.nameScheme_;

    HashMap<String, LinkedList<Integer> > a_ids = new HashMap<String, LinkedList<Integer> >();
    for (int i = 0; i < a.numGenes_; i++) {
      GeneData gene = a.data_[i+a.numGeneHeader_];
      String id = (String) gene.getDataAt(a.id_);
      if (!a_ids.containsKey(id)) {
        a_ids.put(id, new LinkedList<Integer>());
      }
      LinkedList<Integer> val = a_ids.get(id);
      val.add(new Integer(i+a.numGeneHeader_));
    }

    Vector<GeneData> dataV = new Vector<GeneData>();

    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene = a.data_[i];
      dataV.add(gene);
    }

    // Genes in a
    for (int i = 0; i < b.numGenes_; i++) {
      GeneData gene_b = b.data_[i+b.numGeneHeader_];
      String id = (String) gene_b.getDataAt(b.id_);
      if (a_ids.containsKey(id)) {
        LinkedList<Integer> val = a_ids.get(id);
        Iterator<Integer> itr = val.iterator();
        while (itr.hasNext()) {
          Integer e = (Integer) itr.next();
          GeneData gene = a.data_[e.intValue()];
          dataV.add(gene);
        }
      }
    }

    GeneData[] data = new GeneData[dataV.size()];
    numGenes = 0;
    Iterator<GeneData> itr = dataV.iterator();
    while (itr.hasNext()) {
      GeneData gene = (GeneData) itr.next();
      data[numGenes++] = gene;
    }
    numGenes = numGenes - numGeneHeader;

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  /**
   *  @param a - First Data
   *  @param b - Second Data
   *  @return - new data with first data gene names seleted by second data id
   */
  public static Data selectNames(Data a, Data b)
  {
    int numArrays = a.numArrays_;
    int numGenes = 0;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = a.numArrayHeader_;
    GeneNameScheme nameScheme = a.nameScheme_;

    HashMap<String, Integer> b_ids = new HashMap<String,Integer>();
    for (int i = 0; i < b.numGenes_; i++) {
      GeneData gene = b.data_[i+b.numGeneHeader_];
      String id = (String) gene.getDataAt(b.id_);
      //System.out.println(":"+id+":");
      b_ids.put(id, new Integer(i+b.numGeneHeader_));
    }

    Vector<GeneData> dataV = new Vector<GeneData>();

    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene = a.data_[i];
      dataV.add(gene);
    }

    // Genes in a
    for (int i = 0; i < a.numGenes_; i++) {
      String geneName = a.getGenesAt(i+a.numGeneHeader_);
      //System.out.println(":"+geneName+":");
      GeneData gene_a = a.data_[i+a.numGeneHeader_];
      if (b_ids.containsKey(geneName)) {
        dataV.add(gene_a);
      }
    }

    GeneData[] data = new GeneData[dataV.size()];
    numGenes = 0;
    Iterator<GeneData> itr = dataV.iterator();
    while (itr.hasNext()) {
      GeneData gene = (GeneData) itr.next();
      data[numGenes++] = gene;
    }
    numGenes = numGenes - numGeneHeader;

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  /**
   *  @param a - First Data
   *  @param b - Second Data
   *  @return - new data with gene names present in the second data (First column)
   */
  public static Data containNames(Data a, Data b)
  {
    int numArrays = a.numArrays_;
    int numGenes = 0;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = a.numArrayHeader_;
    GeneNameScheme nameScheme = a.nameScheme_;

    HashMap<String, Integer> b_ids = new HashMap<String,Integer>();
    for (int i = 0; i < b.numGenes_; i++) {
      GeneData gene = b.data_[i+b.numGeneHeader_];
      String id = (String) gene.getDataAt(0);
      id = id.toUpperCase();
      b_ids.put(id, new Integer(i+b.numGeneHeader_));
    }

    Vector<GeneData> dataV = new Vector<GeneData>();

    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene = a.data_[i];
      dataV.add(gene);
    }

    // Genes in a
    for (int i = 0; i < a.numGenes_; i++) {
      String geneName = a.getGenesAt(i+a.numGeneHeader_);
      if (geneName != null) {
        geneName = geneName.toUpperCase();
        GeneData gene_a = a.data_[i+a.numGeneHeader_];
        if (b_ids.containsKey(geneName)) {
          dataV.add(gene_a);
        }
      }
    }

    GeneData[] data = new GeneData[dataV.size()];
    numGenes = 0;
    Iterator<GeneData> itr = dataV.iterator();
    while (itr.hasNext()) {
      GeneData gene = (GeneData) itr.next();
      data[numGenes++] = gene;
    }
    numGenes = numGenes - numGeneHeader;

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  /**
   *  @param a - First Data
   *  @param b - Second Data
   *  @return - new Data with First Data concat with second data
   */
  public static Data concatData(Data a, Data b)
  {
    int numArrays = a.numArrays_;
    int numGenes = a.numGenes_ + b.numGenes_;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = a.numArrayHeader_;
    GeneNameScheme nameScheme = a.nameScheme_;

    if (numArrays != b.numArrays_) {
        return null;
    }
    if (numArrayHeader != b.numArrayHeader_) {
        return null;
    }

    GeneData[] data = new GeneData[numGenes+numGeneHeader];
    int index = 0;
    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene = a.data_[i];
      data[index++] = gene;
    }

    // Genes in a
    for (int i = 0; i < a.numGenes_; i++) {
      GeneData gene = a.data_[i+a.numGeneHeader_];
      data[index++] = gene;
    }
    // Genes in b
    for (int i = 0; i < b.numGenes_; i++) {
      GeneData gene = b.data_[i+b.numGeneHeader_];
      data[index++] = gene;
    }

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  /**
   *  @param a - First Data
   *  @param b - Second Data
   *  @return - new Data with First Data concat with second data in columns
   */
  public static Data concatDataColumns(Data a, Data b)
  {
    int numCols = a.getNumColumns() + b.getNumColumns();
    int numRows = a.getNumRows();
    GeneNameScheme nameScheme = a.nameScheme_;

    if (numRows < b.getNumRows()) {
        numRows = b.getNumRows();
    }

    GeneData g_a_n = new GeneData(new Object[a.getNumColumns()]);
    GeneData g_b_n = new GeneData(new Object[b.getNumColumns()]);
    GeneData[] data = new GeneData[numRows];
    for (int i =0; i < numRows; i++) {
      GeneData g_a = null, g_b = null;
      if ( i < a.getNumRows() ) {
        g_a = a.data_[i];
      }
      else {
        g_a = g_a_n;
      }
      if ( i < b.getNumRows() ) {
        g_b = b.data_[i];
      }
      else {
        g_b = g_b_n;
      }
      data[i] = GeneData.merge(g_a, g_b);
    }
    int numArrays = numCols - a.numArrayHeader_;
    int numGenes = numRows - a.numGeneHeader_;

    Data res = new Data(numArrays,numGenes,a.numGeneHeader_,a.numArrayHeader_,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  /**
   *  @param a - First Data
   *  @param selectOrder - selection of genes in order
   *  @return - new Data with selected genes in order
   */
  public static Data selectGenesFromData(Data a, int[] selectOrder)
  {
    int numArrays = a.numArrays_;
    int numGenes = selectOrder.length;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = a.numArrayHeader_;
    GeneNameScheme nameScheme = a.nameScheme_;

    if (numGenes > a.numGenes_) {
        return null;
    }

    GeneData[] data = new GeneData[numGenes+numGeneHeader];
    int index = 0;
    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene = a.data_[i];
      data[index++] = gene;
    }

    // Genes in a
    for (int i = 0; i < selectOrder.length; i++) {
      int num = selectOrder[i];
      if (num < a.numGenes_) {
        GeneData gene = a.data_[num+a.numGeneHeader_];
        data[index++] = gene;
      }
      else {
        System.err.println(" **error : bad selection of genes");
      }
    }
    numGenes = index - numGeneHeader;

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);
    res.setTimepoints(a.getTimepoints());

    return res;
  }

  /**
   *  @param a - First Data
   *  @param selectOrder - selection of arrays in order
   *  @return - new Data with selected arrays in order
   */
  public static Data selectArraysFromData(Data a, int[] selectOrder)
  {
    int[] order = new int[a.numArrayHeader_+selectOrder.length];
    for (int i =0; i < a.numArrayHeader_; i++) {
      order[i] = i;
    }
    for (int i =0; i < selectOrder.length; i++) {
      order[i+a.numArrayHeader_] = selectOrder[i] + a.numArrayHeader_;
    }
    Data res = selectColumnsFromData(a, order);
    res.numArrayHeader_ = a.numArrayHeader_;
    res.numArrays_ -= a.numArrayHeader_;
    return res;
  }

  /**
   *  @param a - First Data
   *  @param selectOrder - selection of columns in order
   *  @return - new Data with selected columns in order
   */
  public static Data selectColumnsFromData(Data a, int[] selectOrder)
  {
    int numArrays = selectOrder.length;
    int numGenes = a.numGenes_;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = 0;
    GeneNameScheme nameScheme = a.nameScheme_;

    if (numGenes > a.numGenes_) {
        return null;
    }

    int[] order = selectOrder;

    GeneData[] data = new GeneData[numGenes+numGeneHeader];
    int index = 0;
    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene = a.data_[i];
      data[index++] = gene.permute(0, order);
    }

    // Genes in a
    for (int i = 0; i < a.numGenes_; i++) {
      GeneData gene = a.data_[i+a.numGeneHeader_];
      data[index++] = gene.permute(0, order);
    }

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  /**
   *  @param a - First Data
   *  @param column - select column
   *  @return - new Data with selected column
   */
  public static Data selectColumnFromData(Data a, int column)
  {
    int numArrays = 1;
    int numGenes = a.numGenes_;
    int numGeneHeader = a.numGeneHeader_;
    int numArrayHeader = 0;
    GeneNameScheme nameScheme = a.nameScheme_;

    int[] order = new int[1];
    order[0] = column;

    GeneData[] data = new GeneData[numGenes+numGeneHeader];
    int index = 0;
    // Header Data
    for (int i =0; i < numGeneHeader; i++) {
      GeneData gene = a.data_[i];
      data[index++] = gene.permute(0, order);
    }

    // Genes in a
    for (int i = 0; i < a.numGenes_; i++) {
      GeneData gene = a.data_[i+a.numGeneHeader_];
      data[index++] = gene.permute(0, order);
    }

    Data res = new Data(numArrays,numGenes,numGeneHeader,numArrayHeader,data);
    res.setGeneNameScheme(nameScheme);

    return res;
  }

  public double average() {
    double res = 0.0;
    int count = 0;
    for (int i = numGeneHeader_; i < numGenes_+numGeneHeader_; i++) {
      Double[] values = data_[i].getVector(numArrayHeader_, getNumColumns()-1);
      for (int j=0; j < values.length; j++) {
        if (values[j] != null) {
          res += values[j].doubleValue();
          count++;
        }
      }
    }
    return res/count;
  }

  public void meanCenter() {
    int start = numArrayHeader_;
    int end = getNumColumns()-1;
    for (int i = numGeneHeader_; i < numGenes_+numGeneHeader_; i++) {
      try {
        GeneData gene = data_[i];
        double center = GeneData.getMean(gene.getData(), start, end);
        gene.performCentering(start, center);
      }
      catch(Exception e) {
        // No Action
      }
    } // end for
  }

  public void zeroCenter() {
    int start = numArrayHeader_;
    int end = getNumColumns()-1;
    for (int i = numGeneHeader_; i < numGenes_+numGeneHeader_; i++) {
      try {
        GeneData gene = data_[i];
        Double center = (Double) gene.getDataAt(start);
        if (center != null) {
          gene.performCentering(start, center.doubleValue());
        }
      }
      catch(Exception e) {
        // No Action
      }
    } // end for
  }

  public void normalize() {
    double avg = average();
    System.out.println("Average = " + avg);
    int start = numArrayHeader_;
    for (int i = numGeneHeader_; i < numGenes_+numGeneHeader_; i++) {
      data_[i].performCentering(start, avg);
    } // end for
  }

  public void addConstant(double value) {
    int start = numArrayHeader_;
    int end = getNumColumns()-1;
    for (int i = numGeneHeader_; i < numGenes_+numGeneHeader_; i++) {
      data_[i].addConstant(start, end, value);
    } // end for
  }

  public void addWeightRow(int index) {
    int start = numArrayHeader_;
    int end = getNumColumns()-1;
    GeneData[] data = new GeneData[getNumRows()+1];
    for (int i = 0; i < index; i++) {
      data[i] = data_[i];
    } // end for
    data[index] = GeneData.getWeight(0, end);
    numGenes_ = numGenes_+numGeneHeader_+1 - index - 1;
    numGeneHeader_ = index + 1;
    for (int i = numGeneHeader_; i < numGenes_+numGeneHeader_; i++) {
      data[i] = data_[i-1];
    } // end for
    data_ = data;
  }

  public void addWeightCol(int index) {
    numArrays_ = numArrays_+numArrayHeader_+1 - index - 1;
    numArrayHeader_ = index + 1;
    for (int i = 0; i < numGenes_+numGeneHeader_; i++) {
      data_[i].addWeightCol(index);
    } // end for
  }

};

