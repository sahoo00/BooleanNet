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

package tools.io;

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.lang.Math;

import tools.microarray.GeneData;
import tools.algorithm.Bimodal;

enum  Operator {PLUS, MINUS, TIMES, DIVIDE, LT, GT, LEQ, GEQ, 
    AND, OR, IMPLIES, NOT,
    PHENOTYPE, SEARCH, CSEARCH, RSEARCH, MIN, MAX};

public class NetworkInfo {
  String pcl_file, bv_file, thr_file, idx_file;
  String ph_file, gsminfodb_file, gsminfo_file;
  String rl_file, bvidx_file;

  public NetworkInfo() {
    reset();
  }

  public static boolean exists(String filename) {
    File f = new File(filename);
    return f.exists();
  }

  public void setPclFile(String f) { pcl_file = f; }
  public void setBvFile(String f) { bv_file = f; }
  public void setThrFile(String f) { thr_file = f; }
  public void setIdxFile(String f) { idx_file = f; }
  public void setPhFile(String f) { ph_file = f; }
  public void setGsminfodbFile(String f) { gsminfodb_file = f; }
  public void setGsminfoFile(String f) { gsminfo_file = f; }
  public void setRlFile(String f) { rl_file = f; }
  public void setBvidxFile(String f) { bvidx_file = f; }

  public void setAllFiles(String prefix) {
    if (exists(prefix + ".pcl")) {
      pcl_file = prefix + ".pcl";
    }
    if (exists(prefix + ".bv")) {
      bv_file = prefix + ".bv";
    }
    if (exists(prefix + ".thr")) {
      thr_file = prefix + ".thr";
    }
    if (exists(prefix + ".idx")) {
      idx_file = prefix + ".idx";
    }
    if (exists(prefix + ".ph")) {
      ph_file = prefix + ".ph";
    }
    if (exists(prefix + ".gsminfodb")) {
      gsminfodb_file = prefix + ".gsminfodb";
    }
    if (exists(prefix + ".gsminfo")) {
      gsminfo_file = prefix + ".gsminfo";
    }
    if (exists(prefix + ".rl")) {
      rl_file = prefix + ".rl";
    }
    if (exists(prefix + ".bvidx")) {
      bvidx_file = prefix + ".bvidx";
    }
  }

  public void printAllFiles() {
    System.out.println(pcl_file);
    System.out.println(bv_file);
    System.out.println(thr_file);
    System.out.println(idx_file);
    System.out.println(ph_file);
    System.out.println(gsminfodb_file);
    System.out.println(gsminfo_file);
    System.out.println(rl_file);
    System.out.println(bvidx_file);
  }

  public void reset() {
    pcl_file = bv_file = thr_file = idx_file = null;
    ph_file = gsminfodb_file = gsminfo_file = null;
  }

  HashMap<String, LinkedList<Long> > idhash_;
  HashMap<String, LinkedList<Long> > namehash_;
  HashMap<Long, Long> revptrhash_;
  HashMap<Long, String> revidhash_;
  HashMap<Long, String> revnamehash_;
  HashMap<Long, Vector<Double> > thrhash_;
  Vector<String> header_;
  HashMap<String, Integer> headerHash_;
  BufferedRandomAccessFile pclReader_;
  int numArrays_, numArrayHeader_;

  public void readIndex() throws IOException {
    readFileIndex(idx_file);
  }

  public void readBvIndex() throws IOException {
    readFileIndex(bvidx_file);
  }

  public HashMap<Long, Long> getLineMap() { return revptrhash_; }

  public void readFileIndex(String file) throws IOException {
    BufferedReader reader = getReader(file);
    String record = reader.readLine(); // Header

    idhash_ = new HashMap<String, LinkedList<Long> >();
    namehash_ = new  HashMap<String, LinkedList<Long> >();
    revptrhash_ = new HashMap<Long, Long>();
    revidhash_ = new HashMap<Long, String>();
    revnamehash_ = new HashMap<Long, String>();

    long index = 0;
    while (record != null) {
      record = reader.readLine();
      if (record == null) {
        break;
      }
      record.trim();
      String[] result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
      String idu = result[0].toUpperCase();
      String nameu = result[1].toUpperCase();
      Long ptr = new Long(Long.parseLong(result[2]));
      String desc = result[3];
      Long idx = new Long(index);
      if (!idhash_.containsKey(idu)) {
        LinkedList<Long> list = new LinkedList<Long>();
        idhash_.put(idu, list);
      }
      if (!namehash_.containsKey(nameu)) {
        LinkedList<Long> list = new LinkedList<Long>();
        namehash_.put(nameu, list);
      }
      idhash_.get(idu).add(idx);
      namehash_.get(nameu).add(idx);
      revptrhash_.put(idx, ptr);
      revidhash_.put(idx, result[0]);
      revnamehash_.put(idx, result[1]);
      index++;
    }
    reader.close();
  }

  public void readThresholds() throws Exception {
    BufferedReader reader = getReader(thr_file);

    thrhash_ = new HashMap<Long, Vector<Double> >();

    String record = reader.readLine(); // Header
    while (record != null) {
      record.trim();
      String[] result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
      Long index = new Long(Long.parseLong(result[0]));
      Vector<Double> list = new Vector<Double>();
      list.add(new Double(Double.parseDouble(result[1]))); // thr1
      list.add(new Double(Double.parseDouble(result[2]))); // stat
      list.add(new Double(Double.parseDouble(result[3]))); // thr0
      list.add(new Double(Double.parseDouble(result[4]))); // thr2
      thrhash_.put(index, list);

      record = reader.readLine();
    }
    reader.close();
  }

  public void readPCLHeader() throws Exception {
    BufferedReader reader = getReader(pcl_file);
    String record = reader.readLine(); // Header
    record.trim();
    String[] result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
    numArrayHeader_ = 3;
    numArrays_ = result.length - numArrayHeader_;
    header_ = new Vector<String>(Arrays.asList(result));
    reader.close();
    buildHeaderHash();
  }

  public void buildHeaderHash() throws Exception {
    headerHash_ = new HashMap<String, Integer>();
    Pattern p1 = Pattern.compile("[Gg][Ss][Mm]([0-9]+)");
    Pattern p2 = Pattern.compile("(FL_[0-9]+)_");
    for (int i =0; i < header_.size(); i++) {
      String str = header_.get(i);
      String key = str;
      Matcher m = p1.matcher(str);
      if (m.matches()) {
        key = m.group(1);
      }
      else {
        m = p2.matcher(str);
        if (m.matches()) {
          key = m.group(1);
        }
      }
      //System.out.println(key + " " + i);
      headerHash_.put(key, new Integer(i));
    }
  }

  public void printHeader() {
    for (int i = 0; i < header_.size(); i++) {
      if (i != 0) System.out.print("\t");
      System.out.print(header_.get(i));
    }
    System.out.println();
    int nCols = numCols();
    System.out.print("EWEIGHT\t\t\t");
    for (int i = numArrayHeader_; i < nCols; i++) {
      System.out.print("\t1");
    }
    System.out.println();
  }

  public void initPCLReader() throws IOException {
    pclReader_ = new BufferedRandomAccessFile(pcl_file, "r");
  }

  public void seekPCL(long ptr) throws IOException {
    if (ptr != pclReader_.getFilePointer()) {
      pclReader_.seek(ptr);
    }
  }

  public String getLineAt(long ptr) throws IOException {
    seekPCL(ptr);
    return pclReader_.getNextLine();
  }

  public int numCols() {
    return numArrays_ + numArrayHeader_;
  }

  public GeneData getGeneDataAt(long ptr) throws IOException {
    String record = getLineAt(ptr);
    String[] result = record.split("\\t", -2);
    int nCols = numArrays_ + numArrayHeader_;
    if (nCols != result.length) {
      System.out.println(" *Warning* Column mismatch - Orig :" + nCols + 
          ", New : " + result.length);
    }
    // System.out.println(record);
    Object[] d = new Object[nCols];
    for (int j = 0; j < Math.min(result.length, nCols); j++) {
      d[j] = result[j];
    }
    GeneData res = new GeneData(d);
    return res;
  }

  public void closePCLReader() throws IOException {
    pclReader_.close();
  }

  public Long getIndexByName(String id) throws IOException {
    if (idhash_.containsKey(id.toUpperCase())) {
      Long index = idhash_.get(id.toUpperCase()).getFirst();
      return index;
    }
    if (namehash_.containsKey(id.toUpperCase())) {
      Long index = namehash_.get(id.toUpperCase()).getFirst();
      return index;
    }
    return null;
  }

  public long getFilePointer(Long index) {
    return revptrhash_.get(index).longValue();
  }

  public Double getThresholdByIndex(Long index) {
    return thrhash_.get(index).get(0);
  }

  public Double getUpperThresholdByIndex(Long index) {
    return thrhash_.get(index).get(3);
  }

  public Double getLowerThresholdByIndex(Long index) {
    return thrhash_.get(index).get(2);
  }

  public GeneData getGeneDataByName(String id) throws IOException {
    Long index = getIndexByName(id);
    if (index != null) {
      long ptr = getFilePointer(index);
      return getGeneDataAt(ptr);
    }
    return null;
  }

  public void init() throws Exception {
    readIndex();
    readThresholds();
    readPCLHeader();
    initPCLReader();
  }

  public void close() throws IOException {
    closePCLReader();
  }

  public static Boolean getBoolean(Object o) throws Exception {
    if (o == null) {
        return null;
    }
    if (o instanceof Boolean) {
      return (Boolean) o;
    }
    if (o instanceof Double) {
        Double res = (Double) o;
        return new Boolean(res.doubleValue() > 0);
    }
    if (o instanceof String) {
        try {
            return new Boolean(Double.parseDouble((String)o) > 0);
        }
        catch (Exception e) {
            return null;
        }
    }
    return null;
  }

  public static Double getDouble(Object o) throws Exception {
    if (o == null) {
        return null;
    }
    if (o instanceof String) {
        try {
            return new Double(Double.parseDouble((String)o));
        }
        catch (Exception e) {
            return null;
        }
    }
    if (o instanceof Double) {
        return (Double) o;
    }
    if (o instanceof Boolean) {
        Boolean res = (Boolean) o;
        if (res.booleanValue()) {
            return new Double(1.0);
        }
        else {
            return new Double(0.0);
        }
    }
    return null;
  }

  public GeneData operateDouble(Operator type, GeneData q1, GeneData q2) throws Exception {
    int nCols = numCols();
    Object[] d = new Object[nCols];
    for (int j = numArrayHeader_; j < nCols; j++) {
      Double v1 = getDouble(q1.getDataAt(j));
      Double v2 = getDouble(q2.getDataAt(j));
      if (v1 == null) continue;
      if (v2 == null) continue;
      if (type == Operator.PLUS) {
        d[j] = new Double(v1.doubleValue() + v2.doubleValue());
      }
      if (type == Operator.MINUS) {
        d[j] = new Double(v1.doubleValue() - v2.doubleValue());
      }
      if (type == Operator.TIMES) {
        d[j] = new Double(v1.doubleValue() * v2.doubleValue());
      }
      if (type == Operator.DIVIDE) {
        if (v2.doubleValue() == 0) continue;
        d[j] = new Double(v1.doubleValue() / v2.doubleValue());
      }
      if (type == Operator.LT) {
        d[j] = new Boolean(v1.doubleValue() < v2.doubleValue());
      }
      if (type == Operator.GT) {
        d[j] = new Boolean(v1.doubleValue() > v2.doubleValue());
      }
      if (type == Operator.LEQ) {
        d[j] = new Boolean(v1.doubleValue() <= v2.doubleValue());
      }
      if (type == Operator.GEQ) {
        d[j] = new Boolean(v1.doubleValue() >= v2.doubleValue());
      }
      if (type == Operator.MIN) {
        d[j] = new Double(Math.min(v1.doubleValue(), v2.doubleValue()));
      }
      if (type == Operator.MAX) {
        d[j] = new Double(Math.max(v1.doubleValue(), v2.doubleValue()));
      }
    }
    GeneData res = new GeneData(d);
    return res;
  }

  public GeneData operateBoolean(Operator type, GeneData q1, GeneData q2) throws Exception {
    int nCols = numCols();
    Object[] d = new Object[nCols];
    for (int j = numArrayHeader_; j < nCols; j++) {
      Boolean v1 = getBoolean(q1.getDataAt(j));
      Boolean v2 = getBoolean(q2.getDataAt(j));
      if (v1 == null) continue;
      if (v2 == null) continue;
      if (type == Operator.AND) {
        d[j] = new Boolean(v1.booleanValue() && v2.booleanValue());
      }
      if (type == Operator.OR) {
        d[j] = new Boolean(v1.booleanValue() || v2.booleanValue());
      }
      if (type == Operator.IMPLIES) {
        d[j] = new Boolean(!v1.booleanValue() || v2.booleanValue());
      }
      if (type == Operator.NOT) {
        d[j] = new Boolean(!v1.booleanValue());
      }
    }
    GeneData res = new GeneData(d);
    return res;
  }

  public GeneData plus(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.PLUS, q1, q2);
  }

  public GeneData minus(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.MINUS, q1, q2);
  }

  public GeneData times(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.TIMES, q1, q2);
  }

  public GeneData divide(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.DIVIDE, q1, q2);
  }

  public GeneData searchPhenotype(Operator type, String s) throws Exception {
    FileReader fr = new FileReader(gsminfo_file);
    BufferedReader reader = new BufferedReader(fr);
    String record = reader.readLine();
    Pattern p1 = Pattern.compile("^\\^SAMPLE = [Gg][Ss][Mm]([0-9]+)");
    Pattern p2 = Pattern.compile("^\\^SAMPLE = (.+)\\s*");
    Pattern ps = Pattern.compile(s);
    String key = null;
    int nCols = numCols();
    Object[] d = new Object[nCols];
    for (int i =0; i < nCols; i++) {
      d[i] = new Boolean(false);
    }
    while (record != null) {
      Matcher m = p1.matcher(record);
      if (m.matches()) {
        key = m.group(1);
      }
      else {
        m = p2.matcher(record);
        if (m.matches()) {
          key = m.group(1);
        }
      }
      if (headerHash_.containsKey(key)) {
        if (type == Operator.SEARCH) { // search : case insensitive
          if (record.toUpperCase().indexOf(s.toUpperCase()) != -1) {
            Integer loc = headerHash_.get(key);
            d[loc.intValue()] = new Boolean(true);
            key = null;
          }
        }
        if (type == Operator.CSEARCH) { // csearch : case sensitive
          if (record.indexOf(s) != -1) {
            Integer loc = headerHash_.get(key);
            d[loc.intValue()] = new Boolean(true);
            key = null;
          }
        }
        if (type == Operator.RSEARCH) { // rsearch : case sensitive
          m = ps.matcher(record);
          if (m.matches()) {
            Integer loc = headerHash_.get(key);
            d[loc.intValue()] = new Boolean(true);
            key = null;
          }
        }
      }
      record = reader.readLine();
    }
    reader.close();
    fr.close();
    return new GeneData(d);
  }

  public GeneData getPhenotype(String s) throws Exception {
    BufferedReader reader = getReader(ph_file);
    int nCols = numCols();
    Object[] d = new Object[nCols];
    for (int i =0; i < nCols; i++) {
      d[i] = new Boolean(false);
    }
    while (true) {
      String record = reader.readLine();
      if (record == null) {
        break;
      }
      record.trim();
      String[] result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
      String id = result[0].trim();
      if (id.equals(s)) {
        if (nCols != result.length) {
          System.out.println(" *Warning* ph Column mismatch - Orig :" + nCols + 
              ", New : " + result.length);
        }
        for (int i =0; i < nCols; i++) {
          d[i] = result[i];
        }
        GeneData q = new GeneData(d);
        int start = numArrayHeader_;
        int end = numCols()-1;
        q.convertDouble(start, end);
        return q;
      }
    }
    reader.close();
    return null;
  }

  public GeneData searchPhenotype(String s) throws Exception {
    return searchPhenotype(Operator.SEARCH, s);
  }

  public GeneData csearchPhenotype(String s) throws Exception {
    return searchPhenotype(Operator.CSEARCH, s);
  }

  public GeneData rsearchPhenotype(String s) throws Exception {
    return searchPhenotype(Operator.RSEARCH, s);
  }

  public static BufferedReader getReader(String filename) throws IOException {
    BufferedReader reader = null;
    if (filename.startsWith("http:")) {
      URL url = new URL(filename);
      if (filename.endsWith(".gz")) {
        reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(url.openStream())));
      }
      else {
        reader = new BufferedReader(new InputStreamReader(url.openStream()));
      }
    }
    else {
      if (filename.endsWith(".gz")) {
        reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename))));
      }
      else {
        FileReader fr = new FileReader(filename);
        reader = new BufferedReader(fr);
      }
    }
    return reader;
  }

  public GeneData getNumber(Double n) {
    Object[] data = new Object[numCols()];
    for (int i =0; i < numCols(); i++) {
      data[i] = n;
    }
    return new GeneData(data);
  }

  public GeneData getThreshold(GeneData q) throws Exception {
    String id = (String) q.getDataAt(0);
    if (id != null) {
      Long index = getIndexByName(id);
      if (index != null) {
        //System.out.println("Found thr");
        Double thr = getThresholdByIndex(index);
        return getNumber(thr);
      }
    }
    int start = numArrayHeader_;
    int end = numCols()-1;
    q.convertDouble(start, end);
    Double[] v= q.getVector(start, end);
    Bimodal b = new Bimodal(v);
    double thr = b.getThreshold();
    return getNumber(new Double(thr));
  }

  public GeneData getUpperThreshold(GeneData q) throws Exception {
    String id = (String) q.getDataAt(0);
    if (id != null) {
      Long index = getIndexByName(id);
      if (index != null) {
        //System.out.println("Found thr");
        Double thr = getUpperThresholdByIndex(index);
        return getNumber(thr);
      }
    }
    int start = numArrayHeader_;
    int end = numCols()-1;
    q.convertDouble(start, end);
    Double[] v= q.getVector(start, end);
    Bimodal b = new Bimodal(v);
    double thr = b.getHighThreshold();
    return getNumber(new Double(thr));
  }

  public GeneData getLowerThreshold(GeneData q) throws Exception {
    String id = (String) q.getDataAt(0);
    if (id != null) {
      Long index = getIndexByName(id);
      if (index != null) {
        //System.out.println("Found thr");
        Double thr = getLowerThresholdByIndex(index);
        return getNumber(thr);
      }
    }
    int start = numArrayHeader_;
    int end = numCols()-1;
    q.convertDouble(start, end);
    Double[] v= q.getVector(start, end);
    Bimodal b = new Bimodal(v);
    double thr = b.getLowThreshold();
    return getNumber(new Double(thr));
  }

  public GeneData getHigh(GeneData q) throws Exception {
    GeneData thr = getUpperThreshold(q);
    return greaterThanEq(q, thr);
  }

  public GeneData getLow(GeneData q) throws Exception {
    GeneData thr = getLowerThreshold(q);
    return lessThanEq(q, thr);
  }

  public GeneData getMed(GeneData q) throws Exception {
    GeneData high = getHigh(q);
    GeneData low = getHigh(q);
    return and(not(high), not(low));
  }

  public GeneData equal(GeneData q1, GeneData q2) throws Exception {
    int nCols = numCols();
    Object[] d = new Object[nCols];
    for (int j = numArrayHeader_; j < nCols; j++) {
      Object v1 = q1.getDataAt(j);
      Object v2 = q2.getDataAt(j);
      if (v1 == null) continue;
      if (v2 == null) continue;
      if (v1.equals(v2)) {
        d[j] = new Boolean(true);
      }
      else {
        d[j] = new Boolean(false);
      }
    }
    GeneData res = new GeneData(d);
    return res;
  }

  public GeneData lessThan(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.LT, q1, q2);
  }

  public GeneData greaterThan(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.GT, q1, q2);
  }

  public GeneData lessThanEq(GeneData q1, GeneData q2) throws Exception {
  return operateDouble(Operator.LEQ, q1, q2);
  }

  public GeneData greaterThanEq(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.GEQ, q1, q2);
  }

  public GeneData and(GeneData q1, GeneData q2) throws Exception {
    return operateBoolean(Operator.AND, q1, q2);
  }

  public GeneData or(GeneData q1, GeneData q2) throws Exception {
    return operateBoolean(Operator.OR, q1, q2);
  }

  public GeneData implies(GeneData q1, GeneData q2) throws Exception {
    return operateBoolean(Operator.IMPLIES, q1, q2);
  }

  public GeneData not(GeneData q) throws Exception {
    return operateBoolean(Operator.NOT, q, q);
  }

  public GeneData min(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.MIN, q1, q2);
  }

  public GeneData max(GeneData q1, GeneData q2) throws Exception {
    return operateDouble(Operator.MAX, q1, q2);
  }

}

