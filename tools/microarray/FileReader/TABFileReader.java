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

package tools.microarray.FileReader;

import tools.microarray.Data;
import tools.microarray.GeneData;
import tools.microarray.ArrayException;
import java.io.*;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.util.LinkedList;
import java.util.ListIterator;
import java.util.HashMap;
import java.lang.Math;

public class TABFileReader {

  String filename_;
  BufferedReader reader_;
  String[] header_;
  int numArrays_;
  int numArrayHeader_;

  BufferedRandomAccessFile randomReader_;
  HashMap<Long, Long> lineMap_; // Map Line number to FilePointer

  int state_;
  int lineno_;

  public TABFileReader(String file) {
    filename_ = file;
    header_ = null;
    numArrays_ = 0;
    numArrayHeader_ = 0;
    state_ = State.INIT;
    lineno_ = 0;
  }

  public int getLineNumber() { return lineno_; }
  public int getNumArrays() { return numArrays_; }
  public int getNumArrayHeader() { return numArrayHeader_; }
  public int getNumColumns() { return numArrayHeader_ + numArrays_; }

  public void startReader() throws IOException {
    if (filename_.startsWith("http:")) {
      URL url = new URL(filename_);
      if (filename_.endsWith(".gz")) {
        reader_ = new BufferedReader(new InputStreamReader(new GZIPInputStream(url.openStream())));
      }
      else {
        reader_ = new BufferedReader(new InputStreamReader(url.openStream()));
      }
    }
    else {
      if (filename_.endsWith(".gz")) {
        reader_ = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename_))));
      }
      else {
        FileReader fr = new FileReader(filename_);
        reader_ = new BufferedReader(fr);
      }
    }
  }

  public void begin() throws IOException {
    lineno_ = 0;
    startReader();
    String record = reader_.readLine();
    if (record == null) {
      state_ = State.ERROR;
      throw new IOException("No header - 1");
    }
    lineno_++;
    String[] result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
    numArrays_ = result.length;
    header_ = result;

    // Dill: Debugging code
    System.out.println("Header");
    for (int i=0; i < result.length; i++) {
      System.out.println(i + ": "+ result[i]);
    }
    state_ = State.OPENED;
  }

  public void beginRandomAccess() throws IOException {
    if (filename_.startsWith("http:") || filename_.endsWith(".gz")) {
      throw new IOException("Can't open RandomAccess on file: " + filename_);
    }
    randomReader_ = new BufferedRandomAccessFile(filename_, "r");
    lineno_ = 0;
    System.out.println("Building Indices...");
    lineMap_ = new HashMap<Long, Long>();
    lineMap_.put(new Long(lineno_), new Long(randomReader_.getFilePointer()));
    String record = randomReader_.getNextLine();
    if (record == null) {
      state_ = State.ERROR;
      throw new IOException("No header - 1");
    }
    lineno_++;
    String[] result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
    numArrays_ = result.length;
    header_ = result;

    // Dill: Debugging code
    System.out.println("Header");
    for (int i=0; i < result.length; i++) {
      System.out.println(i + ": "+ result[i]);
    }
    lineMap_.put(new Long(lineno_), new Long(randomReader_.getFilePointer()));
    while ((record = randomReader_.getNextLine()) != null) {
      lineno_++;
      if ( (lineno_ % 1000) == 0 ) {
        System.out.println(lineno_ + "\t" + randomReader_.getFilePointer());
      }
      lineMap_.put(new Long(lineno_), new Long(randomReader_.getFilePointer()));
    }
    randomReader_.seek(0);
    System.out.println("Done");
    state_ = State.OPENED;
  }

  public boolean hasNext() {
    return state_ == State.OPENED;
  }

  public GeneData getHeader() {
    GeneData res = null;
    int nCols = numArrays_ + numArrayHeader_;
    Object[] d = new Object[nCols];
    for (int j = 0; j < Math.min(header_.length, nCols); j++) {
      d[j] = header_[j];
    }
    res = new GeneData(d);
    return res;
  }

  public GeneData getData() throws IOException {
    GeneData res = null;
    String record = reader_.readLine();
    if (record == null) {
      state_ = State.CLOSED;
      return null;
    }
    lineno_++;
    String[] result = record.split("\\t", -2);
    int nCols = numArrays_ + numArrayHeader_;
    if (nCols != result.length) {
      System.out.println(" *Warning* Column mismatch - Orig :" + nCols + 
          ", New : " + result.length + " at line " + lineno_ );
    }
    // System.out.println(record);
    Object[] d = new Object[nCols];
    for (int j = 0; j < Math.min(result.length, nCols); j++) {
      d[j] = result[j];
    }
    res = new GeneData(d);
    return res;
  }

  // Call beginRandomAccess before this function
  //    lineno = -1 -> readCurrentLine
  public GeneData getDataAt(long lineno) throws IOException {
    GeneData res = null;
    if (lineno != -1 && !lineMap_.containsKey(new Long(lineno))) {
      return null;
    }
    if (lineno != -1) {
      Long ptr = lineMap_.get(new Long(lineno));
      randomReader_.seek(ptr.longValue());
    }
    String record = randomReader_.getNextLine();
    if (record == null) {
      state_ = State.CLOSED;
      return null;
    }
    lineno_++;
    String[] result = record.split("\\t", -2);
    int nCols = numArrays_ + numArrayHeader_;
    if (nCols != result.length) {
      System.out.println(" *Warning* Column mismatch - Orig :" + nCols + 
          ", New : " + result.length + " at line " + lineno_ );
    }
    // System.out.println(record);
    Object[] d = new Object[nCols];
    for (int j = 0; j < Math.min(result.length, nCols); j++) {
      d[j] = result[j];
    }
    res = new GeneData(d);
    return res;
  }

  public static Data readFile(String filename) throws Exception {
    System.out.println("Reading file " + filename);
    FileReader     fr;
    BufferedReader br;

    if (filename.startsWith("http:")) {
      URL url = new URL(filename);
      if (filename.endsWith(".gz")) {
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(url.openStream())));
      }
      else {
        br = new BufferedReader(new InputStreamReader(url.openStream()));
      }
    }
    else {
      if (filename.endsWith(".gz")) {
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename))));
      }
      else {
        fr = new FileReader(filename);
        br = new BufferedReader(fr);
      }
    }

    int lineno = 0;
    int numArrays = 0;
    int numGenes = 0;
    int numGeneHeader = 0;
    int numArrayHeader = 0;

    LinkedList<String[] > data = new LinkedList<String[] >();
    String[] result;

    String record = br.readLine();
    if (record != null) {
      result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
      numArrays = result.length;
      data.add(result);
    }

    while ((record = br.readLine()) != null) {
      lineno++;
      result = record.split("\\t", -2);
      if (numArrays != result.length) {
        System.out.println(" *Warning* Column mismatch - Orig :" + numArrays + 
            ", New : " + result.length + " at line " + lineno );
        // System.out.println(record);
        String[] res = new String[numArrays];
        for (int i=0; i < Math.min(result.length, numArrays); i++) {
          res[i] = result[i];
        }
        result = res;
      }
      data.add(result);
    }

    numGenes = data.size();
    GeneData[] data_ = new GeneData[data.size()];
    ListIterator iter = data.listIterator();
    int i = 0;
    while (iter.hasNext()) {
        result = (String[]) iter.next();
        Object[] d = new Object[result.length];
        for (int j = 0; j < result.length; j++) {
            d[j] = result[j];
        }
        data_[i] = new GeneData(d);
        i++;
    }
    Data res = new Data(numArrays, numGenes, numGeneHeader, numArrayHeader, 
                    data_);
    System.out.println("Done");
    return res;

  }

};

