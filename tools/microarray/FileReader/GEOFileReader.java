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

import tools.microarray.*;
import java.io.*;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.util.*;
import java.lang.Math;

public class GEOFileReader {

  public static int START = -1;
  public static int DATABASE = 0;
  public static int SERIES = 1;
  public static int PLATFORM = 2;
  public static int SAMPLE = 3;
  public static int DONE = 4;

  public static int MAX_ARRAYS = 300;
  public static boolean RANDOM = false;

  public static int TITLE_INDEX = 7;
  public static int SYM_INDEX = 8;

  public static Data readFile(String filename) throws Exception {
    HashSet<String> excludeList = new HashSet<String>();
    return readFile(filename, excludeList);
  }

  public static Data readInfo(String filename) throws Exception {
    HashSet<String> excludeList = new HashSet<String>();
    return readInfo(filename, excludeList);
  }

  public static Data readFile(String filename, HashSet<String> excludeList) throws Exception {
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
    int numGeneHeader = 2;
    int numArrayHeader = 3;

    int state = START;
    String name = "None";
    String org = "Hs";
    String sampleName = "None";
    String sampleTitle = "None";
    Vector<String> arrayName = new Vector<String>();
    Vector<String> id = new Vector<String>();
    HashMap<String, Integer> idHash = new HashMap<String, Integer>();
    Vector<String> geneName = new Vector<String>();
    Vector<String[] > data = new Vector<String[] >();
    Vector<String> sampleId = new Vector<String>();
    HashSet<String> sampleIdselect = null;

    String record = br.readLine();
    while (record != null) {
      if (record.startsWith("^DATABASE")) {
        state = DATABASE;
      } 
      if (record.startsWith("^SERIES")) {
        state = SERIES;
        name = record.split("=", -2)[1];
      } 
      if (record.startsWith("^PLATFORM")) {
        state = PLATFORM;
      } 
      if (record.startsWith("^SAMPLE")) {
        state = SAMPLE;
        sampleName = record.split("=", -2)[1];
        sampleName = sampleName.replaceAll("\\s", "");
        sampleTitle = "None";
      } 
      // No sample ID from the SERIES state
      // if (state == SERIES && record.startsWith("!Series_sample_id")) {
      //   numArrays ++;
      //   String n = record.split("=", -2)[1];
      //   arrayName.add(n);
      // }
      if (state == PLATFORM && record.startsWith("!Platform_organism")) {
        org = record.split("=", -2)[1];
      }
      if (state == PLATFORM && record.startsWith("!Platform_sample_id")) {
        String sid = record.split("=", -2)[1];
        sid = sid.replaceAll("\\s", "");
        sampleId.add(sid);
      }
      if (state == PLATFORM && record.startsWith("!platform_table_begin")) {
        // Read platform table
        record =  br.readLine();  // ID ...
        if (record == null) {
          throw new ArrayException("No platform data - 1");
        }
        String[] res = record.split("\\t", -2);
        int titleIndex = TITLE_INDEX;
        int symbolIndex = SYM_INDEX;
        for (int i =0; i < res.length; i++) {
            if (res[i].equals("Gene Symbol")) {
                symbolIndex = i;
            }
            if (res[i].equals("Gene Title")) {
                titleIndex = i;
            }
        }
        int index = 0;
        while (record != null) {
          record =  br.readLine();
          if (record == null) {
            throw new ArrayException("No platform data - 2");
          }
          if (record.startsWith("!platform_table_end")) {
            break;
          }
          res = record.split("\\t", -2);
          id.add(res[0]);
          geneName.add(res[symbolIndex]+": " + res[titleIndex]);
          idHash.put(res[0], new Integer(index));
          index++;
        }
        // Select sample id randomly
        if (RANDOM && sampleId.size() > MAX_ARRAYS) {
          sampleIdselect = new HashSet<String>();
          tools.Permutation permutation = new tools.Permutation(100);
          int[] perm = permutation.getRandomPermutation(sampleId.size());
          for (int i =0; i < MAX_ARRAYS; i++) {
            String s = sampleId.get(perm[i]);
            System.out.println(s);
            sampleIdselect.add(s);
          }
        }
      }
      if (state == SAMPLE && record.startsWith("!Sample_data_row_count")) {
        String ns = record.split("=", -2)[1];
        ns = ns.replaceAll(" ", "");
        int n = Integer.parseInt(ns);
        if (n != id.size()) {
            System.out.println("# genes sample ("+n+") != #genes platform (" + id.size() + ")");
            // throw new ArrayException("# genes sample ("+n+") != #genes platform (" + id.size() + ")");
        }
      }
      if (state == SAMPLE && record.startsWith("!Sample_title")) {
        sampleTitle = record.split("=", -2)[1];
      }
      if (state == SAMPLE && record.startsWith("!sample_table_begin")) {
        if (sampleIdselect == null && MAX_ARRAYS > 0 && numArrays >= MAX_ARRAYS) {
            break;
        }
        if (excludeList.contains(sampleName) ||
            sampleIdselect != null && !sampleIdselect.contains(sampleName)) {
          System.out.println(numArrays + ":" + sampleName + " Excluded");
          while (record != null) {
            record =  br.readLine();
            if (record == null) {
              break;
            }
            if (record.startsWith("!sample_table_end")) {
              break;
            }
          }
        }
        else {
          System.out.println(numArrays + ":" + sampleName + "\t" + sampleTitle);
          numArrays++;
          arrayName.add(sampleName);
          // Read sample table
          record =  br.readLine();  // ID ...
          String[] values = new String[id.size()];
          if (record == null) {
            throw new ArrayException("No sample data - 1");
          }
          int count = 0;
          while (record != null) {
            record =  br.readLine();
            if (record == null) {
              throw new ArrayException("No sample data - 2");
            }
            if (record.startsWith("!sample_table_end")) {
              break;
            }
            String[] res = record.split("\\t", -2);
            Integer index = idHash.get(res[0]);
            if (index == null) {
              throw new ArrayException("sample id " + res[0] + " is not found");
            }
            if (res.length > 2) {
              int call = 2;
              int value = 1;
              call = 2; value = 1;
              //call = 4; value = 3;
              res[call] = res[call].replaceAll("\\s", "");
              //double pvalue = Double.parseDouble(res[call]);
              if (!res[call].equals("A"))
              //if (pvalue < 0.05)
              {
                values[index.intValue()] = res[value];
              }
            }
            else if (res.length > 1) {
              values[index.intValue()] = res[1];
            }
            else {
              System.out.println("Error in line :" + record);
            }
            count++;
          }
          data.add(values);
        }
      }

      record =  br.readLine();
    }
    state = DONE;
    if (numArrays != arrayName.size()) {
      throw new ArrayException("# arrays sample != # arrays platform");
    }

    numArrayHeader = 3;
    numGeneHeader = 2;
    numGenes = id.size();
    GeneData[] data_ = new GeneData[id.size() + 2];
    // header 
    Object[] h1 = new Object[numArrayHeader+numArrays];
    h1[0] = "AID"; h1[1] = "Name"; h1[2] = "GWEIGHT";
    for (int i =0; i < arrayName.size(); i++) {
        h1[i+3] = arrayName.get(i);
    }
    data_[0] = new GeneData(h1);
    h1 = new Object[numArrayHeader+numArrays];
    h1[0] = "EWEIGHT"; h1[1] = ""; h1[2] = "";
    for (int i =0; i < arrayName.size(); i++) {
        h1[i+3] = "1";
    }
    data_[1] = new GeneData(h1);
    for (int j =0; j < numGenes; j++) {
      if ( (j%100) == 0) {
        System.out.println(j);
      }
      h1 = new Object[numArrayHeader+numArrays];
      h1[0] = id.get(j); h1[1] = geneName.get(j); h1[2] = "1";
      for (int i =0; i < arrayName.size(); i++) {
        h1[i+3] = data.get(i)[j];
      }
      data_[j + numGeneHeader] = new GeneData(h1);
    }

    Data res = new Data(numArrays, numGenes, numGeneHeader, numArrayHeader, 
                    data_);
    res.setName(name);
    GeneNameScheme ns = new GeneNameScheme(1, org, ":", 0);
    ns.setNumMissingPoints(0);
    res.setGeneNameScheme(ns);
    System.out.println("Done");
    return res;

  }

  public static Data readInfo(String filename, HashSet<String> excludeList) throws Exception {
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

    int state = START;
    String name = "None";
    String org = "Hs";
    String sampleName = "NA";
    String sampleTitle = "NA";
    String sampleSource = "NA";
    String sampleDescription = "NA";
    Vector<String[] > data = new Vector<String[] >();
    Vector<String> sampleId = new Vector<String>();
    HashSet<String> sampleIdselect = null;

    String record = br.readLine();
    while (record != null) {
      if (record.startsWith("^DATABASE")) {
        state = DATABASE;
      } 
      if (record.startsWith("^SERIES")) {
        state = SERIES;
        name = record.split("=", -2)[1];
      } 
      if (record.startsWith("^PLATFORM")) {
        state = PLATFORM;
      } 
      if (record.startsWith("^SAMPLE")) {
        state = SAMPLE;
        sampleName = record.split("=", -2)[1];
        sampleName = sampleName.replaceAll("\\s", "");
        sampleTitle = "NA";
        sampleSource = "NA";
        sampleDescription = "NA";
      } 
      // No sample ID from the SERIES state
      // if (state == SERIES && record.startsWith("!Series_sample_id")) {
      //   numArrays ++;
      //   String n = record.split("=", -2)[1];
      //   arrayName.add(n);
      // }
      if (state == PLATFORM && record.startsWith("!Platform_organism")) {
        org = record.split("=", -2)[1];
      }
      if (state == PLATFORM && record.startsWith("!Platform_sample_id")) {
        String sid = record.split("=", -2)[1];
        sid = sid.replaceAll("\\s", "");
        sampleId.add(sid);
      }
      if (state == PLATFORM && record.startsWith("!platform_table_begin")) {
        // Read platform table
        record =  br.readLine();  // ID ...
        if (record == null) {
          throw new ArrayException("No platform data - 1");
        }
        while (record != null) {
          record =  br.readLine();
          if (record == null) {
            break;
          }
          if (record.startsWith("!platform_table_end")) {
            break;
          }
        }
        // Select sample id randomly
        if (RANDOM && sampleId.size() > MAX_ARRAYS) {
          sampleIdselect = new HashSet<String>();
          tools.Permutation permutation = new tools.Permutation(100);
          int[] perm = permutation.getRandomPermutation(sampleId.size());
          for (int i =0; i < MAX_ARRAYS; i++) {
            String s = sampleId.get(perm[i]);
            System.out.println(s);
            sampleIdselect.add(s);
          }
        }
      }
      if (state == SAMPLE && record.startsWith("!Sample_data_row_count")) {
        String ns = record.split("=", -2)[1];
        ns = ns.replaceAll(" ", "");
        int n = Integer.parseInt(ns);
      }
      if (state == SAMPLE && record.startsWith("!Sample_title")) {
        sampleTitle = record.split("=", -2)[1];
      }
      if (state == SAMPLE && record.startsWith("!Sample_source_name_ch1")) {
        sampleSource = record.split("=", -2)[1];
      }
      if (sampleDescription.equals("NA") && 
          state == SAMPLE && record.startsWith("!Sample_description")) {
        sampleDescription = record.split("=", -2)[1];
      }
      if (state == SAMPLE && record.startsWith("!sample_table_begin")) {
        if (sampleIdselect == null && MAX_ARRAYS > 0 && numArrays >= MAX_ARRAYS) {
            break;
        }
        if (excludeList.contains(sampleName) ||
            sampleIdselect != null && !sampleIdselect.contains(sampleName)) {
          System.out.println(numArrays + ":" + sampleName + " Excluded");
          while (record != null) {
            record =  br.readLine();
            if (record == null) {
              break;
            }
            if (record.startsWith("!sample_table_end")) {
              break;
            }
          }
        }
        else {
          System.out.println(numArrays + ":" + sampleName + "\t" + sampleTitle);
          numArrays++;
          // Read sample table
          record =  br.readLine();  // ID ...
          if (record == null) {
            throw new ArrayException("No sample data - 1");
          }
          while (record != null) {
            record =  br.readLine();
            if (record == null) {
              throw new ArrayException("No sample data - 2");
            }
            if (record.startsWith("!sample_table_end")) {
              break;
            }
          }
          String[] values = new String[4];
          values[0] = sampleName;
          values[1] = sampleSource;
          values[2] = sampleTitle;
          values[3] = sampleDescription;
          data.add(values);
        }
      }

      record =  br.readLine();
    }
    state = DONE;

    numArrays = 4;
    numGenes = data.size();
    GeneData[] data_ = new GeneData[data.size()];
    for (int j =0; j < numGenes; j++) {
      if ( (j%100) == 0) {
        System.out.println(j);
      }
      Object[] h1 = new Object[4];
      String[] val = data.get(j);
      for (int i =0; i < val.length; i++) {
        h1[i] = val[i];
      }
      data_[j + numGeneHeader] = new GeneData(h1);
    }

    Data res = new Data(numArrays, numGenes, numGeneHeader, numArrayHeader, 
                    data_);
    res.setName(name);
    GeneNameScheme ns = new GeneNameScheme(1, org, ":", 0);
    ns.setNumMissingPoints(0);
    res.setGeneNameScheme(ns);
    System.out.println("Done");
    return res;

  }

  public static Data readSampleIds(String filename) throws Exception {
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
    int numArrays = 1;
    int numGenes = 0;
    int numGeneHeader = 0;
    int numArrayHeader = 0;

    int state = START;
    Vector<String> id = new Vector<String>();

    String record = br.readLine();
    while (record != null) {
      if (record.startsWith("^PLATFORM")) {
        state = PLATFORM;
      } 
      if (state == PLATFORM && record.startsWith("!Platform_sample_id")) {
        String sid = record.split("=", -2)[1];
        sid = sid.replaceAll("\\s", "");
        id.add(sid);
      }
      if (state == PLATFORM && record.startsWith("!platform_table_begin")) {
        break;
      }
      record =  br.readLine();
    }
    state = DONE;
    numArrayHeader = 0;
    numGeneHeader = 0;
    numGenes = id.size();
    GeneData[] data_ = new GeneData[id.size()];
    for (int j =0; j < numGenes; j++) {
      Object[] h1 = new Object[1];
      h1[0] = id.get(j);
      data_[j] = new GeneData(h1);
    }

    Data res = new Data(numArrays, numGenes, numGeneHeader, numArrayHeader, 
                    data_);
    System.out.println("Done");
    return res;

  }
};

