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

import tools.microarray.ArrayException;
import tools.microarray.GeneSet;
import java.io.*;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.util.LinkedList;
import java.util.ListIterator;
import java.lang.Math;

public class GeneSetFileReader {

  public static GeneSet readFile(String filename) throws Exception {
    if (filename.endsWith(".tab") || filename.endsWith(".tab.gz")) {
      return readTABFile(filename);
    }
    if (filename.endsWith(".gmt") || filename.endsWith(".gmt.gz")) {
      return readGMTFile(filename);
    }
    throw new ArrayException("Unsupported gene set file type in " + filename);
  }

  public static GeneSet readTABFile(String filename) throws Exception {
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
    GeneSet res = new GeneSet();

    String record = br.readLine();
    if (record == null) {
      throw new ArrayException("No header - 1");
    }
    String[] set = record.split("\\t", -2); // -2 : Don't discard trailing nulls
    for (int i = 1; i < set.length; i++) {
      res.addName(set[i]);
    }

    while ((record = br.readLine()) != null) {
      lineno++;
      String[] result = record.split("\\t", -2);
      for (int i = 1; i < result.length; i++) {
        int num = 0;
        try {
            num = Integer.parseInt(result[i]);
        }
        catch(Exception e) {
        }
        if (num > 0) {
          res.add(set[i], result[0]);
        }
      }
    }
    res.print();
    return res;
  }

  public static GeneSet readGMTFile(String filename) throws Exception {
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
    GeneSet res = new GeneSet();
    String record = "";

    while ((record = br.readLine()) != null) {
      lineno++;
      String[] result = record.split("\\t", -2);
      for (int i = 2; i < result.length; i++) {
        res.add(result[0], result[1], result[i]);
      }
    }
    //res.print();
    return res;
  }

};

