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

package tools.goanalysis;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.net.URL;
import java.util.Vector;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Collections;

public class AnnotationBuilder {
  String file_;
  String org_;
  Vector<Annotation> annotations_;
  LHashMap<String,String> idHash_; // DBObjectID to GOID
  LHashMap<String,String> symHash_; // DBObjectSymbol to GOID
  LHashMap<String,String> goHash_; // GOID to DBObjectID
  LHashMapUnique<String,String> symIdHash_; // DBObjectSymbol to DBObjectID
  LHashMapUnique<String,String> idSymHash_; // DBObjectID to DBObjectSymbol

  public AnnotationBuilder(String file, String org) throws Exception {
    file_ = file;
    org_ = org;
    annotations_ = new Vector<Annotation>();
    idHash_ = new LHashMap<String,String>();
    symHash_ = new LHashMap<String,String>();
    goHash_ = new LHashMap<String,String>();
    symIdHash_ = new LHashMapUnique<String,String>();
    idSymHash_ = new LHashMapUnique<String,String>();
    populateAnnotation();
  }

  public void populateAnnotation() throws Exception {
    System.out.println("Reading file " + file_);
    FileReader     fr;
    BufferedReader br;

    if (file_.startsWith("http:")) {
      URL url = new URL(file_);
      if (file_.endsWith(".gz")) {
        br = new BufferedReader(new InputStreamReader(new GZIPInputStream(url.openStream())));
      }
      else {
        br = new BufferedReader(new InputStreamReader(url.openStream()));
      }
    }
    else {
      fr = new FileReader(file_);
      br = new BufferedReader(fr);
    }

    String record = null;
    int lineno = 0;
    while ((record = br.readLine()) != null) {
      lineno++;
      if (record.startsWith("!")) {
        continue;
      }
      String[] result = record.split("\\t", -2);
      if (result.length < 15) {
        System.err.println(" *Warning: lineno " + lineno + " has undefined comulns");
      }
      else {
        result[2] = result[2].split("_", -2)[0];
        Annotation ann = new Annotation(result);
        //annotations_.add(ann);
        idHash_.put(ann.getDBObjectID(), ann.getGOID());
        symHash_.put(ann.getDBObjectSymbol(), ann.getGOID());
        goHash_.put(ann.getGOID(), ann.getDBObjectID());
        symIdHash_.put(ann.getDBObjectSymbol(), ann.getDBObjectID());
        idSymHash_.put(ann.getDBObjectID(), ann.getDBObjectSymbol());
        Enumeration<String> e = ann.getDBObjectSynonyms();
        while (e.hasMoreElements()) {
            String sym = (String) e.nextElement();
            symIdHash_.put(sym, ann.getDBObjectID());
        }
      }
    }
    System.out.println("Done");
  }

  public int getNumGenes() {
    int count = idHash_.keySet().size();
    return count;
  }
  public int getNumAnnotations(String goid) {
    int count = 0;
    if (goHash_.containsKey(goid)) {
        Vector<String> genes = goHash_.get(goid);
        count = genes.size();
    }
    return count;
  }
  public HashSet<String> getGenesFromGoid(String goid) {
    HashSet<String> genes = new HashSet<String>();
    if (goHash_.containsKey(goid)) {
        Vector<String> g = goHash_.get(goid);
        genes.addAll(g);
    }
    return genes;
  }

  /* Given a set of ObjectIDs returns the corresponding ObjectSymbols */
  public HashSet<String> findDBObjectSymbols(HashSet<String> genes) {
    HashSet<String> res = new HashSet<String>();
    Iterator<String> e = genes.iterator();
    while (e.hasNext()) {
      String gene = (String) e.next();
      if (idSymHash_.containsKey(gene)) {
        Vector<String> set = idSymHash_.get(gene);
        res.addAll(set);
      }
    }
    return res;
  }

  /* Given a set of genes returns the corresponding ObjectIDs */
  public HashSet<String> findDBObjectIDs(Vector<String> genes) {
    Vector<String> res = new  Vector<String>();
    Enumeration<String> e = genes.elements();
    while (e.hasMoreElements()) {
      String gene = (String) e.nextElement();
      if (symIdHash_.containsKey(gene)) {
        Vector<String> set = symIdHash_.get(gene);
        res.addAll(set);
      }
      if (idHash_.containsKey(gene)) {
        res.add(gene);
      }
    }
    HashSet<String> set = new HashSet<String>(res);
    return set;
  }

  public Enumeration<String> getGOIDs(HashSet<String> genes) {
    HashSet<String> goids = new HashSet<String>();
    Iterator<String> e = genes.iterator();
    while (e.hasNext()) {
        String gene = (String) e.next();
        System.out.println("-["+gene+"]");
        if (idHash_.containsKey(gene)) {
            Vector<String> ids = idHash_.get(gene);
            goids.addAll(ids);
        }
        if (symHash_.containsKey(gene)) {
            Vector<String> ids = symHash_.get(gene);
            goids.addAll(ids);
        }
    }
    Vector<String> res = new Vector<String>();
    res.addAll(goids);
    return res.elements();
  }

  public static void main(String args[]) throws Exception {
    System.out.println("File: " + args[0]);
    AnnotationBuilder ann = new AnnotationBuilder(args[0], "Mm");
  }

}
