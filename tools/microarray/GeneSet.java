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
import tools.microarray.StepMiner.SMHashMapUnique;
import tools.microarray.FileReader.GeneSetFileReader;
import tools.microarray.FileWriter.GeneSetWriter;

public class GeneSet {

  Vector<String> setNames_;
  HashMap<String, String> descriptions_;
  HashSet<String> allGenes_;
  SMHashMapUnique<String, String> map_;

  public GeneSet() {
    setNames_ = new Vector<String>();
    descriptions_ = new HashMap<String, String>();
    allGenes_ = new HashSet<String>();
    map_ = new SMHashMapUnique<String, String>();
  }

  public HashSet<String> getSet(String name) {
    return map_.get(name);
  }

  public String getDescription(String name) {
    return descriptions_.get(name);
  }

  public Iterator<String> iterator() {
    return setNames_.iterator();
  }

  public int getTotalNum() {
    return allGenes_.size();
  }

  public HashSet<String> getAllGenes() {
    return allGenes_;
  }

  public Vector<String> getAllSets() {
    return setNames_;
  }

  public void addName(String name) {
    addName(name, name);
  }

  public void addName(String name, String desc) {
    if (!descriptions_.containsKey(name)) {
      setNames_.add(name);
      descriptions_.put(name, desc);
    }
  }

  public void add(String name, String gene) {
    add(name, name, gene);
  }

  public void add(String name, String desc, String gene) {
    addName(name, desc);
    map_.put(name, gene);
    allGenes_.add(gene);
  }

  public void add(String name, HashSet<String> set) {
    add(name, name, set);
  }

  public void add(String name, String desc, HashSet<String> set) {
    addName(name, desc);
    Iterator<String> itr = set.iterator();
    while(itr.hasNext()) {
      String gene = (String) itr.next();
      map_.put(name, gene);
      allGenes_.add(gene);
    }
  }

  public void reorder(Integer[] order) {
    if (order.length != setNames_.size()) {
       System.err.println("Error in ordering\n");
       return;
    }
    Vector<String> newSets = new Vector<String>();
    for (int j =0; j < order.length ; j++) {
      String name = setNames_.get(order[j].intValue());
      // System.out.println(name);
      newSets.add(name);
    }
    setNames_ = newSets;
  }

  public static GeneSet merge(GeneSet a, GeneSet b) {
    GeneSet res = new GeneSet();
    Iterator<String> itr = a.iterator();
    while(itr.hasNext()) {
      String str = (String) itr.next();
      String desc = a.getDescription(str);
      HashSet<String> set = a.getSet(str);
      res.add(str, desc, set);
    }
    itr = b.iterator();
    while(itr.hasNext()) {
      String str = (String) itr.next();
      String desc = b.getDescription(str);
      HashSet<String> set = b.getSet(str);
      res.add(str, desc, set);
    }
    return res;
  }

  public static GeneSet deleteGenes(GeneSet a, HashSet<String> genes) {
    GeneSet res = new GeneSet();
    Iterator<String> itr = a.iterator();
    while(itr.hasNext()) {
      String str = (String) itr.next();
      String desc = a.getDescription(str);
      HashSet<String> set = a.getSet(str);
      set.removeAll(genes);
      if (set.size() > 0) {
        res.add(str, desc, set);
      }
    }
    return res;
  }

  public static GeneSet toUpperCase(GeneSet a) {
    GeneSet res = new GeneSet();
    Iterator<String> itr = a.iterator();
    while(itr.hasNext()) {
      String str = (String) itr.next();
      String desc = a.getDescription(str);
      HashSet<String> set = a.getSet(str);
      Iterator<String> itr1 = set.iterator();
      while(itr1.hasNext()) {
        String gene = (String) itr1.next();
        gene = gene.toUpperCase();
        res.add(str, desc, gene);
      }
    }
    return res;
  }

  public String toString() {
    StringBuffer res = new StringBuffer();
    Iterator<String> itr = iterator();
    while(itr.hasNext()) {
      String str = (String) itr.next();
      String desc = getDescription(str);
      res.append(str + ":" + desc + "\n\t");
      HashSet<String> set = getSet(str);
      Iterator<String> itr1 = set.iterator();
      int i = 0;
      while(itr1.hasNext()) {
        String s = (String) itr1.next();
        res.append(s + "\t");
        i++;
        if ( (i % 10) == 0) {
          res.append("\n\t");
        }
      }
      res.append("\n");
    }
    return res.toString();
  }

  public void print() {
    System.out.println(toString());
  }

  public static GeneSet readFile(String file) throws Exception {
    return GeneSetFileReader.readFile(file);
  }

  public void writeSetFile(String file) throws Exception {
    writeFile(this, file);
  }

  public static void writeFile(GeneSet set, String file) throws Exception {
    GeneSetWriter.writeFile(set, file);
  }

  public static void writeFile(GeneSet set, String file, String org) throws Exception {
    GeneSetWriter.writeFile(set, file, org);
  }

  public static void main(String arg[]) throws Exception {
    readFile(arg[0]);
  }

};

