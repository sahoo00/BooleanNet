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

import tools.microarray.Data;

import java.io.*;
import java.util.*;

import tools.goanalysis.GOAnalysis;
import tools.goanalysis.GOTerm;

public class SMGOAnalysis {

  String annFile_;
  String onnFile_;
  String org_;
  Data data_;
  GOAnalysis goa_;
  double pvalueThr_;

  public SMGOAnalysis(Data data, String onnFile, String annFile, String org,
    Double pval) throws StepException {
      try {
        data_ = data;
        annFile_ = annFile;
        onnFile_ = onnFile;
        org_ = org;
        pvalueThr_ = pval.doubleValue();
        goa_ = new GOAnalysis(onnFile, annFile, org, pvalueThr_);
      }
      catch(Exception e) {
        e.printStackTrace();
        throw new StepException("Problems in GO Analysis");
      }
    }

  public void writeHtml(String file, Vector<String> allNodes,
      HashMap<Integer,SMHashMapUnique<Integer,String> > map ) throws Exception {
    System.out.println(" Total number of genes : " + allNodes.size() + 
        " (" + org_ + ")" );
    for (int i = 1; i < 5; i++) {
      SMHashMapUnique<Integer,String> set = (SMHashMapUnique<Integer,String>)
        map.get(new Integer(i));
      List<Integer> key = new ArrayList<Integer>(set.keySet());
      Collections.sort(key);
      Iterator<Integer> itr = key.iterator();
      while (itr.hasNext()) {
        Integer j = (Integer) itr.next();
        HashSet<String>  hashSet= set.get(j);
        int num = 0;
        if (hashSet != null) {
          num = hashSet.size();
        }
        System.out.println("( "+ i + ")Step : " + j + " (" + num + ")");
      }
    }

    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    out.write(GOAnalysis.getHeader());

    for (int i = 1; i < 5; i++) {
      if (i == 1) {
        out.write("<li> <font size=+2> Up regulated genes (GR) </font> <ul>\n");
      }
      if (i == 2) {
        out.write("<li> <font size=+2> Down regulated genes (RG) </font> <ul>\n");
      }
      if (i == 3) {
        out.write("<li> <font size=+2> Up-Down regulated genes (GRG) </font> <ul>\n");
      }
      if (i == 4) {
        out.write("<li> <font size=+2> Down-Up regulated genes (RGR) </font> <ul>\n");
      }
      SMHashMapUnique<Integer,String> set = (SMHashMapUnique<Integer,String>)
        map.get(new Integer(i));
      List<Integer> key = new ArrayList<Integer>(set.keySet());
      Collections.sort(key);
      Iterator<Integer> itr = key.iterator();
      while (itr.hasNext()) {
        Integer j = (Integer) itr.next();
        HashSet<String>  hashSet= set.get(j);
        int num = 0;
        if (hashSet != null) {
          num = hashSet.size();
        }
        out.write("<li> <font size=+2> Step = "+ j +" </font> <ul>\n");
        System.out.println("Step : " + j + " (" + num + ")");
        System.out.println("-------------");
        if (hashSet == null) {
          continue;
        }
        Vector<String> strSet = new Vector<String>(hashSet);
        for (int k = 0; k < strSet.size(); k++) {
          System.out.println("["+strSet.get(k)+"]");
        }
        goa_.printGOTerms(out, org_, strSet, allNodes);
        out.write("</ul> </li>\n");
      }
      out.write("</ul> </li>\n");
    }
    out.write(GOAnalysis.getEnd());

    out.close();
  }

  public void writeHtml(String file, SMHashMapUnique<String,String> map) 
    throws Exception {
    Vector<String> allNodes = new Vector<String>(map.get("All"));
    System.out.println(" Total number of genes : " + allNodes.size() + 
        " (" + org_ + ")" );
    HashMap<String, GOTerm> bestGoTerms = new HashMap<String, GOTerm>();
    HashMap<String, String> bestGoTermKey = new HashMap<String, String>();
    Iterator<String> itr = map.keySet().iterator();
    while (itr.hasNext()) {
      String key = (String) itr.next();
      if (key.equals("All") || key.equals("Rest")) {
        continue;
      }
      HashSet<String>  hashSet= map.get(key);
      Vector<String> strSet = new Vector<String>(hashSet);
      for (int k = 0; k < strSet.size(); k++) {
        System.out.println("["+strSet.get(k)+"]");
      }
      Iterator<GOTerm> itr1 = goa_.getSortedGOTerms(strSet, allNodes);
      while (itr1.hasNext()) {
        GOTerm term = (GOTerm) itr1.next();
        String id = term.getID();
        if (!bestGoTerms.containsKey(id)) {
          bestGoTerms.put(id, term);
          bestGoTermKey.put(id, key);
        }
        GOTerm old = bestGoTerms.get(id);
        if (term.getPvalueFwer() < old.getPvalueFwer()) {
          bestGoTerms.put(id, term);
          bestGoTermKey.put(id, key);
        }
      }
    }
    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    out.write(GOAnalysis.getHeader());
    
    List<GOTerm> list = new ArrayList<GOTerm>(bestGoTerms.values());
    Iterator<GOTerm> itr2 = goa_.getSortedGOTerms(list);
    out.write("<li> <font size=+2> GOTerms </font> <ul>\n");
    out.write("<li> <font size=+2> Best GOTerms </font> <ul>\n");
    goa_.printGOTerms(out, org_, itr2, bestGoTermKey);
    out.write("</ul> </li>\n");
    out.write("</ul> </li>\n");
    out.write(GOAnalysis.getEnd());

    out.close();
  }
}
