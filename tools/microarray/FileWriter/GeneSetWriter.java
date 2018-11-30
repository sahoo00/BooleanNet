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

package tools.microarray.FileWriter;

import tools.microarray.GeneSet;
import tools.goanalysis.GOAnalysis;
import tools.microarray.StepMiner.SMHashMapUnique;
import tools.microarray.ArrayException;
import java.io.*;
import java.util.*;

public class GeneSetWriter {

  public static void writeFile(GeneSet set, String filename, String org) throws Exception {
      writeHTMLFile(set, filename, org);
  }

  public static void writeFile(GeneSet set, String filename) throws Exception {
    int done = 0;
    if (filename.endsWith(".tab")) {
      writeTABFile(set, filename);
      done = 1;
    }
    if (filename.endsWith(".gmt")) {
      writeGMTFile(set, filename);
      done = 1;
    }
    if (filename.endsWith(".gxa")) {
      writeGXAFile(set, filename);
      done = 1;
    }
    if (done == 0) {
      throw new ArrayException("Unsupported gene set file type in " + filename);
    }
  }

  public static void writeHTMLFile(GeneSet set, String file, String org) throws IOException {
    System.out.println("Writing file " + file);
    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    out.write(GOAnalysis.getHeader());
    out.write("<li> <font size=+1> Gene Sets </font> ("+set.getTotalNum()+")<ul>\n");
    out.write("<li> <font size=+1> Set </font> <ul>\n");
    Iterator<String> itr = set.iterator();
    while (itr.hasNext()) {
      String tag = (String) itr.next();
      String desc = set.getDescription(tag);
      HashSet<String>  hashSet= set.getSet(tag);
      int k = hashSet.size();
      out.write("<li> <font size=+1> "+ tag +"</font> " +
          desc +
          "(" + k + ")" +
          "<ul> <li> <table border=0><tr>\n");
      Iterator<String> gitr = hashSet.iterator();
      int count = -1;
      while (gitr.hasNext()) {
        String gene = (String) gitr.next();
        count ++;
        if ( (count % 10) == 0) {
          out.write("</tr><tr>\n");
        }
        String link = GOAnalysis.getLink(org, gene);
        out.write("<td> <a target=\"_blank\" href=\""+ link + "\"> "+ gene + "</a></td>\n");
      }
      out.write("</tr></table> </li></ul> </li>\n");
    }
    out.write("</ul> </li>\n");
    out.write("</ul> </li>\n");
    out.write(GOAnalysis.getEnd());

    out.close();
  }

  public static void writeGMTFile(GeneSet set, String file) throws IOException {
    System.out.println("Writing file " + file);
    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    Iterator<String> itr = set.iterator();
    while (itr.hasNext()) {
      String tag = (String) itr.next();
      String desc = set.getDescription(tag);
      out.write(tag+"\t"+desc+"\t");
      HashSet<String>  hashSet= set.getSet(tag);
      Iterator<String> itr1 = hashSet.iterator();
      while (itr1.hasNext()) {
        String name = (String) itr1.next();
        out.write(name+"\t");
      }
      out.write("\n");
    }
    out.close();
    System.out.println("Done");
  }

  public static void writeTABFile(GeneSet set, String file) throws IOException {
    System.out.println("Writing file " + file);
    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    // Writing Key Set
    out.write("Gene");
    int groupid = 0;
    int totalNum = set.getTotalNum();
    SMHashMapUnique<Integer,Integer> rmap = new
      SMHashMapUnique<Integer,Integer>();
    HashMap<String, Integer> hash = new HashMap<String, Integer>();
    String[] list = new String[totalNum];
    int index = 0;
    Iterator<String> itr = set.iterator();
    while (itr.hasNext()) {
      groupid++;
      String tag = (String) itr.next();
      HashSet<String>  hashSet= set.getSet(tag);
      out.write("\t" + tag);
      Iterator<String> itr1 = hashSet.iterator();
      while (itr1.hasNext()) {
        String name = (String) itr1.next();
        if (!hash.containsKey(name)) {
          Integer gid = new Integer(index);
          hash.put(name,gid);
          list[index++] = name;
        }
        Integer gid = hash.get(name);
        rmap.put(gid, new Integer(groupid));
      }
    }
    out.write("\n");

    for (int i=0; i < list.length; i++) {
      out.write(list[i]);
      HashSet<Integer> genes = rmap.get(new Integer(i));
      int[] bitmap = new int[groupid];
      for (int j=0; j < groupid; j++) {
        bitmap[j] = 0;
      }
      Iterator<Integer> itr1 = genes.iterator();
      while (itr1.hasNext()) {
        Integer name = (Integer) itr1.next();
        bitmap[name.intValue()-1] = 1;
      }
      for (int j=0; j < groupid; j++) {
        out.write("\t" + bitmap[j]);
      }
      out.write("\n");
    }

    out.close();
    System.out.println("Done");
  }

  public static void writeGXAFile(GeneSet set, String file) throws IOException {
    System.out.println("Writing file " + file);
    BufferedWriter out = new BufferedWriter(new FileWriter(file));
    // Writing Key Set
    out.write("<?xml version='1.0' encoding='iso-8859-1'?>\n\n");
    out.write("<GeneXPress>\n<GeneXPressAttributes>\n");
    out.write("<Attributes Id=\"0\" Name=\"StepMiner\">\n");
    int groupid = 0;
    int totalNum = set.getTotalNum();
    SMHashMapUnique<Integer,Integer> rmap = new
      SMHashMapUnique<Integer,Integer>();
    HashMap<String, Integer> hash = new HashMap<String, Integer>();
    String[] list = new String[totalNum];
    int index = 0;
    Iterator<String> itr = set.iterator();
    while (itr.hasNext()) {
      groupid++;
      String tag = (String) itr.next();
      HashSet<String>  hashSet= set.getSet(tag);
      out.write("<Attribute Name=\"" +tag+ "\" Id=\"" +groupid+"\"");
      out.write(" Counts=\"" + totalNum + " " + hashSet.size() + "\"");
      out.write(" Value=\"0 1\" />");
      out.write("\n");
      Iterator<String> itr1 = hashSet.iterator();
      while (itr1.hasNext()) {
        String name = (String) itr1.next();
        if (!hash.containsKey(name)) {
          Integer gid = new Integer(index);
          hash.put(name,gid);
          list[index++] = name;
        }
        Integer gid = hash.get(name);
        rmap.put(gid, new Integer(groupid));
      }
    }
    out.write("</Attributes>\n");
    out.write("</GeneXPressAttributes>\n");
    out.write("<GeneXPressObjects>\n<Objects Type=\"Genes\">\n");

    for (int i=0; i < list.length; i++) {
      out.write("<Gene Id=\"" + i + "\" ORF=\"" + list[i] + "\">\n");
      out.write("<Attributes AttributesGroupId=\"0\" ");
      out.write("Type=\"Full\" Value=\"");
      HashSet<Integer> genes = rmap.get(new Integer(i));
      int[] bitmap = new int[groupid];
      for (int j=0; j < groupid; j++) {
        bitmap[j] = 0;
      }
      Iterator<Integer> itr1 = genes.iterator();
      while (itr1.hasNext()) {
        Integer name = (Integer) itr1.next();
        bitmap[name.intValue()-1] = 1;
      }
      out.write("" + bitmap[0]);
      for (int j=1; j < groupid; j++) {
        out.write(";" + bitmap[j]);
      }
      out.write("\"/>\n</Gene>\n");
    }

    out.write("</Objects>\n</GeneXPressObjects>\n</GeneXPress>\n");
    out.close();
    System.out.println("Done");
  }

};

