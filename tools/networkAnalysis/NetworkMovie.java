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

package tools.networkAnalysis;

import java.awt.*;

import java.util.*;
import java.text.*;
import java.io.*;

import tools.graphs.*;
import tools.graphs.parser.GML;
import tools.microarray.*;
import tools.microarray.FileReader.PCLFileReader;

public class NetworkMovie {

    Graph graph_;
    Data  data_;
    HashSet<String> genes_;
    HashMap<String, Node> nodes_;
    HashMap<String, Integer> map_;
    Color high_, low_, zero_, missing_;
    double contrast_;

    public NetworkMovie(String pclFile, String networkFile) throws Exception {
      Vector<Graph> graphs = GML.parseGMLFile(networkFile);
      graph_ = graphs.get(0);
      data_ = PCLFileReader.readFile(pclFile); 
      data_.convertDoubles();
      GeneNameScheme ns;
      //ns = new GeneNameScheme(0, "Mm", ":", 0);
      ns = new GeneNameScheme(1, "Mm", "||", 3);
      data_.setGeneNameScheme(ns);
      genes_ = new HashSet<String>();
      nodes_ = new HashMap<String, Node>();
      map_ = new HashMap<String, Integer>();
      high_ = Color.red;
      low_ = Color.green;
      zero_ = Color.black;
      missing_ = Color.lightGray;
      contrast_ = 0.5;
    }

    public void buildGenes() {
      class GraphDFSVisitor implements DFSVisitor {
        public boolean visitBefore(Graph g, Node v, Node parent) {
          return true;
        }
        public boolean visit(Graph g, Node v, Node parent) {
          String label = (String) v.getAttribute("label");
          label = NetworkImage.trimDoubleQuotes(label);
          label = label.replaceAll("\\s", "");
          label = label.toUpperCase();
          genes_.add(label);
          nodes_.put(label, v);
          return true;
        }
        public boolean visitAfter(Graph g, Node v, Node parent) {
          return true;
        }
      };
      graph_.traverseDFS(new GraphDFSVisitor());
      int count = 0;
      for (int i = data_.getNumGeneHeader(); i < data_.getNumRows(); i++) {
        String name = data_.getGenesAt(i);
        if (name != null) {
          name = name.replaceAll("\\s", "");
          name = name.toUpperCase();
          if (genes_.contains(name)) {
            map_.put(name, new Integer(i));
            count++;
          }
        }
      }
      System.out.println("Number of genes in the Network = " + genes_.size());
      System.out.println("Number of genes found in the data = " + map_.size());
    }

    public static Color getColor(Color z, Color h, double val) {
      int red = (int) (h.getRed() * val + z.getRed() * ( 1 - val));
      int blue = (int) (h.getBlue() * val + z.getBlue() * ( 1 - val));
      int green = (int) (h.getGreen() * val + z.getGreen() * ( 1 - val));
      return new Color(red, green, blue);
    }

    public Color getColor(Double entry) {
      if (entry == null) {
        return missing_;
      }
      double val = entry.doubleValue();
      val = val / (double) contrast_;
      if (val > 1) {
        val = 1;
      }
      if (val < -1) {
        val = -1;
      }
      if (val > 0) {
        return getColor(zero_, high_, val);
      }
      else {
        return getColor(zero_, low_, -val);
      }
    }

    public static String getHexString(int c) {
      String res = Integer.toHexString(c);
      if (res.length() == 1) {
        res = "0" + res; 
      }
      return res;
    }
    public static String getColorString(Color c) {
        return "#"+getHexString(c.getRed())+getHexString(c.getGreen())+getHexString(c.getBlue());
    }

    @SuppressWarnings(value={"unchecked"}) 
    public void writeMovies(String dir, String name) throws Exception {
      buildGenes();
      int start = data_.getNumArrayHeader();
      int end = data_.getNumColumns() - 1;
      NetworkImage image = new NetworkImage(graph_);
      int index = 0;
      int lastI = start;
      for (int i = start; i <= end; i++) {
        for (int j = 0; j < 30; j++) {
          Iterator<String> itr = nodes_.keySet().iterator();
          while (itr.hasNext()) {
            String label = (String) itr.next();
            Node node = nodes_.get(label);
            HashMap<Object,Object> gr = (HashMap<Object, Object>) node.getAttribute("graphics");
            Integer loc = map_.get(label);
            if (loc != null ) {
              GeneData gene = data_.getGeneData(loc.intValue());
              Double val = (Double) data_.getGeneData(loc.intValue()).getDataAt(i);
              Double lastval = (Double) data_.getGeneData(loc.intValue()).getDataAt(lastI);
              if (val != null) {
                String cstr = getColorString(getColor(getColor(lastval), getColor(val), j/30.0));
                gr.put("fill", cstr);
              }
            }
            else {
              String cstr = getColorString(missing_);
              gr.put("fill", cstr);
            }
          }
          int copy = 2;
          if (i == end && j == 29) {
            copy = 70;
          }
          for (int k = 0; k < copy; k++) {
            DecimalFormat myFormatter = new DecimalFormat("0000");
            String output = myFormatter.format(index);
            System.out.println(dir + File.separator + name+"_"+output+".jpg");
            image.prepareImage();
            DecimalFormat dox = new DecimalFormat("Dox : 0.00");
            String doxstr = dox.format((i-3)/100.0);
            image.drawString(doxstr, 50, 50);
            image.writeImage(dir + File.separator + name+"_"+output+".jpg");
            index++;
            //return;
          }
        }
        lastI = i;
      }
    }

    public static void main(String[] arg) throws Exception {
        NetworkMovie m = new NetworkMovie(arg[0], arg[1]);
        m.writeMovies(arg[2], arg[3]);
    }
}
