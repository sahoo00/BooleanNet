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

import tools.goanalysis.parser.OBO;
import tools.goanalysis.parser.ParseException;
import tools.goanalysis.StanzaMap;
import tools.graphs.DAGGraph;
import tools.graphs.DAGNode;

import java.util.Vector;
import java.util.Enumeration;
import java.util.HashMap;

public class Ontology {
  String file_;
  Vector<StanzaMap> stanzas_;
  DAGGraph graph_;
  HashMap<String,DAGNode> hashNode_;

  public Ontology(String file) throws Exception {
    file_ = file;
    stanzas_ = OBO.parseOBOFile(file);
    buildGraph();
  }

  public DAGGraph getGraph() { return graph_; }
  public DAGNode getNode(String id) { return hashNode_.get(id); }
  public int size() { return stanzas_.size(); }

  public void buildGraph() {
    hashNode_ = new HashMap<String,DAGNode>();
    Enumeration<StanzaMap> e = stanzas_.elements();
    graph_ = new DAGGraph();
    DAGNode others = new DAGNode("Others");
    others.setAttribute("name", "Other genes: No Match");
    others.setAttribute("namespace", "biological_process");
    hashNode_.put("Others", others);
    graph_.addRoot(others);
    while (e.hasMoreElements()) {
      StanzaMap st = (StanzaMap) e.nextElement();
      if (st.isTerm()) {
        String id = st.getID();
        DAGNode n = new DAGNode(id);
        n.setAttribute("name", st.getName());
        n.setAttribute("namespace", st.getNamespace());
        hashNode_.put(id, n);
      }
    }
    e = stanzas_.elements();
    while (e.hasMoreElements()) {
      StanzaMap st = (StanzaMap) e.nextElement();
      if (st.isTerm()) {
        String id = st.getID();
        DAGNode n = getNode(id);
        int numParents = 0;
        Enumeration<String> e1 = st.getParents();
        while (e1.hasMoreElements()) {
          String pid = (String) e1.nextElement();
          if (hashNode_.containsKey(pid)) {
            DAGNode p = getNode(pid);
            p.addChild(n);
            numParents ++;
          }
        }
        if (numParents == 0) {
          graph_.addRoot(n);
        }
      }
    }
  }

  public static void main(String args[]) throws Exception {
    System.out.println("File: " + args[0]);
    Ontology onn = new Ontology(args[0]);
  }

}
