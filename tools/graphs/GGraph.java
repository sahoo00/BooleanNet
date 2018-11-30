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

package tools.graphs;

import java.util.Enumeration;
import java.util.HashSet;
import java.util.Iterator;

public class GGraph implements Graph {

  HashSet<Node> roots_;

  public GGraph() {
    roots_ = new HashSet<Node>();
  }

  public Enumeration<Node> getNodes() { return null;}
  public Enumeration<Edge> getEdges() { return null;}

  public void addRoot(Node v) { roots_.add(v); }

  void traverseDFS_(DFSVisitor v, HashSet<Node> visited, Node n, Node p) {
    if (visited.contains(n)) {
      return;
    }
    visited.add(n);
    v.visitBefore(this, n, p);
    v.visit(this, n, p);
    Enumeration<Node> children = n.getChildren();
    if (children != null) {
      while (children.hasMoreElements()) {
        Node c = (Node) children.nextElement();
        traverseDFS_(v, visited, c, n);
      }
    }
    v.visitAfter(this, n, p);
  }

  public void traverseDFS(DFSVisitor v) {
    HashSet<Node> visited = new HashSet<Node>();
    Iterator<Node> itr = roots_.iterator();
    while (itr.hasNext()) {
      Node n = (Node) itr.next();
      traverseDFS_(v, visited, n, null);
    }
  }

  public void print() {
    Iterator<Node> itr = roots_.iterator();
    while (itr.hasNext()) {
      Node n = (Node) itr.next();
      n.print();
    }
  }

}
