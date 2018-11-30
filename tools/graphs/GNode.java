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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;
import java.util.Iterator;

public class GNode implements Node {

  String id_;
  HashMap<Object,Object> attributes_;

  HashSet<Edge> edges_;

  public GNode(String id) {
    id_ = id;
    attributes_ = new HashMap<Object,Object>();
    edges_ = new HashSet<Edge>();
  }

  public void setAttribute(Object key, Object val) {
    attributes_.put(key, val);
  }

  public void addChild(Node v) { 
    GEdge o = new GEdge(this, v, false);
    edges_.add(o);
    v.addIncomingEdges(o);
  }

  public void addIncomingEdges(Edge e) { edges_.add(e); }

  public Object getAttribute(Object key) { return attributes_.get(key); }
  public String getID() { return id_;}
  public int getNumChildren() { 
    int count = 0;
    Iterator<Edge> itr = edges_.iterator();
    while (itr.hasNext()) {
      Edge e = (Edge) itr.next();
      if (e.isDirected()) {
        if (e.isSource(this)) {
          count++;
        }
      }
      else {
        count++;
      }
    }
    return count;
  }
  public int getNumParents() {
    int count = 0;
    Iterator<Edge> itr = edges_.iterator();
    while (itr.hasNext()) {
      Edge e = (Edge) itr.next();
      if (e.isDirected()) {
        if (e.isDestination(this)) {
          count++;
        }
      }
      else {
        count++;
      }
    }
    return count;
  }

  public boolean isRoot() {
    int countd = 0;
    int countu = 0;
    Iterator<Edge> itr = edges_.iterator();
    while (itr.hasNext()) {
      Edge e = (Edge) itr.next();
      if (e.isDirected()) {
        if (e.isDestination(this)) {
          countd++;
        }
      }
      else {
        countu++;
      }
    }
    return (countd <= 0 && countu <= 1);
  }


  public Enumeration<Node> getParents() { 
    Vector<Node> res = new Vector<Node>();
    Iterator<Edge> itr = edges_.iterator();
    while (itr.hasNext()) {
      Edge e = (Edge) itr.next();
      if (e.isDirected()) {
        if (e.isDestination(this)) {
          Node s = e.getSource();
          res.add(s);
        }
      }
      else {
        Node s = e.getSource(this);
        res.add(s);
      }
    }
    return res.elements();
  }
  public Enumeration<Node> getChildren() {
    Vector<Node> res = new Vector<Node>();
    Iterator<Edge> itr = edges_.iterator();
    while (itr.hasNext()) {
      Edge e = (Edge) itr.next();
      if (e.isDirected()) {
        if (e.isSource(this)) {
          Node s = e.getDestination();
          res.add(s);
        }
      }
      else {
        Node s = e.getSource(this);
        res.add(s);
      }
    }
    return res.elements();
  }
  public Enumeration<Edge> getOutGoingEdges() {
    Vector<Edge> res = new Vector<Edge>();
    Iterator<Edge> itr = edges_.iterator();
    while (itr.hasNext()) {
      Edge e = (Edge) itr.next();
      if (e.isDirected()) {
        if (e.isSource(this)) {
          res.add(e);
        }
      }
      else {
        res.add(e);
      }
    }
    return res.elements();
  }
  public String toString() { 
    StringBuffer res = new StringBuffer();
    res.append("[\n\tid " + id_ + "\n");
    Iterator<Object> itr = attributes_.keySet().iterator();
    while (itr.hasNext()) {
      Object key = itr.next();
      Object val = attributes_.get(key);
      res.append("\t" + key + " " + val.toString() + "\n");
    }
    res.append("]\n");
    return res.toString();
  }

  public void print() { 
    System.out.print(toString());
  }

  void traverseDFS_(DFSVisitor v, HashSet<Node> visited, Node n, Node p) {
    if (visited.contains(n)) {
      return;
    }
    visited.add(n);
    v.visitBefore(null, n, p);
    v.visit(null, n, p);
    Enumeration<Node> children = n.getChildren();
    if (children != null) {
      while (children.hasMoreElements()) {
        Node c = (Node) children.nextElement();
        traverseDFS_(v, visited, c, n);
      }
    }
    v.visitAfter(null, n, p);
  }
  public void traverseDFS(DFSVisitor v) {
    HashSet<Node> visited = new HashSet<Node>();
    traverseDFS_(v, visited, this, null);
  }
}
