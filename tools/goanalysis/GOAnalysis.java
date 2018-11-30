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

import java.util.*;
import java.io.*;
import tools.graphs.*;
import tools.microarray.GeneSet;

public class GOAnalysis {
  String ontologyFile_;
  String annotationFile_;
  String org_;
  double pvalueThr_;

  AnnotationBuilder annBuilder_;
  Ontology  onnBuilder_;

  HashSet<String> foundGenes_; // tmp variable for the DFS visitor
  Vector<GOTerm> resTerm_;
  int totalGenes_;
  HashSet<String> significantGOIDs_; // tmp variable for the DFS visitor
  HashSet<String> significantGOIDgenes_; // tmp variable for the DFS visitor

  public GOAnalysis(String ontologyFile, String annotationFile, String org,
    double pval) throws Exception {
      ontologyFile_ = ontologyFile;
      annotationFile_ = annotationFile;
      org_ = org;
      pvalueThr_ = pval;
      annBuilder_  = new AnnotationBuilder(annotationFile, org);
      onnBuilder_  = new Ontology(ontologyFile);
      calculateGeneSets();
    }

  public void calculateGeneSets() {
    DAGGraph graph = onnBuilder_.getGraph();
    // Traverse the graph
    class GraphDFSVisitor implements DFSVisitor {
        public boolean visitBefore(Graph g, Node v, Node parent) {
            return true;
        }
        public boolean visit(Graph g, Node v, Node parent) {
            String goid = v.getID();
            HashSet set = (HashSet) v.getAttribute("geneSets");
            if (set == null) {
              set = annBuilder_.getGenesFromGoid(goid);
            }
            v.setAttribute("geneSets", set);
            return true;
        }
        public boolean visitAfter(Graph g, Node v, Node parent) {
            Enumeration<Node> e = v.getChildren();
            HashSet<String> res = new HashSet<String>();
            while (e.hasMoreElements()) {
                Node c = (Node) e.nextElement();
                HashSet childSet = (HashSet) c.getAttribute("geneSets");
                Iterator itr = childSet.iterator();
                while (itr.hasNext()) {
                    String gene = (String) itr.next();
                    res.add(gene);
                }
            }
            HashSet set = (HashSet) v.getAttribute("geneSets");
            Iterator itr = set.iterator();
            while (itr.hasNext()) {
                String gene = (String) itr.next();
                res.add(gene);
            }
            v.setAttribute("geneSets", res);
            return true;
        }
    };
    graph.traverseDFS(new GraphDFSVisitor());
  }

  public Vector<GOTerm> getGOTerms(
      Vector<String> geneNames, int totalGenes) {
    HashSet<String> genes = annBuilder_.findDBObjectIDs(geneNames);
    return getGOTermsSignificantMultipleCorrection(genes, totalGenes);
  }

  public static Iterator<GOTerm> getSortedGOTerms(List<GOTerm> list) {
    class GOTermComparator implements Comparator<GOTerm> {
      public GOTermComparator() {
        super();
      }
      public int compare(GOTerm s1, GOTerm s2) {
        if (s1.getPvalue() > s2.getPvalue()) {
          return 1;
        }
        return -1;
      }
    };
    Collections.sort(list, new GOTermComparator());
    Iterator<GOTerm> itr1 = list.iterator();
    return itr1;
  }

  public Iterator<GOTerm> getSortedGOTerms(
      Vector<String> geneNames, Vector<String> allNodes) {
    Vector<GOTerm> terms = getGOTerms(geneNames, allNodes.size());
    List<GOTerm> list = new ArrayList<GOTerm>(terms);
    return getSortedGOTerms(list);
  }

  public Vector<GOTerm> getAllGOTerms() {
    resTerm_ = new Vector<GOTerm>();
    totalGenes_ = annBuilder_.getNumGenes();
    DAGGraph graph = onnBuilder_.getGraph();
    // Traverse the graph
    class GODFSVisitor implements DFSVisitor {
      public boolean visitBefore(Graph g, Node v, Node parent) {
        return true;
      }
      public boolean visit(Graph g, Node v, Node parent) {
        return true;
      }
      public boolean visitAfter(Graph g, Node v, Node parent) {
        String goid = v.getID();
        HashSet set = (HashSet) v.getAttribute("geneSets");
        int totalGenesGOID = set.size();
        int numGenes = totalGenesGOID;
        HashSet<String> intersection = new HashSet<String>();
        Iterator itr = set.iterator();
        while (itr.hasNext()) {
            String name = (String) itr.next();
            intersection.add(name);
        }
        HashSet<String> syms = annBuilder_.findDBObjectSymbols(intersection);
        int numGenesGOID = intersection.size();
        GOTerm term = new GOTerm((DAGNode)v, syms, totalGenes_, 
            totalGenesGOID, numGenes,  numGenesGOID);
        resTerm_.add(term);
        return true;
      }
    };
    graph.traverseDFS(new GODFSVisitor());
    return resTerm_;
  }

  public Vector<GOTerm> getGOTermsSignificant(
      HashSet<String> geneNames, int totalGenes) {
    resTerm_ = new Vector<GOTerm>();
    int numGenes = geneNames.size();
    totalGenes_ = totalGenes;
    if (totalGenes_ < annBuilder_.getNumGenes()) {
      System.out.print(" *Warning: Given Total genes ["+totalGenes_+"] < ");
      totalGenes_ = annBuilder_.getNumGenes();
      System.out.println(totalGenes_);
    }
    foundGenes_ = geneNames;
    DAGGraph graph = onnBuilder_.getGraph();
    // Traverse the graph
    class GODFSVisitor implements DFSVisitor {
      public boolean visitBefore(Graph g, Node v, Node parent) {
        return true;
      }
      public boolean visit(Graph g, Node v, Node parent) {
        return true;
      }
      public boolean visitAfter(Graph g, Node v, Node parent) {
        String goid = v.getID();
        HashSet set = (HashSet) v.getAttribute("geneSets");
        int totalGenesGOID = set.size();
        int numGenes = foundGenes_.size();
        HashSet<String> intersection = new HashSet<String>(foundGenes_);
        intersection.retainAll(set);
        HashSet<String> syms = annBuilder_.findDBObjectSymbols(intersection);
        int numGenesGOID = intersection.size();
        if (numGenesGOID > 1) {
          GOTerm term = new GOTerm((DAGNode)v, syms, totalGenes_, 
              totalGenesGOID, numGenes,  numGenesGOID);
          resTerm_.add(term);
        }
        return true;
      }
    };
    graph.traverseDFS(new GODFSVisitor());
    return resTerm_;
  }

  // Remove irrelavent GO terms
  public List<GOTerm> reduceList(HashSet<String> geneNames, List<GOTerm> list) {
    List<GOTerm> list1 = new ArrayList<GOTerm>();
    Iterator<GOTerm> itr = list.iterator();
    while (itr.hasNext()) {
      GOTerm term = (GOTerm) itr.next();
      DAGNode n = term.getDAGNode();
      HashSet set = (HashSet) n.getAttribute("geneSets");
      HashSet<String> intersection = new HashSet<String>(geneNames);
      intersection.retainAll(set);
      int numGenes = intersection.size();
      int maxCgenes = 0;
      // for comparison with intersection
      Enumeration<Node> e = n.getChildren();
      while (e.hasMoreElements()) {
        Node c = (Node) e.nextElement();
        String c_goid = c.getID();
        HashSet childSet = (HashSet) c.getAttribute("geneSets");
        HashSet<String> cint= new HashSet<String>(geneNames);
        cint.retainAll(childSet);
        if (maxCgenes < cint.size()) {
          maxCgenes = cint.size();
        }
        intersection.removeAll(childSet);
      }
      if (intersection.size() > 0 || numGenes > 2 * maxCgenes) {
        list1.add(term);
      }
    }
    return list1;
  }

  public Vector<GOTerm> getGOTermsSignificantMultipleCorrection (
      HashSet<String> geneNames, int totalGenes) {
    Vector<GOTerm> terms = getGOTermsSignificant(geneNames, totalGenes);
    List<GOTerm> list = new ArrayList<GOTerm>(terms);

    // Descending order
    class GOTermComparatorDes implements Comparator<GOTerm> {
      public int compare(GOTerm s1, GOTerm s2) {
        if (s1.getPvalue() < s2.getPvalue()) {
          return 1;
        }
        return -1;
      }
    };
    // Ascending order
    class GOTermComparatorAsc implements Comparator<GOTerm> {
      public int compare(GOTerm s1, GOTerm s2) {
        if (s1.getPvalue() > s2.getPvalue()) {
          return 1;
        }
        return -1;
      }
    };

    // Calculation of Bonferoni adjusted pvalue
    // p < alpha/n   (First i)
    // Adjusted p-val = p_i * n
    Collections.sort(list, new GOTermComparatorDes());
    List<GOTerm> list1 = new ArrayList<GOTerm>();
    // Calculation of FDR adjusted pvalue
    // p < (i/n) * alpha   (First i)
    // Adjusted p-val = min{k=i..n}{1, n/k * p_k}
    int n = terms.size();
    int i = n;
    double minpval = 1;
    Iterator<GOTerm> itr = list.iterator();
    while (itr.hasNext()) {
      GOTerm term = (GOTerm) itr.next();
      double pval = term.getPvalue();
      double adp = pval * n/i;
      if (minpval > adp) {
        minpval = adp;
      }
      double bonferroniPval = pval * n;
      if (bonferroniPval > 1) {
        bonferroniPval = 1;
      }
      term.setPvalueFdr(minpval);
      term.setPvalueFwer(bonferroniPval);
      if (bonferroniPval < pvalueThr_) {
        list1.add(term);
      }
      i--;
    }

    list1 = reduceList(geneNames, list1);

    // Constructing Other genes
    significantGOIDs_ = new HashSet<String>();
    significantGOIDgenes_ = new HashSet<String>();
    itr = list1.iterator();
    while (itr.hasNext()) {
      GOTerm term = (GOTerm) itr.next();
      significantGOIDs_.add(term.getID());
      significantGOIDgenes_.addAll(term.getGeneHash());
    }
    HashSet<String> syms = annBuilder_.findDBObjectSymbols(geneNames);
    syms.removeAll(significantGOIDgenes_);
    int numGenes = geneNames.size();
    int numGenesGOID = 0;
    int totalGenesGOID = syms.size();
    DAGNode others = onnBuilder_.getNode("Others");
    GOTerm term = new GOTerm(others, syms, totalGenes_, 
        totalGenesGOID, numGenes,  numGenesGOID);
    list1.add(term);

    Collections.sort(list1, new GOTermComparatorAsc());
    Vector<GOTerm> res = new Vector<GOTerm>(list1);

    return res;
  }

  public Vector<GOTerm> getGOTermsOld(
        Vector<String> geneNames, int totalGenes) {
    resTerm_ = new Vector<GOTerm>();
    HashSet<String> genes = annBuilder_.findDBObjectIDs(geneNames);
    int numGenes = genes.size();
    totalGenes_ = totalGenes;
    if (totalGenes_ < annBuilder_.getNumGenes()) {
        System.out.print(" *Warning: Given Total genes ["+totalGenes_+"] < ");
        totalGenes_ = annBuilder_.getNumGenes();
        System.out.println(totalGenes_);
    }
    foundGenes_ = genes;
    significantGOIDs_ = new HashSet<String>();
    significantGOIDgenes_ = new HashSet<String>();
    DAGGraph graph = onnBuilder_.getGraph();
    // Traverse the graph
    class GODFSVisitor implements DFSVisitor {
        public boolean visitBefore(Graph g, Node v, Node parent) {
            return true;
        }
        public boolean visit(Graph g, Node v, Node parent) {
            return true;
        }
        public boolean visitAfter(Graph g, Node v, Node parent) {
            String goid = v.getID();
            HashSet set = (HashSet) v.getAttribute("geneSets");
            int totalGenesGOID = set.size();
            int numGenes = foundGenes_.size();
            HashSet<String> intersection = new HashSet<String>(foundGenes_);
            intersection.retainAll(set);
            HashSet<String> syms = annBuilder_.findDBObjectSymbols(intersection);
            int numGenesGOID = syms.size();
            int numIds = intersection.size();
            if (numGenesGOID > 1) {
              boolean include = true;
              Vector<Integer> score = new Vector<Integer>();
              int sum = numGenesGOID;
              // for comparison with intersection
              HashSet<String> down = new HashSet<String>();
              boolean significantChild = false;
              score.add(new Integer(numGenesGOID));
              Enumeration<Node> e1 = v.getChildren();
              System.out.println("Children contribution:");
              while (e1.hasMoreElements()) {
                  Node c = (Node) e1.nextElement();
                  String c_goid = c.getID();
                  HashSet childSet = (HashSet) c.getAttribute("geneSets");
                  HashSet<String> intersect= new HashSet<String>(foundGenes_);
                  intersect.retainAll(childSet);
                  down.addAll(intersect);
                  int n = intersect.size();
                  sum +=n;
                  score.add(new Integer(n));
                  System.out.print(" " + n);
                  if (significantGOIDs_.contains(c_goid)) {
                    System.out.print("S");
                    significantChild = true;
                  }
                  else if (n > 0) {
                    System.out.print("N");
                  }
              }
              System.out.println();
              double infoscore = 0;
              Enumeration<Integer> e2 = score.elements();
              while (e2.hasMoreElements()) {
                  Integer cVal = (Integer) e2.nextElement();
                  int c = cVal.intValue();
                  if ( c > 0) {
                    infoscore += -(c*1.0/sum) * Math.log(c*1.0/sum)/Math.log(2);
                  }
              }
              if (significantChild && down.size() == numIds) {
                include = false;
              }
              System.out.println("#genes : " + down.size() + " " + numIds +
                  " " + infoscore + " " +
                  Math.log(score.size())/Math.log(2));
              if (include) {
                significantGOIDs_.add(goid);
                significantGOIDgenes_.addAll(intersection);
                GOTerm term = new GOTerm((DAGNode)v, syms, totalGenes_, 
                    totalGenesGOID, numGenes,  numGenesGOID);
                if (term.getPvalue() < pvalueThr_) {
                  resTerm_.add(term);
                  term.print();
                }
              }
            }
            return true;
        }
    };
    graph.traverseDFS(new GODFSVisitor());
    // Construct others
    foundGenes_.removeAll(significantGOIDgenes_);
    HashSet<String> syms = annBuilder_.findDBObjectSymbols(foundGenes_);
    int numGenesGOID = syms.size();
    int totalGenesGOID = syms.size();
    DAGNode others = onnBuilder_.getNode("Others");
    GOTerm term = new GOTerm(others, syms, totalGenes_, 
        totalGenesGOID, numGenes,  0);
    resTerm_.add(term);
    term.print();
    return resTerm_;
  }

  public static String getLink(String org, String gene) {
    String res = "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=cerevisiae&name=" + gene + "&desc=yes&submit=Search";
    if (org == null) {
      return res;
    }
    if (org.equals("Sgd")) {
      res = "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=cerevisiae&name=" + gene + "&desc=yes&submit=Search";
    }
    if (org.equals("Mm")) {
      res = "http://source.stanford.edu/cgi-bin/source/sourceResult?criteria=" + gene + "&choice=Gene&option=Name&organism=Mm";
    }
    if (org.equals("Hs")) {
      res = "http://source.stanford.edu/cgi-bin/source/sourceResult?criteria=" + gene + "&choice=Gene&option=Name&organism=Hs";
    }
    if (org.equals("Pombie")) {
      res = "http://www.genedb.org/genedb/Dispatcher?formType=navBar&organism=pombe&name=" + gene + "&desc=yes&submit=Search";
    }
    if (org.equals("Dm")) {
      res = "http://flybase.org/cgi-bin/uniq.html?context="+ gene + "&species=Dmel&db=fbgn&caller=quicksearch";
    }
    if (org.equals("Affy")) {
      res = "https://www.affymetrix.com/LinkServlet?probeset=" + gene;
    }
    if (org.equals("Card")) {
      res = "http://www.genecards.org/cgi-bin/carddisp.pl?gene=" + gene;
    }
    return res;
  }

  public static String getHeader() {
    String res = ""
      + "<html>\n"
      + "<head>\n"
      + "<SCRIPT SRC=\"http://genepyramid.ucsd.edu/home/public/mktree.js\" LANGUAGE=\"JavaScript\"></SCRIPT>\n"
      + "<LINK REL=\"stylesheet\" HREF=\"http://genepyramid.ucsd.edu/home/public/mktree.css\">\n"
      + "</head>\n"
      + "<body>\n"
      + "<center>\n"
      + "<h1> Analysis using GO Term finder </h1>\n"
      + "<br>\n"
      + "<A href=\"#\" onClick=\"expandTree('tree1'); return false;\">Expand"
      + "All</A>&nbsp;&nbsp;&nbsp;\n"
      + "<A href=\"#\" onClick=\"collapseTree('tree1'); return false;\">Collapse"
      + "All</A>&nbsp;&nbsp;&nbsp;\n"
      + "<A href=\"#\" onClick=\"expandTreeDepth('tree1',2); return false;\">Expand"
      + "Depth 2 </A>\n"
      + "</center>\n"
      + "<ul class=\"mktree\" id=\"tree1\">\n";
    return res;
  }
  public static String getEnd() {
    String res = " </ul> </body> </html> \n";
    return res;
  }

  public void printGOTerms(BufferedWriter out, String org, Iterator<GOTerm> itr1)  throws IOException {
    while (itr1.hasNext()) {
      GOTerm term = (GOTerm) itr1.next();
      out.write("<li> <font size=+1> "+ term.getName() +"</font> " +
          term.getStats() +
          "<ul> <li> <table border=0><tr>\n");
      Iterator<String> gitr = term.getGenes();
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
  }

  public void printGOTerms(BufferedWriter out, String org, Iterator<GOTerm> itr1, HashMap<String, String> map)  throws IOException {
    while (itr1.hasNext()) {
      GOTerm term = (GOTerm) itr1.next();
      String key = map.get(term.getID());
      out.write("<li> <font size=+1> "+ term.getName() +"</font> " + key +
          term.getStats() +
          "<ul> <li> <table border=0><tr>\n");
      Iterator<String> gitr = term.getGenes();
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
  }

  public void printGOTerms(BufferedWriter out, String org, Vector<String> strSet, Vector<String> allNodes) throws IOException {
    Iterator<GOTerm> itr = getSortedGOTerms(strSet, allNodes);
    printGOTerms(out, org, itr);
  }

  public void printGOTermsTAB(BufferedWriter out, String org, GeneSet geneSet, Vector<String> allNodes) throws IOException {
    Iterator<String> itr1 = geneSet.iterator();
    HashMap<String, Vector<GOTerm> > map = new  HashMap<String, Vector<GOTerm> >();
    out.write("Name\tGOID");
    while(itr1.hasNext()) {
      String str1 = (String) itr1.next();
      out.write("\t" + str1);
      HashSet<String> set1 = geneSet.getSet(str1);
      Vector<String> strSet = new Vector<String>(set1);
      for (int k = 0; k < strSet.size(); k++) {
        System.out.println("["+strSet.get(k)+"]");
      }
      Iterator<GOTerm> itr = getSortedGOTerms(strSet, allNodes);
      while (itr.hasNext()) {
        GOTerm term = (GOTerm) itr.next();
        String id = term.getID();
        term.setAttribute("_gset", str1);
        Vector<GOTerm> res = new Vector<GOTerm>();
        if (map.containsKey(id)) {
          res = map.get(id);
        }
        res.add(term);
        map.put(id, res);
      }
    }
    out.write("\n");
    itr1 = map.keySet().iterator();
    while(itr1.hasNext()) {
      String id = (String) itr1.next();
      Vector<GOTerm> res = map.get(id);
      if (res.size() > 0) {
        GOTerm head = res.get(0);
        out.write(head.getName() + "\t" + id);
        HashMap<String, GOTerm> map1 = new  HashMap<String, GOTerm>();
        for (int i =0; i < res.size(); i++) {
          GOTerm t = res.get(i);
          String setName = (String) t.getAttribute("_gset");
          map1.put(setName, t);
        }
        Iterator<String> itr = geneSet.iterator();
        while(itr.hasNext()) {
          String str1 = (String) itr.next();
          if (map1.containsKey(str1)) {
            GOTerm t = map1.get(str1);
            out.write("\t(" + t.get_n() + "," + t.getPvalue() + ")");
          }
          else {
            out.write("\t");
          }
        }
        out.write("\n");
      }
    }
  }

  public void printGOTermsTABAlt(BufferedWriter out, String org, GeneSet geneSet, Vector<String> allNodes) throws IOException {
    Iterator<String> itr1 = geneSet.iterator();
    while(itr1.hasNext()) {
      String str1 = (String) itr1.next();
      out.write(str1 + "\n");
      System.out.println("["+str1+"]");
      HashSet<String> set1 = geneSet.getSet(str1);
      Vector<String> strSet = new Vector<String>(set1);
      /*
      for (int k = 0; k < strSet.size(); k++) {
        System.out.println("["+strSet.get(k)+"]");
      }
      */
      Iterator<GOTerm> itr = getSortedGOTerms(strSet, allNodes);
      while (itr.hasNext()) {
        GOTerm term = (GOTerm) itr.next();
        out.write("\t"+ term.getName() +"\t" + term.getStats());
        Iterator<String> gitr = term.getGenes();
        while (gitr.hasNext()) {
          String gene = (String) gitr.next();
          out.write("\t" + gene);
        }
        out.write("\n");
      }
    }
  }

  public GeneSet createGeneSet() {
    DAGGraph graph = onnBuilder_.getGraph();
    Vector<GOTerm> terms = getAllGOTerms();
    GeneSet res = new GeneSet();
    Iterator<GOTerm> itr = terms.iterator();
    while (itr.hasNext()) {
      GOTerm term = (GOTerm) itr.next();
      String id = term.getID();
      String name = term.getName();
      System.out.println(id + "\t" + name);
      Iterator<String> itr1 = term.getGenes();
      while (itr1.hasNext()) {
        String gene = (String) itr1.next();
        res.add(id, name, gene);
      }
    }
    return res;
  }

  public static void main(String args[]) throws Exception {
    System.out.println("ontology: " + args[0]);
    System.out.println("annotation: " + args[1]);
    System.out.println("output: " + args[2]);
    GOAnalysis ana = new GOAnalysis(args[0], args[1], "Mm", 0.05);
    /*
    String[] list = { "Htatip2", "Vegfc", "Ctcf", "Csnk2a2", "Ccna2", "1700009P03Rik", "Stmn1", "Cdc20", "Camk2d", "Ube1c", "Ets2", "Ccnb2", "Ing1", "Sept7", "Tubb5", "Calm2", "Cdk2", "Rnf2"};
    Vector<String> v = new Vector<String>(Arrays.asList(list));
    Vector<GOTerm> res = ana.getGOTerms(v, 7000);
    */
    GeneSet set = ana.createGeneSet();
    set.writeSetFile(args[2]);
  }

}

