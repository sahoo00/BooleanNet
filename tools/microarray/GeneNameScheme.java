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

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GeneNameScheme {
  int arrayIndex_;
  String org_;
  String sep_;
  int index_;

  String annotationFile_;
  String ontologyFile_;
  int numMissingPoints_;

  public GeneNameScheme() {
    arrayIndex_ = 1;
    org_ = "Mm";
    sep_ = "\\|\\|";
    index_ = -1;
    annotationFile_ = getAnnotationFile(org_);
    ontologyFile_ = getOnnFile();
    numMissingPoints_ = 0;
  }

  public GeneNameScheme(int arrayIndex, String o, String s, int i) {
    arrayIndex_ = arrayIndex;
    org_ = o;
    sep_ = s;
    index_ = i;
    annotationFile_ = getAnnotationFile(org_);
    ontologyFile_ = getOnnFile();
    numMissingPoints_ = 0;
  }

  public int getArrayIndex() { return arrayIndex_; }
  public String getOrg() { return org_; }
  public String getAnnotationFile() { return annotationFile_; }
  public String getOntologyFile() { return ontologyFile_; }
  public void setAnnotationFile(String s) { annotationFile_ = s; }
  public void setOntologyFile(String s) { ontologyFile_ = s; }
  public void setNumMissingPoints(int n) { numMissingPoints_=n; }
  public int getNumMissingPoints() { return numMissingPoints_; }

  public static String getOnnUrl() {
    return "http://www.geneontology.org/ontology/gene_ontology.obo";
  }

  public static String getOnnFile() {
    return "gene_ontology.obo";
  }

  public static String getAnnotationUrl(String org) {
    if (org.equals("Sgd")) {
      return "http://www.geneontology.org/gene-associations/gene_association.sgd.gz";
    }
    if (org.equals("Hs")) {
      return "http://www.geneontology.org/gene-associations/gene_association.goa_human.gz";
    }
    if (org.equals("Mm")) {
      return "http://www.geneontology.org/cgi-bin/downloadGOGA.pl/gene_association.mgi.gz";
    }
    if (org.equals("Pombie")) {
      return "http://www.geneontology.org/cgi-bin/downloadGOGA.pl/gene_association.GeneDB_Spombe.gz";
    }
    return null;
  }

  public static String getAnnotationFile(String org) {
    if (org.equals("Sgd")) {
      return "gene_association.sgd";
    }
    if (org.equals("Hs")) {
      return "gene_association.goa_human";
    }
    if (org.equals("Mm")) {
      return "gene_association.mgi";
    }
    if (org.equals("Pombie")) {
      return "gene_association.GeneDB_Spombe";
    }
    return null;
  }

  public String toString() {
    String res = "" +
        "Annotation File: " + annotationFile_ + "\n" +
        "Ontology File  : " + ontologyFile_ + "\n" +
        "GeneName Spec  : " + org_ + "." + arrayIndex_ + "-" + index_ +
        ", " + sep_ + "\n";
    return res;
  }

  public void print() {
    System.out.println(toString());
  }

  public String getGeneName(String str) {
    String[] res = str.split(sep_, -2);
    String gene = null;
    if (index_ >= 0 && index_ < res.length) {
      gene = res[index_];
      gene = gene.replaceAll("[\\s\"']", "");
    }
    else {
      gene = str.replaceAll("[\\s\"']", "");
    }
    return gene;
  }

  public String getUnigene(String str) {
    Pattern p = Pattern.compile("^.*("+org_+"\\.\\d+).*$");
    Matcher m = p.matcher(str);
    String gene = null;
    if (m.matches()) {
        gene = m.group(1);
    }
    return gene;
  }

  public String getGene(String str) {
    String gene = getGeneName(str);
    if (gene == null || gene.matches("^\\s*$")) {
        gene = getUnigene(str);
    }
    return gene;
  }

}

