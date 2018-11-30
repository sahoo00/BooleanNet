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

import java.util.Vector;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Enumeration;
import java.text.MessageFormat;

import tools.graphs.DAGNode;
import java.math.BigInteger;
import java.math.BigDecimal;

public class GOTerm {
    DAGNode term_;
    HashSet<String> genes_;
    int totalGenes_;
    int totalGenesGOID_;
    int numGenes_;
    int numGenesGOID_;

    static double[] logFact_;

    double pvalue_;
    double pvalueFdr_;
    double pvalueFwer_;

    public GOTerm(DAGNode term, HashSet<String> genes, int totalGenes, 
        int totalGenesGOID, int numGenes,  int numGenesGOID) {
        term_ = term;
        genes_ = genes;
        totalGenes_ = totalGenes;
        totalGenesGOID_ = totalGenesGOID;
        numGenes_ = numGenes;
        numGenesGOID_ = numGenesGOID;
        calculatePvalue();
    }

    public double getPvalue() { return pvalue_; }
    public double getPvalueFdr() { return pvalueFdr_; }
    public double getPvalueFwer() { return pvalueFwer_; }
    public void setPvalueFdr(double v) { pvalueFdr_ = v; }
    public void setPvalueFwer(double v) { pvalueFwer_ = v; }
    public DAGNode getDAGNode() { return term_; }
    public HashSet<String> getGeneHash() { return genes_; }
    public int get_n() { return numGenesGOID_; }

    public String getID() {
        return term_.getID();
    }

    public void setAttribute(Object key, Object val) {
        term_.setAttribute(key, val);
    }

    public Object getAttribute(Object key) {
      return term_.getAttribute(key);
    }

    public int hashCode() {
        int hash = super.hashCode();
        String res = term_.getID();
        if (res == null) {
            return hash;
        }
        return res.hashCode();
    }

    public String getName() {
        return (String) term_.getAttribute("name");
    }

    public String getNamespace() {
        return (String) term_.getAttribute("namespace");
    }

    public Iterator<String>  getGenes() { return genes_.iterator(); }

    public void calculatePvalue() {
        double pval = 0;
        // Setup logFact cache
        if (logFact_ == null || logFact_.length < (totalGenes_+2)) {
          logFact_ = new double[totalGenes_+2];
          logFact_[0] = 0; logFact_[1] = 0;
          for (int i = 2; i < (totalGenes_+2); i++) {
            logFact_[i] = logFact_[i-1] + Math.log(i);
          }
        }
        int min = numGenes_;
        if (min > totalGenesGOID_) {
            min = totalGenesGOID_;
        }
        double lo = logFact_[totalGenes_] - logFact_[numGenes_] - logFact_[totalGenes_-numGenes_];
        for (int i =numGenesGOID_; i <=min; i++) {
            double upp1 = logFact_[totalGenesGOID_] - logFact_[i] - logFact_[totalGenesGOID_-i];
            double upp2 = logFact_[totalGenes_-totalGenesGOID_] - logFact_[numGenes_-i] - logFact_[totalGenes_-totalGenesGOID_-numGenes_+i];
            pval = pval + Math.exp(upp1 + upp2 - lo);
        }
        pvalue_ = pval;
    }

    public String getStats() {
      StringBuffer res = new StringBuffer();
      res.append(" (" + numGenesGOID_ + "/" + numGenes_ + ", ");
      res.append(totalGenesGOID_ + "/" + totalGenes_ + ", ");
      MessageFormat mf = new MessageFormat(
          "{0,number,0.###E0}, {1,number,0.###E0}, {2,number,0.###E0} )");
      Object[] objs = {
        new Double(pvalue_), 
        new Double(pvalueFdr_), 
        new Double(pvalueFwer_)};
      res.append(mf.format(objs));
      // res.append(pvalue_ + ", " + pvalueFdr_ + ", " + pvalueFwer_+ " )" );
      res.append(" (" + getID()+ ")\n");
      return res.toString();
    }

    public String toString() {
        StringBuffer res = new StringBuffer();
        res.append(getName() + " (" + getID() + ", " + getNamespace() + ")\n");
        res.append("\t" + getStats());
        if (genes_ != null) {
            Iterator<String> e = genes_.iterator();
            int i = 0;
            while (e.hasNext()) {
                String gene = (String) e.next();
                i++;
                res.append("\t" + gene);
                if ( (i % 10) == 0) {
                    res.append("\n");
                }
            }
            res.append("\n");
        }
        return res.toString();
    }
    public void print() {
        System.out.println(toString());
    }
}

