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

import java.text.MessageFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


public class Step {

  int num_;
  int numSteps_;
  int[] steps_;
  double[] means_;
  double sstot_;
  double sse_;

  int geneIndex_;
  int label_;

  public Step() {
    num_ = -1;
    numSteps_ = 0;
    steps_ = null;
    means_ = null;
    sstot_ = -1;
    sse_   = -1;
    geneIndex_ = -1;
    label_ = 0;
  }

  public void setGeneIndex(int index) { geneIndex_ = index;}
  public int getGeneIndex() { return geneIndex_;}
  public int getNum() { return num_; }
  public int getNumSteps() { return numSteps_; }
  public int getStep(int segment) { return steps_[segment]; }
  public int getLabel() { return label_;}

  public double getSstot() { return sstot_; }
  public double getSse() { return sse_; }
  public double getSsr() { return sstot_-sse_; }
  public double getMean(int segment) { return means_[segment]; }
  public double getCenter() { return 0; }

  public int getDir(int segment) { 
    if (means_[segment] < means_[segment+1]) {
      return 1;
    }
    else {
      return -1;
    }
  }

  public double getStatistic() throws StepException 
    { throw new StepException("Undefined Step Statistic"); }

  public double getPvalue() throws StepException 
    { throw new StepException("Undefined Step PValue"); }

  public String toString() {
    String res = "Num=" + num_ +
        " Num Steps=" + numSteps_ +
        " Gene Index=" + geneIndex_ +
        "\nSSTOT=" + sstot_ +
        " SSE=" + sse_ +
        " Label=" + label_;
    if ( steps_ != null) {
        res = res + "\nSteps : " + steps_[0] + ", ";
        for (int i =1; i < steps_.length; i++) {
            res = res + steps_[i] ;
        }
    }
    if ( means_ != null) {
        res = res + "\nMeans : " + means_[0] + ", ";
        for (int i =1; i < means_.length; i++) {
            res = res + means_[i] ;
        }
    }
    res = res + "\n";
    return res;
  }

  public static String headString() {
    String res = "num\tnumSteps\tgeneIndex\tpvalue\tsstot\tsse\tlabel\tstep0\tstep1"
        + "\tmean0\tmean1\tmean2";
    return res;
  }

  public String tabString() throws StepException {
    String res = num_ + "\t" + numSteps_ + "\t" + geneIndex_ +
        "\t" + getPvalue() +
        "\t" + sstot_ + "\t" + sse_ + "\t" + label_;
    for (int i =0; i < 2; i++) {
      if (steps_ != null && i < steps_.length ) {
        res = res +  "\t" + steps_[i] ;
      }
      else {
        res = res +  "\t ";
      }
    }
    for (int i =0; i < 3; i++) {
      if ( means_ != null && i < means_.length) {
        res = res +  "\t" + means_[i] ;
      }
      else {
        res = res +  "\t ";
      }
    }
    return res;
  }

  public void print() {
    System.out.print(toString());
  }

  public void performCentering(double offset) {
    if (means_ != null && means_.length > 1) {
        // double offset = (means_[0] + means_[1])/2;
        for (int i =0; i < means_.length; i++) {
            means_[i] = means_[i] - offset;
        }
    }
  }

  public static String formatString(String format, double x) {
    MessageFormat mf = new MessageFormat(
        "{0,number," + format + "}");
    Object[] objs = {new Double(x)};
    return mf.format(objs);
  }

  public String getPvalueStr() throws StepException { 
    double pval = getPvalue();
    if (pval > 0.001) {
        return formatString("0.###", pval);
    }
    return formatString("0.##E0", pval);
  }

  public static HashSet<Integer> getSelectedLabels(String tag) {
    HashSet<Integer> set = new HashSet<Integer>();
    if (tag == null || tag.equals("All")) {
      for (int i =0; i < 7; i++) {
        set.add(new Integer(i));
      }
    }
    if (tag.equals("Step")) {
      for (int i =1; i < 5; i++) {
        set.add(new Integer(i));
      }
    }
    if (tag.equals("Up")) {
      set.add(new Integer(1));
    }
    if (tag.equals("Down")) {
      set.add(new Integer(2));
    }
    if (tag.equals("UpDown")) {
      set.add(new Integer(3));
    }
    if (tag.equals("DownUp")) {
      set.add(new Integer(4));
    }
    if (tag.equals("Rest")) {
      set.add(new Integer(0));
    }
    System.out.print("Selected tags(" + tag + "): ");
    Iterator<Integer> itr = set.iterator();
    while (itr.hasNext()) {
      Integer n = (Integer) itr.next();
      System.out.print("" + n + " ");
    }
    System.out.println();
    return set;
  }

  public HashSet<String> getTags(int numArrayHeader) {
    HashSet<String> set = new HashSet<String>();
    set.add("All");
    return set;
  }

  public static Integer[] orderTags(Vector<String> allSets) {
    Integer[] order = new Integer[allSets.size()];
    for (int j =0; j < order.length ; j++) {
      order[j] = new Integer(j);
    }
    class SetNameComparator implements Comparator<Integer> {
      Vector<String> names_;
      public SetNameComparator(Vector<String> names) {
        names_ = names;
      }
      int getNum(int index, String str) {
        Pattern p = Pattern.compile("^.*-(\\d+)-.*-(\\d+)");
        Matcher m = p.matcher(str);
        String numStr = null;
        if (m.matches()) {
          try {
            numStr = m.group(index+1);
            int num = Integer.parseInt(numStr);
            return num;
          }
          catch(Exception e) {
          }
        }
        else if(index == 0) {
          p = Pattern.compile("^.*-(\\d+)");
          m = p.matcher(str);
          if (m.matches()) {
            try {
              numStr = m.group(1);
              int num = Integer.parseInt(numStr);
              return num;
            }
            catch(Exception e) {
            }
          }
        }
        return -1;
      }
      int getLabel(String str) {
        Pattern p = Pattern.compile("^.*-(\\d+)-.*-(\\d+)$");
        Matcher m = p.matcher(str);
        Pattern p1 = Pattern.compile("^.*-.*-(\\d+)$");
        Matcher m1 = p1.matcher(str);
        Pattern p2 = Pattern.compile("^.*-(\\d+)-.*$");
        Matcher m2 = p2.matcher(str);
        Pattern p3 = Pattern.compile("^.*-(\\d+)$");
        Matcher m3 = p3.matcher(str);
        if (str.startsWith("UpOnly")) return 0;
        if (str.equals("Up")) return 1;
        if (str.startsWith("DownOnly")) return 2;
        if (str.equals("Down")) return 3;
        if (str.startsWith("Up") && m.matches()) return 4;
        if (str.startsWith("Up") && m1.matches()) return 5;
        if (str.startsWith("Up") && m2.matches()) return 6;
        if (str.startsWith("Down") && m.matches()) return 7;
        if (str.startsWith("Down") && m1.matches()) return 8;
        if (str.startsWith("Down") && m2.matches()) return 9;
        if (str.startsWith("Up") && m3.matches()) return 10;
        if (str.startsWith("Down") && m3.matches()) return 11;
        if (str.equals("UpDown")) return 12;
        if (str.equals("DownUp")) return 13;
        if (str.equals("Step")) return 14;
        if (str.equals("Rest")) return 15;
        if (str.equals("All")) return 16;
        return -1;
      }
      public int compare(Integer s1, Integer s2) {
        String c1 = names_.get(s1.intValue());
        String c2 = names_.get(s2.intValue());
        int l1 = getLabel(c1); 
        int l2 = getLabel(c2); 
        if (l1 == l2) {
            int n1 = getNum(0, c1);
            int n2 = getNum(0, c2);
            if (n1 == n2) {
              n2 = getNum(1, c1);
              n1 = getNum(1, c2);
            }
            if (n1 > n2) {
              return 1;
            }
            return -1;
        }
        if (l1 < l2) {
          return -1;
        }
        return 1;
      }
    };
    Arrays.sort(order, new SetNameComparator(allSets));
    return order;
  }

};

