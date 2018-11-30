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

import java.util.*;

public class GeneData implements Cloneable {

  Object[] data_;

  public GeneData(Object[] d) {
    data_ = d;
  }

  public Object[] getData() { return data_; }
  public Object getDataAt(int index) { return data_[index]; }
  public void setDataAt(int index, Object o) { data_[index] = o; }
  public int size() { return data_.length; };

  public Object clone() {
    Object[] objs = new Object[data_.length];
    for (int i = 0; i < objs.length; i++) {
      if (data_[i] != null) {
        objs[i] = data_[i];
      }
    }
    return new GeneData(objs);
  }

  public int getMissingPoints(int start, int end) {
    int count = 0;
    for (int i = start; i <= end; i++) {
      Object entry = data_[i];
      if (entry == null) {
        count++;
      }
    }
    return count;
  }

  public void performCentering(int start, double offset) {
    for (int i = start; i < data_.length; i++) {
      data_[i] = convertDouble(data_[i]);
      Double val = (Double) data_[i];
      if (val != null) {
        data_[i] = new Double(val.doubleValue() - offset);
      }
    }
  }

  public void addConstant(int start, int end, double value) {
    for (int i = start; i <= end; i++) {
      data_[i] = convertDouble(data_[i]);
      Double val = (Double) data_[i];
      if (val != null) {
        data_[i] = new Double(val.doubleValue() + value);
      }
    }
  }

  public Double convertDouble(Object entry) {
    Double res = null;
    if (entry != null && entry instanceof Double) {
      return (Double) entry;
    }
    if (entry != null && entry instanceof Boolean) {
      Boolean r = (Boolean) entry;
      if (r.booleanValue()) {
        return new Double(1.0);
      }
      else {
        return new Double(0.0);
      }
    }
    if (entry != null && entry instanceof String) {
      try {
        res = new Double(Double.parseDouble((String)entry));
      }
      catch(NumberFormatException e) {
        res = null;
      }
    }
    return res;
  }

  public void convertDouble(int start, int end) throws ArrayException {
    if (data_ == null) {
      return;
    }
    if (data_.length <= end) {
      throw new ArrayException("GeneData index out of bound error");
    }
    for (int i = start; i <= end; i++) {
        data_[i] = convertDouble(data_[i]);
    }
  }

  public void reduceLog(int start, int end) throws ArrayException {
    if (data_ == null) {
      return;
    }
    if (data_.length <= end) {
      throw new ArrayException("GeneData index out of bound error");
    }
    for (int i = start; i <= end; i++) {
      // System.out.println(data_[i]);
      Double res = convertDouble(data_[i]);
      if (res != null) {
        res = new Double(Math.log(res.doubleValue())/Math.log(2.0));
        if (res.isNaN()) {
          throw new ArrayException("Unable to take log of " + data_[i] + " at " + i);
        }
      }
      data_[i] = res;
    }
  }

  static String getString_(Object obj) {
    if (obj == null) {
      return " ";
    }
    if (obj instanceof Boolean) {
      Boolean res = (Boolean) obj;
      if (res.booleanValue()) {
        return "1";
      }
      else {
        return "0";
      }
    }
    else {
      return obj.toString();
    }
  }
  /**
   * Join the elements of array using the provided separator.
   */
  public static String join( Object[] col, String sep ) {
    StringBuffer sb = new StringBuffer();
    if ( col.length >= 1 ) {
      sb.append( getString_(col[0]) );
    }
    for ( int i = 1; i < col.length; i++ ) {
      sb.append( sep );
      sb.append( getString_(col[i]) );
    }
    return sb.toString();
  }

  public String toString() {
    return join(data_, "\t") + "\n";
  }

  public void print() {
    System.out.print(toString());
  }

  public GeneData insert(int index, Object[] a) throws ArrayException {
    if (data_ == null || data_.length <= index) {
      throw new ArrayException("GeneData index out of bound error");
    }
    Object[] obj = new Object[a.length+data_.length];
    int count = 0;
    for (int i = 0; i < index; i++) {
      obj[count++] = data_[i];
    }
    for (int i = 0; i < a.length; i++) {
      obj[count++] = a[i];
    }
    for (int i = index; i < data_.length; i++) {
      obj[count++] = data_[i];
    }
    return new GeneData(obj);
  }

  public Integer[] getSortedOrderAsc(int start, int[] perm) {
    int length = perm.length;
    Integer[] order = new Integer[length];
    for (int i =0; i < perm.length ; i++) {
        order[i] = new Integer(perm[i] + start);
    }
    class IndexComparator implements Comparator<Integer> {
      public int compare(Integer si1, Integer si2) {
        int s1 = si1.intValue();
        int s2 = si2.intValue();
        Double d1 = convertDouble(data_[s1]);
        Double d2 = convertDouble(data_[s2]);
        if (d1 == null) {
            return 1; // Swap
        }
        if (d2 == null) {
            return -1; // Don't Swap
        }
        int res = 1;  // Swap
        if (d1.doubleValue() < d2.doubleValue()) {
          res = -1;  // Don't Swap
        }
        return res;
      }
    };
    Arrays.sort(order, new IndexComparator());
    for (int i=0; i < order.length ; i++) {
        order[i] = new Integer(order[i].intValue() - start);
    }
    return order;
  }
  public Integer[] getSortedOrderDes(int start, int[] perm) {
    int length = perm.length;
    Integer[] order = new Integer[length];
    for (int i =0; i < perm.length ; i++) {
        order[i] = new Integer(perm[i] + start);
    }
    class IndexComparator implements Comparator<Integer> {
      public int compare(Integer si1, Integer si2) {
        int s1 = si1.intValue();
        int s2 = si2.intValue();
        Double d1 = convertDouble(data_[s1]);
        Double d2 = convertDouble(data_[s2]);
        if (d1 == null) {
            return -1; // Don't Swap
        }
        if (d2 == null) {
            return 1; // Swap
        }
        int res = -1;  // Don't Swap
        if (d1.doubleValue() < d2.doubleValue()) {
          res = 1;  // Swap
        }
        return res;
      }
    };
    Arrays.sort(order, new IndexComparator());
    for (int i=0; i < order.length ; i++) {
        order[i] = new Integer(order[i].intValue() - start);
    }
    return order;
  }

  public GeneData permute(int start, int[] perm) {
    Object[] obj = new Object[perm.length];
    for (int i=0; i < perm.length ; i++) {
      if (perm[i] >= 0) {
        obj[i] = data_[perm[i] + start];
      }
    }
    return new GeneData(obj);
  }

  public GeneData subset(int start, int end) {
    Object[] obj = new Object[end-start+1];
    int index = 0;
    for (int i=start; i <= end ; i++) {
      obj[index++] = data_[i];
    }
    return new GeneData(obj);
  }
  public GeneData subset(int start, HashSet<Integer> set) {
    Object[] obj = new Object[set.size()];
    int index = 0;
    for (int i=0; i < data_.length ; i++) {
      if (set.contains(new Integer(i-start))) {
        obj[index++] = data_[i];
      }
    }
    return new GeneData(obj);
  }

  public static GeneData merge(GeneData a, GeneData b) {
    Object[] obj1 = a.getData();
    Object[] obj2 = b.getData();
    Object[] obj = new Object[obj1.length + obj2.length];
    for (int i=0; i < obj1.length ; i++) {
      obj[i] = obj1[i];
    }
    for (int i=0; i < obj2.length ; i++) {
      obj[i+obj1.length] = obj2[i];
    }
    return new GeneData(obj);
  }

  public void scramble(Random random, int start, int end) {
    Object[] obj = new Object[data_.length];
    for (int j =0; j < obj.length ; j++) {
      obj[j] = data_[j];
    }
    byte[] bytes = new byte[end - start + 1];
    random.nextBytes(bytes);
    Integer[] order = new Integer[bytes.length];
    for (int j =0; j < order.length ; j++) {
      order[j] = new Integer(j);
    }
    class PermComparator implements Comparator<Integer> {
      byte[] bytes_;
      public PermComparator(byte[] b) {
        bytes_ = b;
      }
      public int compare(Integer s1, Integer s2) {
        return bytes_[s1.intValue()] - bytes_[s2.intValue()];
      }
    };
    Arrays.sort(order, new PermComparator(bytes));
    for (int j =0; j < order.length ; j++) {
      obj[j+start] = data_[order[j].intValue()+start];
    }
    data_ = obj;
  }

  public static double getMax(Object[] data, int start, int end) {
    double max = Double.MIN_VALUE;
    for (int i=start; i <= end; i++) {
      Double entry = (Double) data[i];
      if (entry != null && entry.doubleValue() > max) {
        max = entry.doubleValue();
      }
    }
    return max;
  }

  public static double getMin(Object[] data, int start, int end) {
    double min = Double.MAX_VALUE;
    for (int i=start; i <= end; i++) {
      Double entry = (Double) data[i];
      if (entry != null && entry.doubleValue() < min) {
        min = entry.doubleValue();
      }
    }
    return min;
  }

  public static int getCount(Object[] data, int start, int end) {
    int count = 0;
    for (int i=start; i <= end; i++) {
      Double entry = (Double) data[i];
      if (entry != null) {
        count++;
      }
    }
    return count;
  }

  public static double getSum(Object[] data, int start, int end) {
    double sum = 0.0;
    for (int i=start; i <= end; i++) {
      Double entry = (Double) data[i];
      if (entry != null) {
        sum = sum + entry.doubleValue();
      }
    }
    return sum;
  }

  public static double getMean(Object[] data, int start, int end) throws ArrayException {
    double sum = getSum(data, start, end);
    int count = getCount(data, start, end);
    if (count == 0) {
      throw new ArrayException("getMean: count == 0");
    }
    else {
      return (sum/count);
    }
  }

  public static double getSquareError(Object[] data, int start, int end) throws ArrayException {
    double mean = getMean(data, start, end);
    double sum = 0.0;
    for (int i=start; i <= end; i++) {
      Double entry = (Double) data[i];
      if (entry != null) {
        sum = sum + (entry.doubleValue() - mean) * (entry.doubleValue() - mean);
      }
    }
    return sum;
  }

  public static double findDistance(GeneData a, GeneData b, int start, int end) {
    int count = 0;
    double sum = 0;
    double res = Double.MAX_VALUE;
    for (int j =start; j <= end; j++) {
      Double aa = (Double) a.data_[j];
      Double bb = (Double) b.data_[j];
      if (aa != null && bb != null) {
        double tmp = aa.doubleValue() - bb.doubleValue();
        sum += tmp * tmp;
        count ++;
      }
    }
    if (count != 0) {
        res = Math.sqrt(sum);
    }
    return res;
  }

  public static int countZero(Object[] data, int start, int end) {
    int count = 0;
    for (int i=start; i <= end; i++) {
      Double entry = (Double) data[i];
      if (entry == null || entry.doubleValue() == 0) {
        count++;
      }
    }
    return count;
  }

  public void removeZero(int start, int end) {
    int count = countZero(data_, start, end);
    int len = end - start + 1;
    if (count == len) {
      for (int i=start; i <= end; i++) {
        data_[i] = null;
      }
    }
  }

  public double foldChange(int start, int end) throws ArrayException {
    convertDouble(start, end);
    double max = getMax(data_, start, end);
    double min = getMin(data_, start, end);
    return max/min;
  }

  public double foldChangeLog(int start, int end) throws ArrayException {
    convertDouble(start, end);
    double max = getMax(data_, start, end);
    double min = getMin(data_, start, end);
    return max-min;
  }

  public int countInversion(int start, int[] perm, boolean ascending) {
    int count = 0;

    for (int i =0; i < perm.length ; i++) {
        for (int j = 0; j < i; j++) {
            Double di = (Double) data_[perm[i]+start];
            Double dj = (Double) data_[perm[j]+start];
            if (di != null && dj != null) {
              if (ascending) {
                if (di.doubleValue() < dj.doubleValue()) {
                  count++;
                }
              }
              else {
                if (di.doubleValue() > dj.doubleValue()) {
                  count++;
                }
              }
            }
        }
    }

    return count;
  }

  public Double[] getVector(int start, int end) {
    Double[] obj = new Double[end-start+1];
    int index = 0;
    for (int i=start; i <= end ; i++) {
      obj[index++] = convertDouble(data_[i]);
    }
    return obj;
  }

  public void toUpperCase(int start, int end) {
    for (int i=start; i <= end ; i++) {
      if (data_[i] != null) {
        String str = getString_(data_[i]).toUpperCase();
        data_[i] = str;
      }
    }
  }

  public static GeneData getWeight(int start, int end) {
    Object[] obj = new Object[end-start+1];
    int index = 0;
    for (int i=start; i <= end ; i++) {
      obj[index++] = "1";
    }
    return new GeneData(obj);
  }

  public void addWeightCol(int index) {
    Object[] obj = new Object[data_.length+1];
    for (int i=0; i < index ; i++) {
      obj[i] = data_[i];
    }
    obj[index] = "1";
    for (int i=index+1; i < obj.length ; i++) {
      obj[i] = data_[i-1];
    }
    data_ = obj;
  }

};

