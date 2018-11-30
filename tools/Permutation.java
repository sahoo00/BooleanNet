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

package tools;

import java.util.*;

public class Permutation {

  Random random_;
  byte[] bytes_;

  public Permutation(int seed) {
    random_ = new Random(seed);
  }

  public int[] getRandomPermutation(int num) {
    bytes_ = new byte[num];
    random_.nextBytes(bytes_);
    Integer[] order = new Integer[bytes_.length];
    for (int j =0; j < order.length ; j++) {
      order[j] = new Integer(j);
    }
    class PermComparator implements Comparator<Integer> {
      public int compare(Integer s1, Integer s2) {
        return bytes_[s1.intValue()] - bytes_[s2.intValue()];
      }
    };
    Arrays.sort(order, new PermComparator());
    int[] perm = new int[bytes_.length];
    for (int j =0; j < order.length ; j++) {
      perm[j] = order[j].intValue();
    }
    return perm;
  }

  public static int[] getNullPermutation(int num) {
    int[] perm = new int[num];
    for (int j =0; j < perm.length ; j++) {
      perm[j] = j;
    }
    return perm;
  }

  public static Double[] permute(Double[] x, int[] perm) {
    Double[] res = new Double[perm.length];
    for (int i =0; i < perm.length; i++) {
      res[i] = x[perm[i]];
    }
    return res;
  }

  public static byte[] permute(byte[] x, int[] perm) {
    byte[] res = new byte[perm.length];
    for (int i =0; i < perm.length; i++) {
      res[i] = x[perm[i]];
    }
    return res;
  }

}

