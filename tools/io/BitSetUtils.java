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

package tools.io;

import java.util.Arrays;
import java.util.BitSet;

public class BitSetUtils {

  public static BitSet intToBitSet(int value) {
    char[] binary = Integer.toBinaryString(value).toCharArray();
    BitSet bits = new BitSet(binary.length);
    for (int i = binary.length - 1, nbit = 0; i >= 0; i--, nbit++)
      if (binary[i] == '1') bits.set(nbit);
    return bits;
  }

  public static int bitSetToInt(BitSet bits, int size) {
    if (bits == null) return 0;
    if (bits.isEmpty()) return 0;

    if (size <= 0) throw new IllegalArgumentException("\"size\" cannot be less than or equal to zero.");

    char[] binary = new char[size];
    Arrays.fill(binary, '0');

    for (int i = 0, nbit = (size - 1); i < size; i++, nbit--)
      if (bits.get(nbit)) binary[i] = '1';

    return Integer.parseInt(new String(binary), 2);
  }

  /**
   *  Bitvector file format
   *    0 - low
   *    1 - intermediate
   *    2 - high
   *    Argument type == 0 : returns bitSet where char is 2
   *    Argument type == 1 : returns bitSet where char is not 1 and not blank
   */
  public static BitSet stringToBitSet(String str, int type) {
    BitSet res = new BitSet(str.length());
    //System.out.println(str.length());
    for (int i =0; i < str.length(); i++) {
      char c = str.charAt(i);
      res.clear(i);
      if (type == 0 && c == '2') {
        res.set(i);
      }
      if (type == 1 && !(c == '1' || c == ' ')) {
        res.set(i);
      }
    }
    //System.out.println(res.size());
    return res;
  }

}
