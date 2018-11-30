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

import java.util.HashMap;
import java.util.Vector;
import java.util.Map;

public class LHashMap<K extends Comparable<K>,V> extends HashMap<K,Vector<V> > implements Map<K,Vector<V> > {
  public LHashMap() {
    super();
  }
  public Vector<V> put(K key,V val) {
    if (key == null) {
      return null;
    }
    Vector<V> curr = (Vector<V>) get(key);
    if (curr == null) {
      curr = new Vector<V>();
    }
    curr.add(val);
    return super.put(key, curr);
  }
  public void putAllMap(HashMap<? extends K,? extends V> m) {
    if (m == null) {
      return;
    }
    java.util.Set<? extends K> set = m.keySet();
    java.util.Iterator<? extends K> itr = set.iterator();
    while (itr.hasNext()) {
      K key = (K) itr.next();
      V val = (V) m.get(key);
      Vector<V> curr = (Vector<V>) get(key);
      if (curr == null) {
        curr = new Vector<V>();
      }
      curr.add(val);
      super.put(key, curr);
    }
  }
  public void putAllUnique(LHashMap<? extends K,? extends V> m) {
    java.util.Set<? extends K> set = m.keySet();
    java.util.Iterator<? extends K> itr = set.iterator();
    while (itr.hasNext()) {
      K key = (K) itr.next();
      Vector<? extends V> val = (Vector<? extends V>) m.get(key);
      Vector<V> curr = (Vector<V>) get(key);
      if (curr == null) {
        curr = new Vector<V>();
      }
      curr.addAll(val);
      super.put(key, curr);
    }
  }
}

