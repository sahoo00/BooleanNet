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

import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.HashMap;
import java.util.Vector;
import java.util.Map;
import java.util.Enumeration;

public class StanzaMap {
    LHashMap<String,String> map_;
    String objType_;

    public StanzaMap(LHashMap<String,String> map,  String objType) {
        map_ = map;
        objType_ = objType;
    }

    public String getObjType() { return objType_;}
    public boolean isTerm() { return objType_.equals("Term"); }
    public String getID() { 
        Vector<String> res = map_.get("id"); 
        String id = res.elementAt(0).trim();
        return id;
    }
    public String getName() { 
        Vector<String> res = map_.get("name"); 
        return res.elementAt(0).trim();
    }
    public String getNamespace() { 
        Vector<String> res = map_.get("namespace"); 
        if (res == null || res.size() <= 0) {
            return "biological_process";
        }
        return res.elementAt(0).trim();
    }
    public Enumeration<String> getParents() { 
      Vector<String> isa = map_.get("is_a"); 
      Vector<String> res = new Vector<String>();
      if (isa != null) {
        Enumeration e = isa.elements();
        while (e.hasMoreElements()) {
          String isa1 = (String) e.nextElement();
          //System.out.println(":"+isa1.trim()+":");
          res.add(isa1.trim());
        }
      }
      Vector<String> rel = map_.get("relationship"); 
      if (rel != null) {
        Enumeration e = rel.elements();
        while (e.hasMoreElements()) {
          String reln = (String) e.nextElement();
          reln = reln.trim();
          if (reln.matches("^.*part_of.*$")) {
            //System.out.println(":"+reln+":");
            Pattern p = Pattern.compile("^.*(GO\\:\\d+).*$");
            Matcher m = p.matcher(reln);
            String id = null;
            if (m.matches()) {
              id = m.group(1);
              //System.out.println(":"+id+":");
              res.add(id);
            }
          }
        }
      }
      return res.elements();
    }

    public String toString() {
        StringBuffer res = new StringBuffer();
        res.append("["+objType_+"]\n");
        java.util.Set<String> set = map_.keySet();
        java.util.Iterator<String> itr = set.iterator();
        while (itr.hasNext()) {
            String key = (String) itr.next();
            Vector<String> val = (Vector<String>) map_.get(key);
            java.util.Enumeration<String> e = val.elements();
            while (e.hasMoreElements()) {
                String v = (String) e.nextElement();
                res.append(key+": "+v); 
            }
        }
        return res.toString();
    }
    public void print() {
        System.out.println(toString());
    }
}

