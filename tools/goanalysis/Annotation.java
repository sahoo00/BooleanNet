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

import java.util.Enumeration;
import java.util.Arrays;
import java.util.Vector;

public class Annotation {
  String[] ann_;

  public Annotation(String[] ann) {
    ann_ = ann;
  }

  public String getDB() { return ann_[0]; }
  public String getDBObjectID() { return ann_[1]; }
  public String getDBObjectSymbol() { return ann_[2]; }
  public String getQualifier() { return ann_[3]; }
  public String getGOID() { return ann_[4]; }
  public String getDBRef() { return ann_[5]; }
  public String getEvidence() { return ann_[6]; }
  public String getWithFrom() { return ann_[7]; }
  public String getAspect() { return ann_[8]; }
  public String getDBObjectName() { return ann_[9]; }
  public String getDBObjectSynonym() { return ann_[10]; }
  public String getDBObjectType() { return ann_[11]; }
  public String getTaxon() { return ann_[12]; }
  public String getDate() { return ann_[13]; }
  public String getAssignedBy() { return ann_[14]; }
  public String getAnnotationAt(int index) { return ann_[index];}

  public Enumeration<String> getDBObjectSynonyms() { 
    String syn = ann_[10].trim();
    if (syn.equals("")) {
      Vector<String> res = new Vector<String>();
      return res.elements();
    }
    else {
      Vector<String> res = new Vector<String>(Arrays.asList(syn.split("\\|")));
      return res.elements();
    }
  }

}
