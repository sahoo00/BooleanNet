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

import java.io.*;
import java.util.*;
import java.util.regex.*;
import java.net.URL;
import java.util.zip.GZIPInputStream;
import java.text.MessageFormat;

public class HomologData {

  public static int EUGENE = 1;

  String filename_;
  int type_;

  BufferedReader reader_;

  HashMap<String, EUGeneID> map_;

  public HomologData() {
    reader_ = null;
    filename_ = null;
    type_ = 0;
    map_ = new HashMap<String, EUGeneID>();
  }

  public HomologData(String file, int type) {
    reader_ = null;
    filename_ = file;
    type_ = type;
    map_ = new HashMap<String, EUGeneID>();
  }

  public void setFilename(String file, int type) {
    filename_ = file;
    type_ = type;
  }

  public void startReader() throws IOException {
    if (filename_.startsWith("http:")) {
      URL url = new URL(filename_);
      if (filename_.endsWith(".gz")) {
        reader_ = new BufferedReader(new InputStreamReader(new GZIPInputStream(url.openStream())));
      }
      else {
        reader_ = new BufferedReader(new InputStreamReader(url.openStream()));
      }
    }
    else {
      if (filename_.endsWith(".gz")) {
        reader_ = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(filename_))));
      }
      else {
        FileReader fr = new FileReader(filename_);
        reader_ = new BufferedReader(fr);
      }
    }
  }

  public void parse() throws IOException {
    if (type_ == EUGENE) {
      parseEugene();
    }
  }

  public static boolean regex_match(String pattern, String str) {
    Pattern p = Pattern.compile(pattern);
    Matcher m = p.matcher(str);
    return m.matches();
  }

  public void parseEugene() throws IOException {
    startReader();
    String record = reader_.readLine();
    EUGeneID id = null;
    for (int i =0; record != null; i++) {
      if (record.startsWith("#")) {
        record = reader_.readLine();
        continue;
      }
      boolean space = regex_match("^\\s*$", record);
      if (id == null && !space) {
        //System.out.println(record);
        id = new EUGeneID(record);
      }
      else if (space) {
        if (id != null) {
          String gene = id.getGene();
          map_.put(gene, id);
        }
        id = null;
      }
      else {
        id.add(record);
      }
      record = reader_.readLine();
    }
    close();
  }

  public void close() throws IOException {
    if (reader_ != null) {
      reader_.close();
    }
  }

  public LinkedList<String> getHomolog(String gene, String sOrg, String dOrg) {
    LinkedList<String> res = new LinkedList<String>();
    res.add(gene);
    res.add(gene.toUpperCase());
    if ( map_.containsKey(gene) ) {
      EUGeneID id = map_.get(gene);
      if (id.getOrg().equals(sOrg)) {
        LinkedList<String> tmp = id.getHomolog(gene, dOrg);
        res.addAll(tmp);
      } 
    }
    return res;
  }

  public void print() {
    Iterator<String> itr = map_.keySet().iterator();
    while (itr.hasNext()) {
      String gene = itr.next();
      System.out.println("==========------->" + gene + "<---------========");
      EUGeneID id = map_.get(gene);
      id.print();
    }
  }

  public static void main(String[] args) throws Exception {
    if (args[0].equals("print")) {
      HomologData data = new HomologData(args[1], HomologData.EUGENE);    
      data.parse();
      data.print();
    }
    if (args[0].equals("test")) {
      HomologData data = new HomologData(args[1], HomologData.EUGENE);    
      data.parse();
      LinkedList<String> list = data.getHomolog(args[2], args[3], args[4]);
      Iterator<String> itr = list.iterator();
      while (itr.hasNext()) {
        String gene = itr.next();
        System.out.println(gene);
      }
    }
  }
}

class EUGeneScore {
  public String org_;
  public String gene_;
  public int blast_score_;
  public double blast_prob_;
  public int blast_ident_;

  public EUGeneScore(String[] result) {
    parse(result);
  }

  public EUGeneScore(String record) {
    String[] result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
    parse(result);
  }

  void parse(String[] result) {
    org_ = result[2];
    gene_ = result[3];
    gene_ = gene_.replaceAll("\\s", "");
    try {
      blast_score_ = Integer.parseInt(result[4]);
      blast_prob_ = Double.parseDouble(result[5]);
      blast_ident_ = Integer.parseInt(result[6]);
    }
    catch (Exception e) {
      blast_score_ = 0;
      blast_prob_ = 1.0;
      blast_ident_ = 0;
    }
  }

  public void print() {
    System.out.println(org_ +"\t" + gene_ + "\t" + blast_score_);
  }

}

class EUGeneID {
  String id_;
  EUGeneScore score_;
  HashMap<String, LinkedList<EUGeneScore> > map_;

  public EUGeneID(String record) {
    String[] result = record.split("\\t", -2); // -2 : Don't discard trailing nulls
    id_ = result[0];
    map_ = new HashMap<String, LinkedList<EUGeneScore> >();
    EUGeneScore score = new EUGeneScore(result);
    score_ = score;
  }

  public String getGene() { return score_.gene_; }
  public String getOrg() { return score_.org_; }

  public void add(String record) {
    EUGeneScore score = new EUGeneScore(record);
    add(score);
  }

  public void add(EUGeneScore score) {
    if (!map_.containsKey(score.org_)) {
      LinkedList<EUGeneScore> list = new LinkedList<EUGeneScore>();
      map_.put(score.org_, list);
    }
    LinkedList<EUGeneScore> list = map_.get(score.org_);
    list.add(score);
  }

  public LinkedList<String> getHomolog(String gene, String dOrg) {
    LinkedList<String> res = new LinkedList<String>();
    LinkedList<EUGeneScore> list = map_.get(dOrg);
    if (list != null) {
      Iterator<EUGeneScore> itr = list.iterator();
      while (itr.hasNext()) {
        EUGeneScore s = itr.next();
        res.add(s.gene_);
        break;
      }
    }
    return res;
  }

  public void print() {
    System.out.println(id_ + ":");
    score_.print();
    Iterator<String> itr = map_.keySet().iterator();
    while (itr.hasNext()) {
      String org = itr.next();
      System.out.println(org + "-->");
      LinkedList<EUGeneScore> list = map_.get(org);
      Iterator<EUGeneScore> itr1 = list.iterator();
      while (itr1.hasNext()) {
        EUGeneScore s = itr1.next();
        s.print();
      }
    }
  }

}


