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

package tools.algorithm;

import java.io.*;
import java.util.*;

/* wclique.c exact algorithm for finding one maximum-weight 
   clique in an arbitrary graph,
   10.2.2000, Patric R. J. Ostergard, 
   patric.ostergard@hut.fi */

/* compile: gcc wclique.c -o wclique -O2 */

/* usage: wclique infile */

/* infile format: see http://www.tcs.hut.fi/~pat/wclique.html */

public class MaxWeightClique {

  int Vnbr, Enbr;
  int[][] bit;
  int[] wt;
  int[] mask;

  public static int INT_SIZE = 32;
  public static int MAX_VERTEX = 2000;/* maximum number of vertices */
  public static int MAX_WEIGHT = 1000000;/* maximum weight of vertex */


  public MaxWeightClique() {
    /* initialize mask */
    mask = new int[INT_SIZE];
    mask[0] = 1;
    for(int i=1;i<INT_SIZE;i++) {
      mask[i] = mask[i-1]<<1;
    }
  }

  public boolean is_edge(int a, int b){
    return (bit[a][b/INT_SIZE]&(mask[b%INT_SIZE])) != 0;
  }

  public void readGraph(String filename) throws Exception {
    int i,j,k;
    int weight,degree,entry;
    FileReader     fr;
    BufferedReader br;
    fr = new FileReader(filename);
    br = new BufferedReader(fr);

    String record = br.readLine();
    String[] result = record.split(" ", -2); // -2 : Don't discard trailing nulls
    if (result.length < 2) {
      System.out.println("Header should contain nV nE");
      throw new Exception("Header should contain nV nE");
    }
    Vnbr = Integer.parseInt(result[0]);
    Enbr = Integer.parseInt(result[1]);

    wt = new int[Vnbr];
    bit = new int[Vnbr][];
    for(i=0;i<Vnbr;i++) {    /* empty graph table */
      bit[i] = new int[Vnbr];
      for(j=0;j<Vnbr/INT_SIZE+1;j++) {
        bit[i][j] = 0;
      }
    }
    for(i=0;i<Vnbr;i++) {
      record = br.readLine();
      if (record == null) {
        System.out.println("Line should contain Weight Degree");
        throw new Exception("Line should contain Weight Degree");
      }
      result = record.split(" ", -2); // -2 : Don't discard trailing nulls
      if (result.length < 2) {
        System.out.println("Line should contain Weight Degree");
        throw new Exception("Line should contain Weight Degree");
      }
      weight = Integer.parseInt(result[0]);
      degree = Integer.parseInt(result[1]);
      wt[i] = weight;
      for(j=2;j<(degree+2);j++) {
        entry = Integer.parseInt(result[j]);
        bit[i][entry/INT_SIZE] |= mask[entry%INT_SIZE]; /* record edge */
      }
    }
  }

  public void writeGraph(String filename) throws Exception {
    BufferedWriter out = new BufferedWriter(new FileWriter(filename));
    out.write(Vnbr + " " + Enbr + "\n");
    for(int i=0;i<Vnbr;i++) {
      int degree = 0;
      for(int j=0;j<Vnbr;j++) {
        if (is_edge(i,j)) {
          degree++;
        }
      }
      out.write(wt[i] + " "+ degree);
      for(int j=0;j<Vnbr;j++) {
        if (is_edge(i,j)) {
          out.write(" " + j);
        }
      }
      out.write("\n");
    }
    out.close();
  }

  int[] pos;    /* reordering function */
  int[] set;    /* current clique */
  int[] rec;    /* best clique so far */
  int record;   /* weight of best clique */
  int rec_level;/* # of vertices in best clique */
  int[] clique; /* table for pruning */

  public void solve() {
    int i,j,k;
    int[] nwt = new int[Vnbr];
    boolean[] used = new boolean[Vnbr];

    pos = new int[Vnbr];
    set = new int[Vnbr];
    rec = new int[Vnbr];
    clique = new int[Vnbr];

    /* order vertices */
    for(i=0;i<Vnbr;i++) {
      nwt[i] = 0;
      for(j=0;j<Vnbr;j++)
        if (is_edge(i,j)) nwt[i] += wt[j];
    }
    for(i=0;i<Vnbr;i++) {
      used[i] = false;
    }

    int count = 0;
    int p = 0;
    do {
      int min_wt = MAX_WEIGHT+1; int max_nwt = -1; 
      for(i=Vnbr-1;i>=0;i--)
        if((!used[i])&&(wt[i]<min_wt))
          min_wt = wt[i];
      for(i=Vnbr-1;i>=0;i--) {
        if(used[i]||(wt[i]>min_wt)) continue;
        if(nwt[i]>max_nwt) {
          max_nwt = nwt[i];
          p = i;
        }
      }
      pos[count++] = p;
      used[p] = true;
      for(j=0;j<Vnbr;j++)
        if((!used[j])&&(j!=p)&&(is_edge(p,j)))
          nwt[j] -= wt[p];
    } while(count<Vnbr);

    /* main routine */
    record = 0;
    int wth = 0;
    for(i=0;i<Vnbr;i++) {
      wth += wt[pos[i]];
      sub(i,pos,0,0,wth);
      clique[pos[i]] = record;
    }
    printRecord();
  }

  void printRecord() {
    System.out.print("Record: ");
    for(int i=0;i<rec_level;i++) 
      System.out.print (rec[i] + " ");
    System.out.println();
  }

  int sub(int ct,int[] table,int level,int weight,int l_weight) {
    int i,j,k;
    int best;
    int curr_weight,left_weight;
    int[] newtable = new int[Vnbr];
    int p1 = 0;
    int p2 = 0;

    if(ct<=0) { /* 0 or 1 elements left; include these */
      if(ct==0) { 
        set[level++] = table[0];
        weight += l_weight;
      }
      if(weight>record) {
        record = weight;
        rec_level = level;
        for (i=0;i<level;i++) rec[i] = set[i];
      }
      return 0;
    }
    for(i=ct;i>=0;i--) {
      if((level==0)&&(i<ct)) return 0;
      k = table[i];
      if((level>0)&&(clique[k]<=(record-weight))) return 0;  /* prune */
      set[level] = k;
      curr_weight = weight+wt[k];
      l_weight -= wt[k];
      if(l_weight<=(record-curr_weight)) return 0; /* prune */
      p1 = 0;
      p2 = 0;
      left_weight = 0;   
      while (p2<i) { 
        j = table[p2++];
        if(is_edge(j,k)) {
          newtable[p1++] = j;
          left_weight += wt[j];
        }
      }
      if(left_weight<=(record-curr_weight)) continue;
      sub(p1-1,newtable,level+1,curr_weight,left_weight);
    }
    return 0;
  }

  public static void main(String args[]) throws Exception {
    MaxWeightClique c = new MaxWeightClique();
    c.readGraph(args[1]);
    c.writeGraph(args[0]);
    c.solve();
  }
}
