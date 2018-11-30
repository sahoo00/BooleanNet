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

public class BitMatrixSubNetworkFile implements NetworkFile {

  static byte VERSION_MAJOR = 1;
  static byte VERSION_MINOR = 2;
  static byte MAJIC = 0x55;

  BitMatrixFile file_;
  long offset_;
  byte magic_;
  byte major_;
  byte minor_;
  long startPtr_;
  int num_;
  int numBits_;

  Vector<String> low_;
  Vector<String> high_;
  Vector<String> balancedx_;
  Vector<String> balancedy_;

  public BitMatrixSubNetworkFile(String filename, int mode) throws IOException {
    file_ = new BitMatrixFile(filename, mode);
    offset_ = 0;
  }

  public BitMatrixSubNetworkFile(String filename) throws IOException {
    file_ = new BitMatrixFile(filename);
    offset_ = 0;
  }

  public int getType() { return NetworkFile.FILE_1_2; }

  public void writeHeader() throws IOException {
    file_.writeByte(MAJIC);
    file_.writeByte(VERSION_MAJOR);
    file_.writeByte(VERSION_MINOR);
    offset_ = file_.getFilePointer();
    // Index to low, high, balanced and matrix
    for (int i =0; i < 10; i++) {
      file_.writeInt(0);
    }
  }

  public void writeList(String name, Vector<String> list) throws IOException {
    long ptr = file_.getFilePointer();
    file_.writeIntAt(offset_, (int)ptr);
    offset_ = file_.getFilePointer();
    file_.seek(ptr);
    file_.writeLengthPrefixString(name);
    file_.writeInt(list.size());
    Enumeration<String> e = list.elements();
    while (e.hasMoreElements()) {
      file_.writeLengthPrefixString(e.nextElement());
    }
  }

  public void startMatrix(int num, int numBits) throws IOException {
    if (file_.getMode() != BitMatrixFile.READMODE) {
      long ptr = file_.getFilePointer();
      file_.writeIntAt(offset_, (int)ptr);
      file_.writeInt(2 * num);
      file_.writeInt(1);
      offset_ = file_.getFilePointer();
      file_.seek(ptr);
    }
    file_.startMatrix(2 * num, 1);
  }

  public void setBitMatrix(int a, int b, int code) throws IOException {
    switch(code) {
      case 1: // low -> high
        file_.setBitMatrix(2 * a, 2 * b, 1);
        break;
      case 2: // low -> low
        file_.setBitMatrix(2 * a, 2 * b + 1, 1);
        break;
      case 3: // high -> high
        file_.setBitMatrix(2 * a + 1, 2 * b, 1);
        break;
      case 4: // high -> low
        file_.setBitMatrix(2 * a + 1, 2 * b + 1, 1);
        break;
      case 5:
        file_.setBitMatrix(2 * a, 2 * b + 1, 1);
        file_.setBitMatrix(2 * a + 1, 2 * b, 1);
        break;
      case 6:
        file_.setBitMatrix(2 * a, 2 * b, 1);
        file_.setBitMatrix(2 * a + 1, 2 * b + 1, 1);
        break;
      default: break;
    }
  }

  public int readCode(int i, int j) throws IOException {
    int b0 = file_.readCode(2 * i, 2 * j); 
    int b1 = file_.readCode(2 * i, 2 * j+1); 
    int b2 = file_.readCode(2 * i+1, 2 * j); 
    int b3 = file_.readCode(2 * i+1, 2 * j+1); 
    int total = b0 + b1 + b2 + b3;
    if (total == 1) {
      if (b0 == 1) { return 1; }
      if (b1 == 1) { return 2; }
      if (b2 == 1) { return 3; }
      if (b3 == 1) { return 4; }
    }
    if (total == 2) {
      if (b1 == 1 && b2 == 1) { return 5; }
      if (b0 == 1 && b3 == 1) { return 6; }
    }
    return 0;
  }

  public static int contraPositive(int code) {
    switch(code) {
      case 2: return 3;
      case 3: return 2;
      default: return code;
    }
  }

  /**
    If filling lower triangles : balancedx should be a prefix set of balancedy
    */
  public void fillLowerTriangle() throws IOException {
    if (file_.getMode() == BitMatrixFile.WRITEMODE) {
      System.out.println("Finishing lower triangular matrix");
      for (int i =0; i < balancedx_.size(); i++) {
        System.out.println(i);
        for (int j = i+1; j < balancedx_.size(); j++) {
          int code = readCode(i, j);
          if (code > 0) {
            setBitMatrix(j, i, contraPositive(code));
          }
        }
      }
    }
  }

  public long getMatrixEnd() throws IOException {
    return file_.getMatrixStart() + balancedx_.size()*2*file_.getNumBytes();
  }

  public void fillStats() throws IOException {
    if (file_.getMode() == BitMatrixFile.WRITEMODE) {
      System.out.println("Filling network statistics");
      Vector< int[]> stats = new Vector< int[]>();
      for (int i =0; i < balancedx_.size(); i++) {
        System.out.println(i);
        int[] st = new int[10];
        for (int j = 0; j < st.length; j++) {
            st[j] = 0;
        }
        for (int j = 0; j < balancedy_.size(); j++) {
          int code = readCode(i, j);
          st[code] ++;
        }
        stats.add(st);
      }
      file_.finishMatrix();
      long ptr = getMatrixEnd();
      file_.writeIntAt(3 + 7 * 4, (int)ptr);
      file_.seek(ptr);
      for (int i =0; i < balancedx_.size(); i++) {
        int[] st = stats.get(i);
        for (int j = 0; j < st.length; j++) {
          file_.writeInt(st[j]);
        }
      }
    }
  }

  public void close() throws IOException {
    file_.finishMatrix();
    file_.close();
  }

  public void print(String a) {
    System.out.print(a);
  }

  public void println(String a) {
    System.out.println(a);
  }

  public Vector<String> readList(int num) throws IOException {
    Vector<String> res = new Vector<String>();
    long ptr = file_.readIntAt(3 + num * 4);
    file_.seek(ptr);
    String name = file_.readLengthPrefixString();
    int size = file_.readInt();
    for (int i =0; i < size; i++) {
      String l = file_.readLengthPrefixString();
      res.add(l);
    }
    return res;
  }

  public void printList(String name, Vector<String> list) throws IOException {
    int size = list.size();
    println(name + " (" + size + "):");
    for (int i =0; i < size; i++) {
      String l = list.get(i);
      print(l + ", ");
      if ((i % 5) == 0) {
        println("");
      }
    }
    println("");
  }

  public void readMatrixHeader() throws IOException {
    file_.seek(0);
    magic_ = file_.readByte();
    major_ = file_.readByte();
    minor_ = file_.readByte();
    low_ = readList(0);
    high_ = readList(1);
    balancedx_ = readList(2);
    balancedy_ = readList(3);
    startPtr_ = file_.readIntAt(3 + 4 * 4);
    num_ = file_.readInt();
    numBits_ = file_.readInt();
    file_.seek(startPtr_);
    file_.startMatrix(num_, numBits_);
  }

  public void readMatrixFile() throws IOException {
    readMatrixHeader();
    println("Magic : " + magic_);
    println("Major : " + major_);
    println("Minor : " + minor_);
    printList("low", low_);
    printList("high", high_);
    printList("balancedx", balancedx_);
    printList("balancedy", balancedy_);
    readMatrix(balancedx_, balancedy_);
  }

  public void writePair(int a, int b, int code, double val, String pair) throws IOException {
  }

  public void printStats() throws IOException {
    readMatrixHeader();
    int st = file_.readIntAt(3 + 7 * 4);
    if (st != 0 && num_ != 0) {
      /* Check end pointer */
      long endptr = getMatrixEnd();
      file_.seek(endptr);
      System.out.println("Stats : " + endptr);
      println("AID\tnorel\tlohi\tlolo\thihi\thilo\teqv\topp");
      for (int i =0; i < balancedx_.size(); i++) {
        print(balancedx_.get(i));
        for (int j =0; j < 7; j++) {
          print("\t"+file_.readInt());
        }
        for (int j =0; j < 3; j++) {
          file_.readInt();
        }
        println("");
      }
    }
  }

  public void readMatrix(Vector<String> balancedx, Vector<String> balancedy) throws IOException {
    long ptr = file_.readIntAt(3 + 4 * 4);
    int num = file_.readInt();
    int numbits = file_.readInt();
    file_.seek(ptr);
    file_.startMatrix(num, numbits);
    for (int i =0; i < balancedx.size(); i++) {
      for (int j = 0; j < balancedy.size(); j++) {
        int code = readCode(i, j);
        if (code > 0) {
          print("Found : " + code + "\t" + i + "\t" + j);
          println("\t" + balancedx.get(i) + "\t" + balancedy.get(j));
        }
      }
    }
    printStats();
  }

  public static void main(String args[]) throws Exception {
    if (args.length < 2) {
      System.out.println("Arguments: <cmd> <file>");
      System.exit(1);
    }
    if (args[0].equals("print")) {
      BitMatrixSubNetworkFile file = new BitMatrixSubNetworkFile(args[1],
          BitMatrixFile.READMODE);
      file.readMatrixFile();
    }
    if (args[0].equals("printStats")) {
      BitMatrixSubNetworkFile file = new BitMatrixSubNetworkFile(args[1],
          BitMatrixFile.READMODE);
      file.printStats();
    }
    if (args[0].equals("fill")) {
      BitMatrixFile.BLOCKSIZE = 50000;
      BitMatrixFile.CACHESIZE = 50000;
      BitMatrixSubNetworkFile file = new BitMatrixSubNetworkFile(args[1],
          BitMatrixFile.WRITEMODE);
      file.readMatrixHeader();
      file.fillLowerTriangle();
      file.close();
    }
    if (args[0].equals("fillStats")) {
      BitMatrixFile.BLOCKSIZE = 50000;
      BitMatrixFile.CACHESIZE = 50000;
      BitMatrixSubNetworkFile file = new BitMatrixSubNetworkFile(args[1],
          BitMatrixFile.WRITEMODE);
      file.readMatrixHeader();
      file.fillStats();
      file.close();
    }
  }

}

