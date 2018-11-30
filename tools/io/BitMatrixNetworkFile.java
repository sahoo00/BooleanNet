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

public class BitMatrixNetworkFile implements NetworkFile {

  static byte VERSION_MAJOR = 1;
  static byte VERSION_MINOR = 0;
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
  Vector<String> balanced_;

  public BitMatrixNetworkFile(String filename, int mode) throws IOException {
    file_ = new BitMatrixFile(filename, mode);
    offset_ = 0;
  }

  public BitMatrixNetworkFile(String filename) throws IOException {
    file_ = new BitMatrixFile(filename);
    offset_ = 0;
  }

  public int getType() { return NetworkFile.FILE_1_0; }

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

  public void startMatrix(int num, int numbits) throws IOException {
    if (file_.getMode() != BitMatrixFile.READMODE) {
      long ptr = file_.getFilePointer();
      file_.writeIntAt(offset_, (int)ptr);
      file_.writeInt(num);
      file_.writeInt(numbits);
      offset_ = file_.getFilePointer();
      file_.seek(ptr);
    }
    file_.startMatrix(num, numbits);
  }

  public void setBitMatrix(int a, int b, int code) throws IOException {
    file_.setBitMatrix(a, b, code);
  }

  public int readCode(int i, int j) throws IOException {
    return file_.readCode(i, j); 
  }

  public static int contraPositive(int code) {
    switch(code) {
      case 2: return 3;
      case 3: return 2;
      default: return code;
    }
  }

  public void fillLowerTriangle() throws IOException {
    if (file_.getMode() == BitMatrixFile.WRITEMODE) {
      System.out.println("Finishing lower triangular matrix");
      for (int i =0; i < file_.getNum(); i++) {
        for (int j = i+1; j < file_.getNum(); j++) {
          int code = readCode(i, j);
          setBitMatrix(j, i, contraPositive(code));
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

  public void readMatrix(Vector<String> balanced) throws IOException {
    for (int i =0; i < num_; i++) {
      for (int j = 0; j < num_; j++) {
        int code = readCode(i, j);
        if (code > 0) {
          print("Found : " + code + "\t" + i + "\t" + j);
          println("\t" + balanced.get(i) + "\t" + balanced.get(j));
        }
      }
    }
  }

  public void readMatrixHeader() throws IOException {
    file_.seek(0);
    magic_ = file_.readByte();
    major_ = file_.readByte();
    minor_ = file_.readByte();
    low_ = readList(0);
    high_ = readList(1);
    balanced_ = readList(2);
    startPtr_ = file_.readIntAt(3 + 3 * 4);
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
    printList("balanced", balanced_);
    readMatrix(balanced_);
  }

  public void writePair(int a, int b, int code, double val, String pair) throws IOException {
  }

  public static void main(String args[]) throws Exception {
    if (args.length < 1) {
      System.out.println("Arguments: <file>");
      System.exit(1);
    }
    BitMatrixNetworkFile file = new BitMatrixNetworkFile(args[0],
        BitMatrixFile.READMODE);
    file.readMatrixFile();
  }

}

