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

public class BitMatrixNetworkSimple extends BitMatrixNetworkFile {

  static byte VERSION_MAJOR = 1;
  static byte VERSION_MINOR = 1;
  static byte MAJIC = 0x55;

  public BitMatrixNetworkSimple(String filename, int mode) throws IOException {
    super(filename, mode);
  }

  public BitMatrixNetworkSimple(String filename) throws IOException {
    super(filename);
  }

  public int getType() { return NetworkFile.FILE_1_1; }

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

  public void fillLowerTriangle() throws IOException {
    if (file_.getMode() == BitMatrixFile.WRITEMODE) {
      System.out.println("Finishing lower triangular matrix");
      for (int i =0; i < file_.getNum()/2; i++) {
        System.out.println(i);
        for (int j = i+1; j < file_.getNum()/2; j++) {
          int code = readCode(i, j);
          if (code > 0) {
            setBitMatrix(j, i, contraPositive(code));
          }
        }
      }
    }
  }

  public void fillStats() throws IOException {
    if (file_.getMode() == BitMatrixFile.WRITEMODE) {
      System.out.println("Filling network statistics");
      Vector< int[]> stats = new Vector< int[]>();
      for (int i =0; i < file_.getNum()/2; i++) {
        System.out.println(i);
        int[] st = new int[10];
        for (int j = 0; j < st.length; j++) {
            st[j] = 0;
        }
        for (int j = 0; j < file_.getNum()/2; j++) {
          int code = readCode(i, j);
          st[code] ++;
        }
        stats.add(st);
      }
      file_.finishMatrix();
      long ptr = file_.getMatrixEnd();
      file_.writeIntAt(3 + 6 * 4, (int)ptr);
      file_.seek(ptr);
      for (int i =0; i < file_.getNum()/2; i++) {
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

  public void printStats() throws IOException {
    readMatrixHeader();
    int st = file_.readIntAt(3 + 6 * 4);
    if (st != 0 && num_ != 0) {
      long endptr = file_.getMatrixEnd();
      file_.seek(endptr);
      System.out.println("Stats : " + endptr);
      println("AID\tnorel\tlohi\tlolo\thihi\thilo\teqv\topp");
      for (int i =0; i < num_/2; i++) {
        print(balanced_.get(i));
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

  public void readMatrix(Vector<String> balanced) throws IOException {
    long ptr = file_.readIntAt(3 + 3 * 4);
    int num = file_.readInt();
    int numbits = file_.readInt();
    file_.seek(ptr);
    file_.startMatrix(num, numbits);
    for (int i =0; i < num/2; i++) {
      for (int j = 0; j < num/2; j++) {
        int code = readCode(i, j);
        if (code > 0) {
          print("Found : " + code + "\t" + i + "\t" + j);
          println("\t" + balanced.get(i) + "\t" + balanced.get(j));
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
      BitMatrixNetworkSimple file = new BitMatrixNetworkSimple(args[1],
          BitMatrixFile.READMODE);
      file.readMatrixFile();
    }
    if (args[0].equals("printStats")) {
      BitMatrixNetworkSimple file = new BitMatrixNetworkSimple(args[1],
          BitMatrixFile.READMODE);
      file.printStats();
    }
    if (args[0].equals("fill")) {
      BitMatrixFile.BLOCKSIZE = 50000;
      BitMatrixFile.CACHESIZE = 50000;
      BitMatrixNetworkSimple file = new BitMatrixNetworkSimple(args[1],
          BitMatrixFile.WRITEMODE);
      file.readMatrixHeader();
      file.fillLowerTriangle();
      file.close();
    }
    if (args[0].equals("fillStats")) {
      BitMatrixFile.BLOCKSIZE = 50000;
      BitMatrixFile.CACHESIZE = 50000;
      BitMatrixNetworkSimple file = new BitMatrixNetworkSimple(args[1],
          BitMatrixFile.WRITEMODE);
      file.readMatrixHeader();
      file.fillStats();
      file.close();
    }
  }

}

