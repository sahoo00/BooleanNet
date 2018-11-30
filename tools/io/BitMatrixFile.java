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

public class BitMatrixFile {

  String filename_;
  RandomAccessFile file_;
  BinaryFile bfile_;

  long matrix_start_;
  int num_;
  int num_bits_;
  int num_bytes_;

  public static int BLOCKSIZE = 5000;
  public static int CACHESIZE = 5000;
  public static int READMODE = 0;
  public static int WRITEMODE = 1;

  HashMap<Integer, byte[]> cache_;
  int mode_;

  public BitMatrixFile(String filename, int mode) throws IOException {
    filename_ = filename;
    if (mode == WRITEMODE) {
      file_= new RandomAccessFile(filename,"rw");
    }
    else {
      file_= new RandomAccessFile(filename,"r");
    }
    bfile_= new BinaryFile(file_);
    cache_ = new HashMap<Integer, byte[]>();
    mode_ = mode;
  }

  public BitMatrixFile(String filename) throws IOException {
    filename_ = filename;
    file_= new RandomAccessFile(filename,"rw");
    bfile_= new BinaryFile(file_);
    cache_ = new HashMap<Integer, byte[]>();
    mode_ = WRITEMODE;
  }

  public BinaryFile getBinaryFile() { return bfile_; }
  public RandomAccessFile getRandomAccessFile() { return file_; }
  public int getNum() { return num_; }
  public int getNumBits() { return num_bits_; }
  public int getNumBytes() { return num_bytes_; }
  public int getMode() { return mode_; }
  public long getMatrixStart() { return matrix_start_;}
  public long getMatrixEnd() { return matrix_start_+num_*num_bytes_;}

  public int readInt() throws IOException {
    long dword = bfile_.readDWord();
    return (int) dword;
  }

  public int readShort() throws IOException {
    int word = bfile_.readWord();
    return word;
  }

  public byte readByte() throws IOException {
    short b = bfile_.readByte();
    return (byte) b;
  }

  public void seek(long ptr) throws IOException {
    file_.seek(ptr);
  }

  public long getFilePointer() throws IOException {
    return file_.getFilePointer();
  }

  public int readIntAt(long ptr) throws IOException {
    seek(ptr);
    return readInt();
  }

  public int readShortAt(long ptr) throws IOException {
    seek(ptr);
    return readShort();
  }

  public byte readByteAt(long ptr) throws IOException {
    seek(ptr);
    return readByte();
  }

  public String readString(int length)  throws IOException {
    String res =  bfile_.readFixedString(length);
    return res;
  }

  public String readStringAt(long ptr, int length)  throws IOException {
    seek(ptr);
    String res =  bfile_.readFixedString(length);
    return res;
  }

  public void writeInt(int a) throws IOException {
    bfile_.writeDWord((long)a);
  }

  public void writeShort(short a) throws IOException {
    bfile_.writeWord((int)a);
  }

  public void writeByte(byte a) throws IOException {
    bfile_.writeByte((short)a);
  }

  public void writeIntAt(long ptr, int a) throws IOException {
    seek(ptr);
    bfile_.writeDWord((long)a);
  }

  public void writeShortAt(long ptr, short a) throws IOException {
    seek(ptr);
    bfile_.writeWord((int)a);
  }

  public void writeByteAt(long ptr, byte a) throws IOException {
    seek(ptr);
    bfile_.writeByte((short)a);
  }

  public void writeString(String str)  throws IOException {
    bfile_.writeFixedString(str, str.length());
  }

  public void writeStringAt(long ptr, String str)  throws IOException {
    seek(ptr);
    bfile_.writeFixedString(str, str.length());
  }

  public void writeLengthPrefixString(String str)  throws IOException {
    writeInt(str.length());
    bfile_.writeFixedString(str, str.length());
  }

  public String readLengthPrefixString()  throws IOException {
    int length = readInt();
    String res =  bfile_.readFixedString(length);
    return res;
  }

  public void close() throws IOException {
    file_.close();
  }

  public void startMatrix(int num, int numbits) throws IOException {
    matrix_start_ = getFilePointer();
    num_ = num;
    num_bits_ = numbits;
    num_bytes_ = (num_ * num_bits_)/8 + 1;
    System.out.println("StartMatrix : [" + num_ + ", " + num_bits_ +", " + num_bytes_ + "]");
  }

  public int hashBlock(int a, int b) {
    int hash = a * num_bytes_ + (b * num_bits_/8);
    return hash / BLOCKSIZE;
  }

  public void flushCache() throws IOException {
    if (mode_ == WRITEMODE) {
      Iterator<Integer> itr = cache_.keySet().iterator();
      while (itr.hasNext()) {
        Integer a = itr.next();
        byte[] buffer = cache_.get(a);
        long pos = matrix_start_ + a.intValue() * BLOCKSIZE;
        file_.seek(pos);
        file_.write(buffer, 0, buffer.length);
      }
    }
    cache_.clear();
    System.out.println("Cache Clear");
  }

  public void loadBlock(Integer hash) throws IOException {
    // If cache size full - flush them to the file
    if (cache_.size() >= CACHESIZE) {
      flushCache();
    }
    if (!cache_.containsKey(hash)) {
      long pos = matrix_start_ + hash * BLOCKSIZE;
      byte[] buffer = new byte[BLOCKSIZE];
      file_.seek(pos);
      int num = file_.read(buffer, 0, BLOCKSIZE);
      if (num < 0) {
        num = 0;
      }
      for (int i = num; i < BLOCKSIZE; i++) {
        buffer[i] = 0;
      }
      cache_.put(hash, buffer);
    }
  }

  public void setBitMatrix(int a, int b, int code) throws IOException {
    Integer hash = new Integer(hashBlock(a, b));
    byte[] buffer;
    loadBlock(hash);
    buffer = cache_.get(hash);
    int byte_offset = (a * num_bytes_ + (b * num_bits_/8)) % BLOCKSIZE;
    int bit_offset = (b * num_bits_) % 8;
    int mask = (1 << num_bits_) - 1;
    if ((bit_offset + num_bits_) > 8) {
      if ((byte_offset+1) < BLOCKSIZE) {
        //System.out.println("Case 1");
        int val = (buffer[byte_offset]&0xff)|(buffer[byte_offset+1] << 8);
        val = ((val & ~(mask << bit_offset)) | ((code & mask) << bit_offset));
        buffer[byte_offset] = (byte)(val & 0xff);
        buffer[byte_offset+1] = (byte)(val >> 8);
      }
      else {
        //System.out.println("Case 2");
        int val = buffer[byte_offset] & 0xff;
        loadBlock(hash+1);
        buffer = cache_.get(hash+1);
        val = val | (buffer[0] << 8);
        val = ((val & ~(mask << bit_offset)) | ((code & mask) << bit_offset));
        buffer[0] = (byte)(val >> 8);
        loadBlock(hash);
        buffer = cache_.get(hash);
        buffer[byte_offset] = (byte)(val & 0xff);
      }
    }
    else {
      int val = buffer[byte_offset] & 0xff;
      val =((val & ~(mask << bit_offset))|((code & mask) << bit_offset));
      buffer[byte_offset] = (byte) (val & 0xff);
    }
    if (false && code > 0) {
      long pos = matrix_start_ + hash.intValue() * BLOCKSIZE;
      System.out.println(code + "\t" + a + "\t" + b + "\t" + pos + "\t" + mask);
      System.out.print(byte_offset + "\t" + bit_offset);
      System.out.print("\t0x" + Integer.toHexString(buffer[byte_offset]&0xff));
      if ((byte_offset+1) < BLOCKSIZE) {
        System.out.print("\t0x" + Integer.toHexString(buffer[byte_offset+1]&0xff));
      }
      System.out.println();
    }
  }

  public int readCode(int a, int b) throws IOException {
    Integer hash = new Integer(hashBlock(a, b));
    byte[] buffer;
    loadBlock(hash);
    buffer = cache_.get(hash);
    int byte_offset = (a * num_bytes_ + (b * num_bits_/8)) % BLOCKSIZE;
    int bit_offset = (b * num_bits_) % 8;
    int mask = (1 << num_bits_) - 1;
    int code = 0;
    if ((bit_offset + num_bits_) > 8) {
      if ((byte_offset+1) < BLOCKSIZE) {
        //System.out.println("Case 1");
        int val = (buffer[byte_offset]&0xff)|(buffer[byte_offset+1] << 8);
        code = (val >> bit_offset) & mask;
        if (false && code > 0) {
          System.out.println("Val 0x" + Integer.toHexString(val&0xffff));
        }
      }
      else {
        //System.out.println("Case 2");
        int val = buffer[byte_offset] & 0xff;
        loadBlock(hash+1);
        buffer = cache_.get(hash+1);
        val = val | (buffer[0] << 8);
        code = (val >> bit_offset) & mask;
      }
    }
    else {
      int val = buffer[byte_offset] & 0xff;
      code = (val >> bit_offset) & mask;
    }
    if (false && code > 0) {
      long pos = matrix_start_ + hash.intValue() * BLOCKSIZE;
      System.out.println(code + "\t" + a + "\t" + b + "\t" + pos + "\t" + mask);
      System.out.print(byte_offset + "\t" + bit_offset);
      System.out.print("\t0x" + Integer.toHexString(buffer[byte_offset]&0xff));
      if ((byte_offset+1) < BLOCKSIZE) {
        System.out.print("\t0x" + Integer.toHexString(buffer[byte_offset+1]&0xff));
      }
      System.out.println();
    }
    return code;
  }

  public void finishMatrix() throws IOException {
    flushCache();
  }
}

