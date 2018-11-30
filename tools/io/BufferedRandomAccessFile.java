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

public class BufferedRandomAccessFile extends RandomAccessFile {

  private static int BUF_SIZE = 8096;

  byte buffer[];
  int buf_end = 0;
  int buf_pos = 0;
  long real_pos = 0;

  public BufferedRandomAccessFile(String filename, String mode) 
    throws IOException{
      super(filename,mode);
      invalidate();
      buffer = new byte[BUF_SIZE];
    }

  public BufferedRandomAccessFile(String filename, String mode, int bufsize) 
    throws IOException{
      super(filename,mode);
      invalidate();
      BUF_SIZE = bufsize;
      buffer = new byte[BUF_SIZE];    
    }

  public final int read() throws IOException{
    if(buf_pos >= buf_end) {
      if(fillBuffer() < 0)
        return -1;
    }
    if(buf_end == 0) {
      return -1;
    } else {
      return buffer[buf_pos++];
    }
  }
  private int fillBuffer() throws IOException {
    int n = super.read(buffer, 0, BUF_SIZE);
    if(n >= 0) {
      real_pos +=n;
      buf_end = n;
      buf_pos = 0;
    }
    return n;
  }
  private void invalidate() throws IOException {
    buf_end = 0;
    buf_pos = 0;
    real_pos = super.getFilePointer();
  }

  public int read(byte b[], int off, int len) throws IOException {
    int leftover = buf_end - buf_pos;
    if(len <= leftover) {
      System.arraycopy(buffer, buf_pos, b, off, len);
      buf_pos += len;
      return len;
    }
    for(int i = 0; i < len; i++) {
      int c = this.read();
      if(c != -1)
        b[off+i] = (byte)c;
      else {
        return i;
      }
    }
    return len;
  }
  public long getFilePointer() throws IOException{
    long l = real_pos;
    return (l - buf_end + buf_pos) ;
  }
  public void seek(long pos) throws IOException {
    int n = (int)(real_pos - pos);
    if(n >= 0 && n <= buf_end) {
      buf_pos = buf_end - n;
    } else {
      super.seek(pos);
      invalidate();
    }
  }
  /**
   * return a next line in String 
   */
  public final String getNextLine() throws IOException {
    String str = null;
    if(buf_end-buf_pos <= 0) {
      if(fillBuffer() < 0) {
        return null;
        //throw new IOException("error in filling buffer! " + 
        //    buf_pos + "-" + buf_end + ":" + real_pos);
      }
    }
    int lineend = -1;
    for(int i = buf_pos; i < buf_end; i++) {
      if(buffer[i] == '\n') {
        lineend = i;
        break;
      }
    }
    if(lineend < 0) {
      StringBuffer input = new StringBuffer(256);
      int c;
      while (((c = read()) != -1) && (c != '\n')) {
        input.append((char)c);
      }
      if ((c == -1) && (input.length() == 0)) {
        return null;
      }
      return input.toString();
    }
    if(lineend > 0 && buffer[lineend-1] == '\r')
      str = new String(buffer, buf_pos, lineend - buf_pos -1);
    else str = new String(buffer, buf_pos, lineend - buf_pos);
    buf_pos = lineend +1;
    return str;
  }

}

