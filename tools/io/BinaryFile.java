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

/**
 * @author Jeff Heaton(http://www.jeffheaton.com)
 * @version 1.0
 */
public class BinaryFile
{

  /**
   * Use this constant to specify big-endian integers.
   */
  public static final short BIG_ENDIAN = 1;

  /**
   * Use this constant to specify litte-endian constants.
   */
  public static final short LITTLE_ENDIAN = 2;

  /**
   * The underlying file.
   */
  protected RandomAccessFile _file;

  /**
   * Are we in LITTLE_ENDIAN or BIG_ENDIAN mode.
   */
  protected short _endian;

  /**
   * Are we reading signed or unsigned numbers.
   */
  protected boolean _signed;

  /**
   * The constructor.  Use to specify the underlying file.
   *
   * @param f The file to read/write from/to.
   */
  public BinaryFile(RandomAccessFile f)
  {
    _file = f;
    _endian = LITTLE_ENDIAN;
    _signed = false;
  }

  /**
   * Set the endian mode for reading integers.
   *
   * @param i Specify either LITTLE_ENDIAN or BIG_ENDIAN.
   * @exception java.lang.Exception Will be thrown if this method is 
   * not passed either BinaryFile.LITTLE_ENDIAN or BinaryFile.BIG_ENDIAN.
   */
  public void setEndian(short i) throws Exception
  {
    if ((i == BIG_ENDIAN) || (i == LITTLE_ENDIAN))
      _endian = i;
    else
      throw (new Exception(
          "Must be BinaryFile.LITTLE_ENDIAN or BinaryFile.BIG_ENDIAN"));
  }

  /**
   * Returns the endian mode.  Will be either BIG_ENDIAN or LITTLE_ENDIAN.
   *
   * @return BIG_ENDIAN or LITTLE_ENDIAN to specify the current endian mode.
   */
  public int getEndian()
  {
    return _endian;
  }

  /**
   * Sets the signed or unsigned mode for integers.  true for signed, false for
   * unsigned.
   *
   * @param b True if numbers are to be read/written as signed, false if
   * unsigned.
   */
  public void setSigned(boolean b)
  {
    _signed = b;
  }

  /**
   * Returns the signed mode.
   *
   * @return Returns true for signed, false for unsigned.
   */
  public boolean getSigned()
  {
    return _signed;
  }

  /**
   * Reads a fixed length ASCII string.
   *
   * @param length How long of a string to read.
   * @return The number of bytes read.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public String readFixedString(int length) throws java.io.IOException
  {
    String rtn = "";

    for (int i = 0; i < length; i++)
      rtn += (char) _file.readByte();
    return rtn;
  }

  /**
   * Writes a fixed length ASCII string.  Will truncate the string if it does
   * not fit in the specified buffer.
   *
   * @param str The string to be written.
   * @param length The length of the area to write to.  Should be larger than
   * the length of the string being written.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void writeFixedString(String str, int length)
      throws java.io.IOException
  {
    int i;

    // trim the string back some if needed

    if (str.length() > length)
      str = str.substring(0, length);

    // write the string

    for (i = 0; i < str.length(); i++)
      _file.write(str.charAt(i));

    // buffer extra space if needed

    i = length - str.length();
    while ((i--) > 0)
      _file.write(0);
  }

  /**
   * Reads a string that stores one length byte before the string.  
   * This string can be up to 255 characters long.  Pascal stores strings this
   * way.
   *
   * @return The string that was read.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public String readLengthPrefixString() throws java.io.IOException
  {
    short len = readUnsignedByte();
    return readFixedString(len);
  }

  /**
   * Writes a string that is prefixed by a single byte that specifies the
   * length of the string.  This is how Pascal usually stores strings.
   *
   * @param str The string to be written.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void writeLengthPrefixString(String str) throws java.io.IOException
  {
    writeByte((byte) str.length());
    for (int i = 0; i < str.length(); i++)
      _file.write(str.charAt(i));
  }

  /**
   * Reads a fixed length string that is zero(NULL) terminated.  This is a type
   * of string used by C/C++.  For example char str[80].
   *
   * @param length The length of the string.

   * @return The string that was read.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public String readFixedZeroString(int length) throws java.io.IOException
  {
    String rtn = readFixedString(length);
    int i = rtn.indexOf(0);
    if (i != -1)
      rtn = rtn.substring(0, i);
    return rtn;
  }

  /**
   * Writes a fixed length string that is zero terminated.  This is the format
   * generally used by C/C++ for string storage.
   *
   * @param str The string to be written.
   * @param length The length of the buffer to receive the string.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void writeFixedZeroString(String str, int length)
      throws java.io.IOException
  {
    writeFixedString(str, length);
  }

  /**
   * Reads an unlimited length zero(null) terminated string.
   *
   * @return The string that was read.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public String readZeroString() throws java.io.IOException
  {
    String rtn = "";
    char ch;

    do
    {
      ch = (char) _file.read();
      if (ch != 0)
        rtn += ch;
    } while (ch != 0);
    return rtn;
  }

  /**
   * Writes an unlimited zero(NULL) terminated string to the file.
   *
   * @param str The string to be written.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void writeZeroString(String str) throws java.io.IOException
  {
    for (int i = 0; i < str.length(); i++)
      _file.write(str.charAt(i));
    writeByte((byte) 0);
  }

  /**
   * Internal function used to read an unsigned byte.  External classes should
   * use the readByte function.
   *
   * @return The byte, unsigned, as a short.
   * @exception java.io.IOException If an IO exception occurs.
   */
  protected short readUnsignedByte() throws java.io.IOException
  {
    return (short) (_file.readByte() & 0xff);
  }

  /**
   * Reads an 8-bit byte.  Can be signed or unsigned depending on the signed
   * property.
   *
   * @return A byte stored in a short.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public short readByte() throws java.io.IOException
  {
    if (_signed)
      return (short) _file.readByte();
    else
      return (short) _file.readUnsignedByte();
  }

  /**
   * Writes a single byte to the file.
   *
   * @param b The byte to be written.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void writeByte(short b) throws java.io.IOException
  {
    _file.write(b & 0xff);
  }

  /**
   * Reads a 16-bit word.  Can be signed or unsigned depending on the signed
   * property.  
   * Can be little or big endian depending on the endian property.
   *
   * @return A word stored in an int.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public int readWord() throws java.io.IOException
  {
    short a, b;
    int result;

    a = readUnsignedByte();
    b = readUnsignedByte();

    if (_endian == BIG_ENDIAN)
      result = ((a << 8) | b);
    else
      result = (a | (b << 8));

    if (_signed)
      if ((result & 0x8000) == 0x8000)
        result = -(0x10000 - result);

    return result;
  }

  /**
   * Write a word to the file.
   *
   * @param w The word to be written to the file.
   * @exception java.io.IOException If an IO exception occurs.
   */

  public void writeWord(int w) throws java.io.IOException
  {
    if (_endian == BIG_ENDIAN)
    {
      _file.write((w & 0xff00) >> 8);
      _file.write(w & 0xff);
    } else
    {
      _file.write(w & 0xff);
      _file.write((w & 0xff00) >> 8);
    }
  }

  /**
   * Reads a 32-bit double word.  Can be signed or unsigned 
   * depending on the signed property.  Can be little or big endian depending
   * on the endian property.
   *
   * @return A double world stored in a long.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public long readDWord() throws java.io.IOException
  {
    short a, b, c, d;
    long result;

    a = readUnsignedByte();
    b = readUnsignedByte();
    c = readUnsignedByte();
    d = readUnsignedByte();

    if (_endian == BIG_ENDIAN)
      result = ((a << 24) | (b << 16) | (c << 8) | d);
    else
      result = (a | (b << 8) | (c << 16) | (d << 24));

    if (_signed)
      if ((result & 0x80000000L) == 0x80000000L)
        result = -(0x100000000L - result);

    return result;
  }

  /**
   * Writes a double word to the file.
   *
   * @param d The double word to be written to the file.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void writeDWord(long d) throws java.io.IOException
  {
    if (_endian == BIG_ENDIAN)
    {
      _file.write((int) (d & 0xff000000) >> 24);
      _file.write((int) (d & 0xff0000) >> 16);
      _file.write((int) (d & 0xff00) >> 8);
      _file.write((int) (d & 0xff));
    } else
    {
      _file.write((int) (d & 0xff));
      _file.write((int) (d & 0xff00) >> 8);
      _file.write((int) (d & 0xff0000) >> 16);
      _file.write((int) (d & 0xff000000) >> 24);
    }
  }

  /**
   * Reads 64-bit
   *
   * @return A a long.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public long read64Word() throws java.io.IOException
  {
    long a, b, c, d, e, f, g, h;
    long result;

    a = readUnsignedByte();
    b = readUnsignedByte();
    c = readUnsignedByte();
    d = readUnsignedByte();
    e = readUnsignedByte();
    f = readUnsignedByte();
    g = readUnsignedByte();
    h = readUnsignedByte();

    if (_endian == BIG_ENDIAN)
      result = ((a << 56) | (b << 48) | (c << 40) | (d << 32) | (e << 24) | (f << 16) | (g << 8) | (h) );
    else
      result = (a | (b << 8) | (c << 16) | (d << 24) | (e << 32) | (f << 40) | (g << 48) | (h << 56) );

    return result;
  }

  /**
   * Writes 64bit to the file.
   *
   * @param d The 64bit to be written to the file.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void write64Word(long d) throws java.io.IOException
  {
    if (_endian == BIG_ENDIAN)
    {
      _file.write((int) ((d >> 56) & 0xff));
      _file.write((int) ((d >> 48) & 0xff));
      _file.write((int) ((d >> 40) & 0xff));
      _file.write((int) ((d >> 32) & 0xff));
      _file.write((int) ((d >> 24) & 0xff));
      _file.write((int) ((d >> 16) & 0xff));
      _file.write((int) ((d >>  8) & 0xff));
      _file.write((int) ((d      ) & 0xff));
    } else
    {
      _file.write((int) ((d      ) & 0xff));
      _file.write((int) ((d >>  8) & 0xff));
      _file.write((int) ((d >> 16) & 0xff));
      _file.write((int) ((d >> 24) & 0xff));
      _file.write((int) ((d >> 32) & 0xff));
      _file.write((int) ((d >> 40) & 0xff));
      _file.write((int) ((d >> 48) & 0xff));
      _file.write((int) ((d >> 56) & 0xff));
    }
  }

  /**
   * Allows the file to be aligned to a specified byte boundary.  
   * For example, if a 4(double word) is specified, the file pointer will be 
   * moved to the next double word boundary.
   *
   * @param a The byte-boundary to align to.
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void align(int a) throws java.io.IOException
  {
    if ((_file.getFilePointer() % a) > 0)
    {
      long pos = _file.getFilePointer() / a;
      _file.seek((pos + 1) * a);
    }
  }

  /**
   * Read len bytes into the buffer
   *
   * @param b Buffer
   * @param off offset into the buffer
   * @param len Length of bytes to be read
   * @exception java.io.IOException If an IO exception occurs.
   */
  public int readBuffer(byte[] b, int off, int len) throws java.io.IOException
  {
    return _file.read(b, off, len);
  }

  /**
   * Write len bytes from the buffer
   *
   * @param b Buffer
   * @param off offset into the buffer
   * @param len Length of bytes to be written
   * @exception java.io.IOException If an IO exception occurs.
   */
  public void writeBuffer(byte[] b, int off, int len) throws java.io.IOException
  {
    _file.write(b, off, len);
  }

}
