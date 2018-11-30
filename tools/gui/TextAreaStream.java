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

package tools.gui;

import java.awt.*;
import java.io.*;
import javax.swing.*;
import javax.swing.text.*;

public class TextAreaStream extends OutputStream {
  private JTextArea           _area;
  private JScrollPane         _pane;

  public TextAreaStream() {
    _area   = new JTextArea(100, 30);
    _pane   = new JScrollPane( _area );
    _area.append("Hello!\n Welcome to StepMiner log\n");

    _pane.setPreferredSize( new Dimension( 500, 500 ) );
    _area.setEditable( false );
  }

  public JComponent getOutput() {
    return _pane;
  }

  public void write( int b ) {
    System.out.write(b);
    String out;
    try {
      out = Integer.toString( b );
    }
    catch( Exception e ) {
      return;
    }
    _area.append( out );
    _area.invalidate();
  }

  public void write( byte[] buf ) throws IOException {
    System.out.write(buf);
    write( buf, 0, buf.length );
  }

  public void write( byte[] buf, int off, int len ) {
    System.out.write(buf, off, len);
    _area.append( new String( buf, off, len ) );
    _area.invalidate();
  }

  public void clear() {
    _area.setText( "" );
  }

  public void flush() throws IOException {
    System.out.flush();
  }

  public void close() throws IOException {
    System.out.close();
  }

  public void saveFile(String filename) {
    try {
      BufferedWriter out = new BufferedWriter(new FileWriter(filename));
      out.write(_area.getText());
      out.close();
    }
    catch(Exception e) {
        e.printStackTrace();
    }
  }
}

