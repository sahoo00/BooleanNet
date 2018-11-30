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

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;

public class DocumentFrame extends JInternalFrame 
implements InternalFrameListener {

  SMRootFrame root;
  String filename, name;
  DocumentPanel doc;

  public DocumentFrame(String filename, String name, SMRootFrame r) {
    // System.out.println(filename + ":" + name);
    super(name,
        true, //resizable
        true, //closable
        true, //maximizable
        true);//iconifiable
    setSize(300,300);
    addInternalFrameListener(this);
    setVisible(true);
    this.root = r;
    this.filename = filename;
    this.name = name;
    loadDocument();
  }

  public void loadDocument() {
    doc = new DocumentPanel(filename, this);
    Container pane = getContentPane();
    pane.add(doc);
    pack();
  }

  public void saveFile(String f) {
    if (doc != null) {
      doc.saveFile(f);
    }
  }

  public void saveHeatmap(String f) {
    if (doc != null) {
      doc.saveHeatmap(f);
    }
  }

  public void runStepMiner() {
    if (doc != null) {
      doc.runStepMiner();
    }
  }

  public void runGOAnalysis(String f) {
    if (doc != null) {
      doc.runGOAnalysis(f);
    }
  }

  public void maximize() {
    try {
      setMaximum(true);
    }
    catch(Exception e) {
    }
  }

  public SMSettings getSettings() {
    return root.getSettings();
  }

  public Frame getRootFrame() {
    return root.getRootFrame();
  }

  public void updateHeatmap() {
    doc.updateHeatmap();
  }

  public void closeWindow() {
    root.removeIntFrameName(name);
    dispose();
  }

  public void   internalFrameActivated(InternalFrameEvent e) {}
  public void   internalFrameClosed(InternalFrameEvent e){}
  public void   internalFrameClosing(InternalFrameEvent e) {
    closeWindow();
  }
  public void   internalFrameDeactivated(InternalFrameEvent e){}
  public void   internalFrameDeiconified(InternalFrameEvent e){}
  public void   internalFrameIconified(InternalFrameEvent e){}
  public void   internalFrameOpened(InternalFrameEvent e){}
}

