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
import java.awt.image.*;

import tools.microarray.Data;
import tools.microarray.GeneData;
import tools.microarray.FileReader.PCLFileReader;
import tools.microarray.FileWriter.PCLFileWriter;
import tools.microarray.StepMiner.StepMiner;
import tools.microarray.ArrayOrder;
import tools.microarray.Impute;
import tools.microarray.GeneNameScheme;

import java.util.*;
import java.io.*;

public class DocumentMap extends JPanel
implements ActionListener {

  DocumentPanel parent;
  Data data;
  DocumentImage heatmap;
  BufferedImage image;
  SMSettings settings;

  public DocumentMap(Data d, DocumentPanel p, SMSettings s) {
    super(new BorderLayout());
    setVisible(true);
    this.data = d;
    this.parent = p;
    settings = s;
    heatmap = new DocumentImage(data);
    heatmap.update(settings);
    //heatmap.setMaximumHeight(1000);
    image = heatmap.getImage();
  }

  public void update(Data d) {
    this.data = d;
    heatmap = new DocumentImage(d);
    heatmap.update(settings);
    //heatmap.setMaximumHeight(500);
    image = heatmap.getImage();
    repaint();
  }

  public void saveHeatmap(String file) throws IOException {
    heatmap.save(file);
  }

  public void paintComponent(Graphics g) {
    super.paintComponent(g);
    paintHeatmap(g);
  }

  public Dimension getMaximumSize() {
    return new Dimension(heatmap.getWidth(), heatmap.getHeight());
  }

  public Dimension getPreferredSize() {
    return new Dimension(heatmap.getWidth(), heatmap.getHeight());
  }

  public Dimension getMinimumSize() {
    return new Dimension(heatmap.getWidth(), heatmap.getHeight());
  }

  public void paintHeatmap(Graphics g) {
    g.drawImage(image, 0, 0, null);
  }

  public void actionPerformed(ActionEvent e) {
    AbstractButton source = (AbstractButton)(e.getSource());
    String event = source.getText();
    System.out.println("DocumentMap event : " + event);
  }

}


