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

import tools.microarray.Data;
import tools.microarray.GeneData;
import tools.microarray.FileReader.PCLFileReader;
import tools.microarray.FileWriter.PCLFileWriter;
import tools.microarray.StepMiner.StepMiner;
import tools.microarray.ArrayOrder;
import tools.microarray.Impute;
import tools.microarray.GeneNameScheme;

import java.util.*;

public class DocumentPanel extends JPanel
implements ActionListener {

  String filename;
  DocumentFrame parent;
  Data data;
  StepMiner sm;
  JTextField[] timepoints;
  JLabel[] labels;
  DocumentMap heatmap;
  JScrollPane heatmapScrollPane;
  SMSettings settings;

  public DocumentPanel(String filename, DocumentFrame p) {
    super(new BorderLayout());
    setVisible(true);
    this.filename = filename;
    this.parent = p;
    settings = p.getSettings();
    loadFile();
    addWidgets();
  }

  public void loadFile() {
    data = null;
    try {
      data = PCLFileReader.readFile(filename);
      data.convertDoubles();
    }
    catch(Exception e) {
      showDialog("Couldn't load PCL file: " + filename, "Error");
      e.printStackTrace();
      parent.closeWindow();
    }
  }

  public void saveFile(String file) {
    try {
      PCLFileWriter.writeFile(data, file, null);
    }
    catch(Exception e) {
      showDialog("Couldn't save PCL file: " + file, "Error");
      e.printStackTrace();
    }
  }

  public void saveHeatmap(String file) {
    try {
      heatmap.saveHeatmap(file);
    }
    catch(Exception e) {
      showDialog("Couldn't save Heatmap to  " + file, "Error");
      e.printStackTrace();
    }
  }

  public void runStepMiner() {
    try {
      long start = System.currentTimeMillis();
      data.setGeneNameScheme(getGeneNameScheme());
      String type = (String) settings.getOption("Type");
      Double pvalue = (Double) settings.getOption("pvalue");
      sm = new StepMiner(data);
      if (type.startsWith("OneStep")) {
        sm.setOneStepAnalysis();
      }
      if (type.startsWith("OneStepFdr")) {
        sm.setFdrAnalysis(true);
      }
      if (type.startsWith("TwoStepFdr")) {
        sm.setFdrAnalysis(true);
      }
      //sm.setStepCentering(false);
      sm.setPvalueThr(pvalue.doubleValue());
      sm.performAnalysis();
      data = sm.getStepOrderedData();
      updateHeatmap();
      // Get elapsed time in milliseconds
      long elapsedTimeMillis = System.currentTimeMillis()-start;
    
      // Get elapsed time in seconds
      float elapsedTimeSec = elapsedTimeMillis/1000F;
      String status = sm.getStats();
      status += "\nElapsed Time : " + elapsedTimeSec + " sec";
      status += "\n" + sm.getFdrStats();
      status += "\nSave the pcl file using File->Save";
      showDialog(status, "StepMiner results");
    }
    catch(Exception e) {
      e.printStackTrace();
      showDialog("StepMiner Error : Please check the settings and PCL file", "Error");
    }
  }

  public GeneNameScheme getGeneNameScheme() {
    Integer geneIndex = (Integer) settings.getOption("geneIndex");
    String org = (String) settings.getOption("Organism");
    String splitString = (String) settings.getOption("splitString");
    Integer splitIndex = (Integer) settings.getOption("splitIndex");
    Integer numMissing = (Integer) settings.getOption("numMissing");
    
    GeneNameScheme ns = new GeneNameScheme(geneIndex.intValue(), 
        org, splitString, splitIndex.intValue());
    ns.setNumMissingPoints(numMissing.intValue());
    ns.print();
    return ns;
  }

  public void runGOAnalysis(String f) {
    try {
      if (sm == null) {
        runStepMiner();
      }
      sm.setGeneNameScheme(getGeneNameScheme());
      Double goPvalue = (Double) settings.getOption("goPvalue");
      sm.performGOAnalysis(f, goPvalue);
      showDialog("Open " + f + " in a web browser", "Notice");
    }
    catch(Exception e) {
      e.printStackTrace();
      showDialog("GOAnalysis Error : Please check the settings and PCL file", "Error");
    }
  }

  public Frame getRootFrame() {
    return parent.getRootFrame();
  }

  // Show StepMiner Dialogs
  public void showDialog(String message, String title) {
    SMRootFrame.getLogger().info(title + ":" + message);
    JOptionPane.showMessageDialog(getRootFrame(), message, title,
        JOptionPane.INFORMATION_MESSAGE);
  }

  public void updateHeatmap() {
    Rectangle rect = heatmapScrollPane.getViewportBorderBounds();
    int width = (int) rect.getWidth();
    int height = (int) rect.getHeight();
    settings.setOption("width", new Integer(width));
    settings.setOption("height", new Integer(height));
    heatmap.update(data);
    heatmapScrollPane.setViewportView(heatmap);
  }

  public Double[] getTimepoints() {
    Double[] times = new Double[timepoints.length];
    for (int i =0; i < times.length; i++) {
      try {
        times[i] = new Double(Double.parseDouble(timepoints[i].getText()));
      }
      catch(Exception e) {
        times[i] = null;
      }
    }
    return times;
  }

  public void addWidgets() {
    if (data == null) {
      return;
    }
    GeneData header = data.getGeneData(0);
    Object[] objs = header.getData();
    int cols =  data.getNumColumns();
    int arrH =  data.getNumArrayHeader();

    timepoints = new JTextField[cols - arrH];
    labels = new JLabel[cols - arrH];

    JPanel group = new JPanel();
    group.setBorder(BorderFactory.createTitledBorder("Header"));
    group.setLayout(new BoxLayout(group, BoxLayout.PAGE_AXIS));
    class PaneGroup extends JPanel {
        public PaneGroup(Component a, Component b) {
          setLayout(new BoxLayout(this, BoxLayout.LINE_AXIS));
          JPanel p = new JPanel();
          Dimension d = new Dimension(100, 25);
          p.setPreferredSize(d);
          p.setMaximumSize(d);
          p.setMinimumSize(d);
          p.add(a);
          add(p);
          p = new JPanel();
          p.setLayout(new BoxLayout(p, BoxLayout.LINE_AXIS));
          d = new Dimension(400, 25);
          p.setMaximumSize(d);
          p.add(b);
          add(p);
        }
    };
    JLabel t = new JLabel("Time points");
    JLabel l = new JLabel("Header");
    group.add(new PaneGroup(t,l));
    for (int i =0; i < arrH; i++) {
        t = new JLabel("NA");
        l = new JLabel("" + objs[i]);
        group.add(new PaneGroup(t,l));
    }
    for (int i =arrH; i < cols; i++) {
        int num = (i - arrH);
        String time = "" + num;
        try {
          double val = Double.parseDouble("" + objs[i]);
          time = "" + val;
        }
        catch(Exception e) {
        }
        JTextField tf = new JTextField(time , 3);
        l = new JLabel("" + objs[i]);
        timepoints[num] = tf;
        labels[num] = l;
        group.add(new PaneGroup(tf,l));
    }
    JButton apply = new JButton("Apply");
    apply.addActionListener(this);
    group.add(apply);
    setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
    data.setTimepoints(getTimepoints());
    heatmap = new DocumentMap(data, this, settings);
    heatmapScrollPane = new JScrollPane(heatmap);
    JScrollPane groupScroll = new JScrollPane(group);
    groupScroll.setMinimumSize(new Dimension(200, 200));
    JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                    groupScroll, heatmapScrollPane);
    add(split);
  }

  /**
   * Rearrange the data according to user specification of time points.
   *    Sort the time points
   *    Permute the data according to the sorting order
   */
  public void updateData() {
    double[] times = new double[timepoints.length];
    Integer[] order = new Integer[timepoints.length];
    for (int i =0; i < times.length; i++) {
      order[i] = new Integer(i);
      try {
        times[i] = Double.parseDouble(timepoints[i].getText());
      }
      catch(Exception e) {
        times[i] = Double.MAX_VALUE;
      }
    }
    class timeComparator implements Comparator<Integer> {
      double[] times_;
      public timeComparator(double[] times) {
        times_ = times;
      }
      public int compare(Integer s1, Integer s2) {
        double c1 = times_[s1.intValue()];
        double c2 = times_[s2.intValue()];
        if (c1 > c2) {
          return 1;
        }
        return -1;
      }
    };
    Arrays.sort(order, new timeComparator(times));

    String[] n_timepoints = new String[order.length];
    String[] n_labels = new String[order.length];
    int[] perm = new int[order.length];

    for (int j =0; j < order.length ; j++) {
      perm[j] = order[j].intValue();
      n_timepoints[j] = timepoints[order[j].intValue()].getText();
      n_labels[j] = labels[order[j].intValue()].getText();
    }
    int start = data.getNumArrayHeader();
    int end = data.getNumColumns() - 1;
    for (int j =0; j < order.length ; j++) {
      timepoints[j].setText(n_timepoints[j]);
      labels[j].setText(n_labels[j]);
      if (times[order[j].intValue()] == Double.MAX_VALUE) {
        int e = j + start - 1;
        if (e < end) {
          // System.out.println(n_timepoints[j] + ":" + n_labels[j]);
          end = e;
        }
      }
    }
    data = Data.selectArraysFromData(data, perm);
    data.setRange(start + ":" + end);
    Double[] ts = getTimepoints();
    data.setTimepoints(ts);
    updateHeatmap();
  }
 
  public void actionPerformed(ActionEvent e) {
    AbstractButton source = (AbstractButton)(e.getSource());
    String event = source.getText();
    System.out.println("DocumentPanel event : " + event);
    if (event.equals("Apply")) {
      updateData();
    }
  }

}


