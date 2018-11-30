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

import java.text.NumberFormat;
import javax.swing.text.NumberFormatter;
import java.beans.*;

import java.util.Properties;

public class SMSettings {

  Properties options;

  public SMSettings() {
    options = new Properties();
    createDefaults();
  }

  public Object getOption(String key) {
    return options.get(key);
  }

  public void setOption(String key, Object val) {
    options.put(key, val);
  }

  public void createDefaults() {

    options.put("Type", "TwoStep");
    options.put("numMissing", new Integer(0));
    options.put("pvalue", new Double(0.05));

    options.put("Organism", "Hs");
    options.put("geneIndex", new Integer(1));
    options.put("splitString", "\\|\\|");
    options.put("splitIndex", new Integer(0));
    options.put("onnFile", "gofull.obo");
    options.put("annFile", "gene_association.goa_human");
    options.put("goPvalue", new Double(0.05));

    options.put("width", new Integer(SMRootFrame.Frame_Width));
    options.put("height", new Integer(SMRootFrame.Frame_Height));
    options.put("xwidth.fixed", new Double(20));
    options.put("ywidth.fixed", new Double(10));
    options.put("xwidth.fill", Boolean.TRUE);
    options.put("ywidth.fill", Boolean.TRUE);
    options.put("border", new Double(5));
    options.put("contrast", new Double(2));
    options.put("color.high", Color.red);
    options.put("color.low", Color.green);
    options.put("color.zero", Color.black);
    options.put("color.missing", Color.lightGray);
  }

  public JPanel getGlobalPanel() {
    return new GlobalPanelSettings();
  }

  public JPanel getGOPanel() {
    return new GOPanelSettings();
  }

  public JPanel getMapPanel(SMRootFrame root) {
    return new MapPanelSettings(root);
  }

  public void showDialog(SMRootFrame root, String title) {
    new SMSettingsDialog(root, title, true);
  }

  class SMSettingsDialog extends JDialog implements ActionListener {
    JButton okButton;
    JButton cancelButton;
    Properties old;
    public SMSettingsDialog(SMRootFrame root, String title, boolean modal) {
      super(root.getRootFrame(), title, modal);
      old = (Properties) options.clone();
      Container pane = getContentPane();
      pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));
      pane.add(getGlobalPanel());
      pane.add(getGOPanel());
      pane.add(getMapPanel(root));
      JPanel commands = new JPanel(new FlowLayout());
      okButton = new JButton("OK");
      okButton.addActionListener(this);
      commands.add(okButton);
      cancelButton = new JButton("Cancel");
      cancelButton.addActionListener(this);
      commands.add(cancelButton);
      pane.add(commands);
      pack();
      setVisible(true);
    }
    public void actionPerformed(ActionEvent e) {
      if (e.getSource().equals(okButton)) {
        dispose();
      }
      if (e.getSource().equals(cancelButton)) {
        options = old;
        dispose();
      }
    }
  }

  class ButtonIcon implements Icon {
    int width, height;
    Color color;
    ButtonIcon(int w, int h, Color c) {
      width = w;
      height = h;
      color = c;
    }
    public int getIconHeight() {
      return height;  
    }
    public int getIconWidth() {
      return width;   
    }
    public void setColor(Color c) {
      color = c;
    }
    public void paintIcon(Component c, Graphics g, int x, int y) {
      g.setColor(color);
      g.fillRect(x, y, width, height);
      g.setColor(Color.blue);
      g.drawRect(x, y, width, height);
    }
  }

  class MapPanelSettings extends JPanel 
    implements ActionListener,DocumentListener {
      JTextField xw; JCheckBox xfill;
      JTextField yw; JCheckBox yfill;
      JTextField border;
      JTextField contrast;
      JButton high, low, zero, missing;
      SMRootFrame rootFrame;

      public MapPanelSettings(SMRootFrame root) {
        super();
        rootFrame = root;
        setBorder(BorderFactory.createTitledBorder("Heatmap"));
        setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));

        JPanel group1 = new JPanel();
        group1.setLayout(new FlowLayout());
        xw = addTextField(group1, "xwidth.fixed", "Fixed X :");
        xfill = new JCheckBox("Fill X");
        Boolean val = (Boolean) options.get("xwidth.fill");
        xfill.setSelected(val.booleanValue());
        xw.setEditable(!val.booleanValue());
        xfill.addActionListener(this);
        group1.add(xfill);
        group1.add(Box.createHorizontalStrut(50));
        yw = addTextField(group1, "ywidth.fixed", "Fixed Y :");
        yfill = new JCheckBox("Fill Y");
        val = (Boolean) options.get("ywidth.fill");
        yfill.setSelected(val.booleanValue());
        yw.setEditable(!val.booleanValue());
        yfill.addActionListener(this);
        group1.add(yfill);

        add(group1);

        JPanel group2 = new JPanel();
        group2.setLayout(new FlowLayout());
        border = addTextField(group2, "border", "Border :");
        group2.add(Box.createHorizontalStrut(50));
        contrast = addTextField(group2, "contrast", "Contrast :");

        add(group2);

        JPanel group3 = new JPanel();
        group3.setLayout(new FlowLayout());
        high = addButton(group3, "color.high", "high");
        low = addButton(group3, "color.low", "low");
        zero = addButton(group3, "color.zero", "zero");
        missing = addButton(group3, "color.missing", "missing");

        add(group3);
      }

      JButton addButton(Container pane, String prop, String label) {
        JButton button = new JButton(label);
        Color c = (Color) options.get(prop);
        button.setIcon(new ButtonIcon(10,10, c));
        button.addActionListener(this);
        pane.add(button);
        return button;
      }

      JTextField addTextField(Container pane, String prop, String label) {
        JTextField textField = new JTextField();
        textField.setColumns(3);
        textField.setText("" + options.get(prop));
        textField.addActionListener(this);
        textField.getDocument().addDocumentListener(this);
        pane.add(new JLabel(label));
        pane.add(textField);
        return textField;
      }

      public void updateNumber(ActionEvent e, double pval, String label, JTextField c) {
        if (e.getSource().equals(c)) {
          try {
            pval = Double.parseDouble(c.getText());
          }
          catch(Exception e1) {
          }
          options.put(label, new Double(pval));
        }
      }

      public void updateDNumber(DocumentEvent e, double pval, String label, JTextField c) {
        if (e.getDocument().equals(c.getDocument())) {
          try {
            pval = Double.parseDouble(c.getText());
          }
          catch(Exception e1) {
          }
          options.put(label, new Double(pval));
        }
      }

      public Color chooseColor(Color old) {
        Color c = JColorChooser.showDialog(this, "Choose color", old);
        return c;
      }

      public void updateColor(ActionEvent e, JButton b, String prop) {
        if (e.getSource().equals(b)) {
          Color old = (Color) options.get(prop);
          Color c = chooseColor(old);
          if (c != null) {
            b.setIcon(new ButtonIcon(10, 10, c));
            options.put(prop, c);
          }
        }
      }

      public void actionPerformed(ActionEvent e) {
        updateNumber(e, 20, "xwidth.fixed", xw);
        updateNumber(e, 10, "ywidth.fixed", yw);
        updateNumber(e, 5, "border", border);
        updateNumber(e, 3, "contrast", contrast);
        if (e.getSource().equals(xfill)) {
          if (xfill.isSelected()) {
            options.put("xwidth.fill", Boolean.TRUE);
            xw.setEditable(false);
          }
          else {
            options.put("xwidth.fill", Boolean.FALSE);
            xw.setEditable(true);
          }
        }
        if (e.getSource().equals(yfill)) {
          if (yfill.isSelected()) {
            options.put("ywidth.fill", Boolean.TRUE);
            yw.setEditable(false);
          }
          else {
            options.put("ywidth.fill", Boolean.FALSE);
            yw.setEditable(true);
          }
        }
        updateColor(e, high, "color.high");
        updateColor(e, low, "color.low");
        updateColor(e, zero, "color.zero");
        updateColor(e, missing, "color.missing");
        if (e.getSource().equals(contrast) ||
            e.getSource().equals(high) ||
            e.getSource().equals(low) ||
            e.getSource().equals(zero) ||
            e.getSource().equals(missing)) {
          rootFrame.updateHeatmap();
        }
      }
      /**
       *  Respond to any change in document
       */
      public void changedUpdate(DocumentEvent e) {
        updateDNumber(e, 20, "xwidth.fixed", xw);
        updateDNumber(e, 10, "ywidth.fixed", yw);
        updateDNumber(e, 5, "border", border);
        updateDNumber(e, 3, "contrast", contrast);
        if (e.getDocument().equals(contrast.getDocument()) ||
            e.getDocument().equals(border.getDocument())) {
          rootFrame.updateHeatmap();
        }
      }
      public void insertUpdate(DocumentEvent e) {
        changedUpdate(e);
      }
      public void removeUpdate(DocumentEvent e) {
        changedUpdate(e);
      }
    }
  
  class GOPanelSettings extends JPanel 
    implements ActionListener,DocumentListener {
      JComboBox org;
      JTextField geneIndex;
      JTextField splitIndex;
      JTextField splitString;
      JTextField goPvalue;

      String[] organism;

      public  GOPanelSettings() {
        super();
        setBorder(BorderFactory.createTitledBorder("GOAnalysis"));
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

        String[] orgs = { "Hs", "Mm", "Sgd", "Pombie"};
        organism = orgs;
        int index = 0;
        String selected = options.getProperty("Organism");
        System.out.println(selected);
        for (int i =0; i < orgs.length; i++) {
          if (selected.equals(orgs[i])) {
            index = i;
          }
        }
        //Create the combo box, select item at index 3.
        //Indices start at 0, so 3 specifies the TwoStep.
        org = new JComboBox(orgs);
        org.setSelectedIndex(index);
        org.addActionListener(this);
        add(new JLabel("Organism :"));
        add(org);

        add(Box.createHorizontalStrut(10));
        geneIndex = addTextField("geneIndex");
        add(Box.createHorizontalStrut(10));
        splitString = addTextField("splitString");
        add(Box.createHorizontalStrut(10));
        splitIndex = addTextField("splitIndex");
        add(Box.createHorizontalStrut(10));
        goPvalue = addTextField("goPvalue");

      }

      JTextField addTextField(String prop) {
        JTextField textField = new JTextField();
        textField.setColumns(3);
        textField.setText("" + options.get(prop));
        textField.addActionListener(this);
        textField.getDocument().addDocumentListener(this);
        add(new JLabel(prop + ":"));
        add(textField);
        return textField;
      }
      /**
       * Responds to the user choosing a new unit from the combo box.
       */
      public void actionPerformed(ActionEvent e) {
        if (e.getSource().equals(org)) {
          int i = org.getSelectedIndex();
          options.put("Organism", organism[i]);
        }
        if (e.getSource().equals(splitString)) {
          options.put("splitString", splitString.getText());
        }
        if (e.getSource().equals(geneIndex)) {
          int num = 0;
          try {
            num = Integer.parseInt(geneIndex.getText());
          }
          catch(Exception e1) {
          }
          options.put("geneIndex", new Integer(num));
        }
        if (e.getSource().equals(splitIndex)) {
          int num = 0;
          try {
            num = Integer.parseInt(splitIndex.getText());
          }
          catch(Exception e1) {
          }
          options.put("splitIndex", new Integer(num));
        }
        if (e.getSource().equals(goPvalue)) {
          double pval = 0.05;
          try {
            pval = Double.parseDouble(goPvalue.getText());
          }
          catch(Exception e1) {
          }
          options.put("goPvalue", new Double(pval));
        }
      }
      /**
       *  Respond to any change in document
       */
      public void changedUpdate(DocumentEvent e) {
        if (e.getDocument().equals(splitString.getDocument())) {
          options.put("splitString", splitString.getText());
        }
        if (e.getDocument().equals(geneIndex.getDocument())) {
          int num = 0;
          try {
            num = Integer.parseInt(geneIndex.getText());
          }
          catch(Exception e1) {
          }
          options.put("geneIndex", new Integer(num));
        }
        if (e.getDocument().equals(splitIndex.getDocument())) {
          int num = 0;
          try {
            num = Integer.parseInt(splitIndex.getText());
          }
          catch(Exception e1) {
          }
          options.put("splitIndex", new Integer(num));
        }
        if (e.getDocument().equals(goPvalue.getDocument())) {
          double pval = 0.05;
          try {
            pval = Double.parseDouble(goPvalue.getText());
          }
          catch(Exception e1) {
          }
          options.put("goPvalue", new Double(pval));
        }
      }
      public void insertUpdate(DocumentEvent e) {
        changedUpdate(e);
      }
      public void removeUpdate(DocumentEvent e) {
        changedUpdate(e);
      }
  }

  class GlobalPanelSettings extends JPanel 
    implements ActionListener,DocumentListener {
      JComboBox typeList;
      JTextField numMissing;
      JTextField pvalue;

      String[] typeStrings;

      public GlobalPanelSettings() {
        super();
        setBorder(BorderFactory.createTitledBorder("Global"));
        setLayout(new BoxLayout(this, BoxLayout.X_AXIS));

        String[] tS = { "TwoStep", "OneStep", "OneStepFdr", "TwoStepFdr"};
        typeStrings = tS;
        int index = 0;
        String selected = options.getProperty("Type");
        System.out.println(selected);
        for (int i =0; i < tS.length; i++) {
          if (selected.equals(tS[i])) {
            index = i;
          }
        }
        //Create the combo box, select item at index 3.
        //Indices start at 0, so 3 specifies the TwoStep.
        typeList = new JComboBox(typeStrings);
        typeList.setSelectedIndex(index);
        typeList.addActionListener(this);
        add(new JLabel("Type:"));
        add(typeList);

        numMissing = new JTextField();
        numMissing.setColumns(3);
        numMissing.setText("" + options.get("numMissing"));
        numMissing.addActionListener(this);
        numMissing.getDocument().addDocumentListener(this);
        add(Box.createHorizontalStrut(10));
        add(new JLabel("numMissing:"));
        add(numMissing);

        pvalue = new JTextField();
        pvalue.setColumns(4);
        pvalue.setText("" + options.get("pvalue"));
        pvalue.addActionListener(this);
        pvalue.getDocument().addDocumentListener(this);
        add(Box.createHorizontalStrut(10));
        add(new JLabel("pvalue :"));
        add(pvalue);
      }

      /**
       * Responds to the user choosing a new unit from the combo box.
       */
      public void actionPerformed(ActionEvent e) {
        if (e.getSource().equals(typeList)) {
          int i = typeList.getSelectedIndex();
          options.put("Type", typeStrings[i]);
        }
        if (e.getSource().equals(numMissing)) {
          int num = 0;
          try {
            num = Integer.parseInt(numMissing.getText());
          }
          catch(Exception e1) {
          }
          options.put("numMissing", new Integer(num));
        }
        if (e.getSource().equals(pvalue)) {
          double pval = 0.05;
          try {
            pval = Double.parseDouble(pvalue.getText());
          }
          catch(Exception e1) {
          }
          options.put("pvalue", new Double(pval));
        }
      }
      /**
       *  Respond to any change in document
       */
      public void changedUpdate(DocumentEvent e) {
        if (e.getDocument().equals(pvalue.getDocument())) {
          double pval = 0.05;
          try {
            pval = Double.parseDouble(pvalue.getText());
          }
          catch(Exception e1) {
          }
          options.put("pvalue", new Double(pval));
        }
        if (e.getDocument().equals(numMissing.getDocument())) {
          int num = 0;
          try {
            num = Integer.parseInt(numMissing.getText());
          }
          catch(Exception e1) {
          }
          options.put("numMissing", new Integer(num));
        }
      }
      public void insertUpdate(DocumentEvent e) {
        changedUpdate(e);
      }
      public void removeUpdate(DocumentEvent e) {
        changedUpdate(e);
      }
    }

}
