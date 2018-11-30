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

import java.net.URL;
import java.io.File;
import java.util.HashMap;
import java.util.logging.*;

public class SMRootFrame extends WindowAdapter
implements Runnable,ActionListener,ItemListener {

  SMRootFrame next; // Doubly Linked list of Root Frames
  SMRootFrame prev; // Doubly Linked list of Root Frames
  JFrame root;
  JDesktopPane desktop;
  JFileChooser fc; // File chooser
  String newline = "\n";
  HashMap<String, Object> docMap; // Retrieve the DocumentPanel from Name
                                  //    Give a unique name to each document
  boolean loaded;
  SMSettings settings;

  static String F_OPEN = "Open";
  static String F_NEW = "New Window";
  static String F_SAVE = "Save";
  static String F_CLOSE = "Close";
  static String F_EXIT = "Exit";
  static String V_LOG = "Log";
  static String S_SETTINGS = "Change Settings";
  static String E_POSTSCRIPT = "Export to Postscript";
  static String E_IMAGE = "Export to Image";
  static String H_ABOUT = "About";
  static String A_STEPMINER = "StepMiner";
  static String A_GOANA = "GOAnalysis";

  private static Logger logger = Logger.getLogger("tools.gui");
  private static TextAreaStream log;

  public static int Frame_Width = 740;
  public static int Frame_Height = 480;

  public static PNGFileFilter pngFilter = new PNGFileFilter();
  public static TXTFileFilter txtFilter = new TXTFileFilter();
  public static PCLFileFilter pclFilter = new PCLFileFilter();
  public static HTMLFileFilter htmlFilter = new HTMLFileFilter();

  public static String VERSION = "1.0 Alpha";

  static {
    log = new TextAreaStream();
    TextAreaStreamHandler sh = new TextAreaStreamHandler(log, new SimpleFormatter());
    sh.setLevel(Level.ALL);
    logger.addHandler(sh);
    logger.setLevel(Level.ALL);
    logger.fine("Logger Started");
  }

  public static Logger getLogger() {
    return logger;
  }

  public SMRootFrame() {
    loaded = false;
    next = null;
    prev = null;
    logger.info("New Root Frame");
    GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
    GraphicsDevice gsd[] = ge.getScreenDevices();
    GraphicsConfiguration gc[] = gsd[0].getConfigurations();
    root = new JFrame("StepMiner", gc[0]);
    //root.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    root.setSize(Frame_Width, Frame_Height);

    // Displaying the window in the center
    Rectangle gcBounds = gc[0].getBounds();
    // System.out.println(gcBounds);
    logger.fine(gcBounds.toString());
    int xoffs = (gcBounds.width - Frame_Width)/2;
    int yoffs = (gcBounds.height - Frame_Height)/2;
    if (xoffs < 0) xoffs = 0;
    if (yoffs < 0) yoffs = 0;
    root.setLocation(xoffs, yoffs);

    EventQueue.invokeLater(this);
  }

  // Create a new window
  //  Maintains a doubly linked list of root frames
  public SMRootFrame(SMRootFrame p) {
    loaded = false;
    logger.info("New Root Frame");
    prev = p;
    next = p.next;
    p.next = this;
    if (next != null) {
      next.prev = this;
    }
    root = new JFrame("StepMiner", p.getGraphicsConfiguration());
    //root.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    root.setSize(Frame_Width, Frame_Height);

    // Displaying the window in the center
    Rectangle gcBounds = p.getBounds();
    // System.out.println(gcBounds);
    int xoffs = gcBounds.x + 30;
    int yoffs = gcBounds.y + 30;
    root.setLocation(xoffs, yoffs);

    EventQueue.invokeLater(this);
  }

  public boolean isLoaded() {
    return loaded;
  }

  public synchronized void afterReady() throws InterruptedException {
    if (!isLoaded()) {
      wait();
    }
  }

  // Get the number of root frames currently present
  public int getNumRootFrames() {
    SMRootFrame top = this;
    while (top.prev != null) {
      top = top.prev;
    }
    int count = 0;
    while (top != null) {
      top = top.next;
      count ++;
    }
    return count;
  }

  GraphicsConfiguration getGraphicsConfiguration() {
    return root.getGraphicsConfiguration();
  }

  Rectangle getBounds() {
    return root.getBounds();
  }

  public SMSettings getSettings() {
    return settings;
  }

  public Frame getRootFrame() {
    return root;
  }

  synchronized public void run() {
    addWidgets();
    root.addWindowListener(this);
    root.pack();
    root.setVisible(true);
    loaded = true;
    notifyAll();
  }

  public void addWidgets() {
    setupMenuBar();
    setupToolBar();

    fc = new JFileChooser(".");
    //Uncomment one of the following lines to try a different
    //file selection mode.  The first allows just directories
    //to be selected (and, at least in the Java look and feel,
    //shown).  The second allows both files and directories
    //to be selected.  If you leave these lines commented out,
    //then the default mode (FILES_ONLY) will be used.
    //
    //fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    //fc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
  }

  public void setupToolBar() {
    JPanel top = new JPanel(new BorderLayout());

    JToolBar toolBar = new JToolBar("toolbar");
    toolBar.setFloatable(true);
    toolBar.setRollover(true);

    addButtons(toolBar);

    desktop = new JDesktopPane();
    docMap = new HashMap<String, Object>();
    settings = new SMSettings();
    ImageIcon icon = createImageIcon("images/stepAna.ico");
    if (icon != null) {
       root.setIconImage(icon.getImage());
    }

    top.setPreferredSize(new Dimension(Frame_Width, Frame_Height));
    top.add(toolBar, BorderLayout.PAGE_START);
    top.add(desktop, BorderLayout.CENTER);
    root.setContentPane(top);
  }

  protected void addButtons(JToolBar toolBar) {
    JButton button = null;

    //first button
    button = makeNavigationButton("New24", F_NEW,
        "New Window", "New Window");
    toolBar.add(button);

    //second button
    button = makeNavigationButton("Open24", F_OPEN,
        "Open a file", "Open");
    toolBar.add(button);

    button = makeNavigationButton("Save24", F_SAVE,
        "Save to a file", "Save");
    toolBar.add(button);

    button = makeNavigationButton("Close24", F_CLOSE,
        "Close this window", "Close");
    toolBar.add(button);

    button = makeNavigationButton("Run24", A_STEPMINER,
        "Run StepMiner", "StepMiner");
    toolBar.add(button);

    button = makeNavigationButton("RunGO24", A_GOANA,
        "Run GO Analysis", "GOAnalysis");
    toolBar.add(button);

  }

  protected JButton makeNavigationButton(String imageName,
      String actionCommand,
      String toolTipText,
      String altText) {
    //Look for the image.
    String imgLocation = "images/" + imageName + ".gif";
    URL imageURL = SMRootFrame.class.getResource(imgLocation);

    //Create and initialize the button.
    JButton button = new JButton();
    button.setActionCommand(actionCommand);
    button.setToolTipText(toolTipText);
    button.addActionListener(this);

    if (imageURL != null) {                      //image found
      button.setIcon(new ImageIcon(imageURL, altText));
    } else {                                     //no image found
      button.setText(altText);
      // System.err.println("Resource not found: " + imgLocation);
    }

    return button;
  }

  public void setupMenuBar() {
    // MenuBar
    JMenuBar menuBar;
    JMenu menu, submenu;
    JMenuItem menuItem;
    JRadioButtonMenuItem rbMenuItem;
    JCheckBoxMenuItem cbMenuItem;

    //Create the menu bar.
    menuBar = new JMenuBar();

    //Build the first menu.
    menu = new JMenu("File");
    menu.setMnemonic(KeyEvent.VK_F);
    menu.getAccessibleContext().setAccessibleDescription(
        "File loading, saving, exporting...");
    menuBar.add(menu);

    //a group of JMenuItems
    menuItem = new JMenuItem(F_OPEN, KeyEvent.VK_O);
    menuItem.setAccelerator(KeyStroke.getKeyStroke(
          KeyEvent.VK_1, ActionEvent.ALT_MASK));
    menuItem.getAccessibleContext().setAccessibleDescription("Open a file");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    menuItem = new JMenuItem(F_NEW, KeyEvent.VK_N);
    menuItem.setAccelerator(KeyStroke.getKeyStroke(
          KeyEvent.VK_2, ActionEvent.ALT_MASK));
    menuItem.getAccessibleContext().setAccessibleDescription("Open a new window");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    menuItem = new JMenuItem(F_SAVE, KeyEvent.VK_S);
    menuItem.setAccelerator(KeyStroke.getKeyStroke(
          KeyEvent.VK_3, ActionEvent.ALT_MASK));
    menuItem.getAccessibleContext().setAccessibleDescription("Save to a file");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    menuItem = new JMenuItem(F_CLOSE, KeyEvent.VK_C);
    menuItem.setAccelerator(KeyStroke.getKeyStroke(
          KeyEvent.VK_4, ActionEvent.ALT_MASK));
    menuItem.getAccessibleContext().setAccessibleDescription("Close this window");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    menuItem = new JMenuItem(F_EXIT, KeyEvent.VK_E);
    menuItem.setAccelerator(KeyStroke.getKeyStroke(
          KeyEvent.VK_5, ActionEvent.ALT_MASK));
    menuItem.getAccessibleContext().setAccessibleDescription("Close all windows and exit");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    //Build the View menu.
    menu = new JMenu("View");
    menu.setMnemonic(KeyEvent.VK_V);
    menu.getAccessibleContext().setAccessibleDescription(
        "View informations...");
    menuBar.add(menu);

    menuItem = new JMenuItem(V_LOG);
    menuItem.getAccessibleContext().setAccessibleDescription("View LOG");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    //Build the Settings menu.
    menu = new JMenu("Settings");
    menu.setMnemonic(KeyEvent.VK_S);
    menu.getAccessibleContext().setAccessibleDescription(
        "View informations...");
    menuBar.add(menu);

    menuItem = new JMenuItem(S_SETTINGS);
    menuItem.getAccessibleContext().setAccessibleDescription("View/Change Settings");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    //Build the Analyze menu.
    menu = new JMenu("Analyze");
    menu.setMnemonic(KeyEvent.VK_A);
    menu.getAccessibleContext().setAccessibleDescription(
        "Analyze timecourse microarray data");
    menuBar.add(menu);

    menuItem = new JMenuItem(A_STEPMINER);
    menuItem.getAccessibleContext().setAccessibleDescription("Run StepMiner");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    menuItem = new JMenuItem(A_GOANA);
    menuItem.getAccessibleContext().setAccessibleDescription("Run GO Analysis");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    //Build the Export menu.
    menu = new JMenu("Export");
    menu.setMnemonic(KeyEvent.VK_E);
    menu.getAccessibleContext().setAccessibleDescription(
        "Analyze timecourse microarray data");
    menuBar.add(menu);

    menuItem = new JMenuItem(E_POSTSCRIPT);
    menuItem.getAccessibleContext().setAccessibleDescription("Export heatmap as Postscript (.ps)");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    menuItem = new JMenuItem(E_IMAGE);
    menuItem.getAccessibleContext().setAccessibleDescription("Export heatmap as Image (.png)");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    menuBar.add(Box.createHorizontalGlue());

    //Build the Help menu.
    menu = new JMenu("Help");
    menu.setMnemonic(KeyEvent.VK_H);
    menu.getAccessibleContext().setAccessibleDescription(
        "Help, Documentations..");
    menuBar.add(menu);

    menuItem = new JMenuItem(H_ABOUT);
    menuItem.getAccessibleContext().setAccessibleDescription("About StepMiner");
    menu.add(menuItem);
    menuItem.addActionListener(this);

    root.setJMenuBar(menuBar);
  }

  /** Returns an ImageIcon, or null if the path was invalid. */
  protected static ImageIcon createImageIcon(String path) {
    java.net.URL imgURL = SMRootFrame.class.getResource(path);
    if (imgURL != null) {
      return new ImageIcon(imgURL);
    } else {
      File f = new File(path);
      if (f.exists()) {
        System.err.println("File exists");
      }
      System.err.println("Couldn't load file: " + path);
      return null;
    }
  }

  // Returns just the class name -- no package info.
  protected String getClassName(Object o) {
    String classString = o.getClass().getName();
    int dotIndex = classString.lastIndexOf(".");
    return classString.substring(dotIndex+1);
  }

  // Close the application only when all the windows (root frames) are closed.
  public void closeThisWindow() {
    int num = getNumRootFrames();
    // System.out.println("Closing " + num);
    if (num <= 1) {
      System.exit(0);
    }
    // Removing this from the doubly linked list
    if (prev != null) {
      prev.next = next;
    }
    if (next != null) {
      next.prev = prev;
    }
    root.dispose();
    // Set everything to null so that garbage collect
    //  can free some memory
    next = null; prev = null; root = null; desktop = null; fc = null;
  }

  public void windowClosing(WindowEvent e) {
    closeThisWindow();
  }

  File openFileDialog() {
    fc.setFileFilter(pclFilter);
    fc.setSelectedFile(new File(""));
    int returnVal = fc.showOpenDialog(desktop);

    if (returnVal == JFileChooser.APPROVE_OPTION) {
      File file = fc.getSelectedFile();
      return file;
    } else {
      return null;
    }
  }

  File saveFileDialog() {
    fc.setSelectedFile(new File(""));
    int returnVal = fc.showSaveDialog(desktop);

    if (returnVal == JFileChooser.APPROVE_OPTION) {
      File file = fc.getSelectedFile();
      return file;
    } else {
      return null;
    }
  }

  String getUniqueName(String filename) {
    int number = 1;
    String name = filename;
    while (docMap.containsKey(name)) {
      name = filename  + "[" + number + "]";
      number ++;
    }
    return name;
  }

  synchronized public void openFile(String filename) {
    try {
      File f = new File(filename);
      if (f.exists()) {
        // System.out.println(f.getAbsoluteFile().getParentFile());
        // System.out.println(fc);
        if (fc != null) {
          fc.setCurrentDirectory(f.getAbsoluteFile().getParentFile());
        }
        String name = getUniqueName(f.getName());
        DocumentFrame tmp = new DocumentFrame(filename, name, this);
        desktop.add(tmp);
        tmp.maximize();
        docMap.put(name, tmp);
      }
      else {
        showDialog("Filename : " + filename + " doesn't exist", "Error");
      }
    }
    catch(Exception e) {
      e.printStackTrace();
      showDialog("Couldn't open " + filename , "Error");
    }
  }

  public boolean hasSelectedDocument() {
     DocumentFrame selected = (DocumentFrame) desktop.getSelectedFrame();
     return (selected != null);
  }

  public void saveFile(String filename) {
     DocumentFrame selected = (DocumentFrame) desktop.getSelectedFrame();
     if (selected != null) {
        selected.saveFile(filename);
     }
     else {
        System.err.println("Select a document before saving");
        showDialog("Select a document before saving", "Error");
     }
  }

  public void saveHeatmap() {
    DocumentFrame selected = (DocumentFrame) desktop.getSelectedFrame();
    if (selected != null) {
      fc.setFileFilter(pngFilter);
      File file = saveFileDialog();
      if (file != null) {
        selected.saveHeatmap(file.getAbsolutePath());
      }
    }
    else {
      System.err.println("Select a document before exporting");
      showDialog("Select a document before exporting", "Error");
    }
  }

  public void runStepMiner() {
     DocumentFrame selected = (DocumentFrame) desktop.getSelectedFrame();
     if (selected != null) {
        fc.setFileFilter(pclFilter);
        File file = saveFileDialog();
        if (file != null) {
          selected.runStepMiner();
          selected.saveFile(file.getAbsolutePath());
        }
     }
     else {
        System.err.println("Select a document before running");
        showDialog("Select a document before running", "Error");
     }
  }

  public void runGOAnalysis() {
     DocumentFrame selected = (DocumentFrame) desktop.getSelectedFrame();
     if (selected != null) {
        fc.setFileFilter(htmlFilter);
        File file = saveFileDialog();
        if (file != null) {
          selected.runGOAnalysis(file.getAbsolutePath());
        }
     }
     else {
        System.err.println("Select a document before running");
        showDialog("Select a document before running", "Error");
     }
  }

  // Show StepMiner Dialogs
  public void showDialog(String message, String title) {
    logger.info(title + ":" + message);
    JOptionPane.showMessageDialog(root, message, title,
        JOptionPane.INFORMATION_MESSAGE);
  }

  // Show StepMiner logs
  public void showLog() {
    try {
      class SMLogDialog extends JDialog implements ActionListener {
        JButton okButton;
        JButton saveButton;
        public SMLogDialog(Frame root, String title, boolean modal) {
          super(root, title, modal);
          Container pane = getContentPane();
          pane.setLayout(new BoxLayout(pane, BoxLayout.Y_AXIS));
          pane.add(log.getOutput());
          JPanel commands = new JPanel(new FlowLayout());
          okButton = new JButton("OK");
          okButton.addActionListener(this);
          commands.add(okButton);
          saveButton = new JButton("Save");
          saveButton.addActionListener(this);
          commands.add(saveButton);
          pane.add(commands);
          pack();
          setVisible(true);
        }
        public void actionPerformed(ActionEvent e) {
          if (e.getSource().equals(okButton)) {
            dispose();
          }
          if (e.getSource().equals(saveButton)) {
            fc.setFileFilter(txtFilter);
            File file = saveFileDialog();
            if (file != null) {
              log.saveFile(file.getAbsolutePath());
            }
          }
        }
      }
      new SMLogDialog(root, "Log", true);
    }
    catch(Exception e) {
      e.printStackTrace();
    }
  }

  public void updateHeatmap() {
    JInternalFrame[] frames = desktop.getAllFrames();
    for (int i =0; i < frames.length; i++) {
      if (frames[i] instanceof DocumentFrame) {
        DocumentFrame selected = (DocumentFrame) frames[i];
        selected.updateHeatmap();
      }
    }
  }

  public void showSettings() {
    settings.showDialog(this, "Settings");
    updateHeatmap();
  }

  public void showAbout() {
    JOptionPane.showMessageDialog(root,
        "StepMiner is created by Debashis Sahoo and David Dill\n" +
        "StepMiner analyzes time course microarray data\n" +
        "Version: " + VERSION);
  }

  public void removeIntFrameName(String intFrameName) {
    docMap.remove(intFrameName);
  }

  public void actionPerformed(ActionEvent e) {
    AbstractButton source = (AbstractButton)(e.getSource());
    String s = "Action event detected."
      + newline
      + "    Event source: " + source.getText()
      + " (an instance of " + getClassName(source) + ")";
    System.out.println(s);
    String event = source.getText();
    if (event.equals(F_EXIT)) {
      System.exit(0);
    }
    if (event.equals(F_CLOSE)) {
      closeThisWindow();
    }
    if (event.equals(F_NEW)) {
      new SMRootFrame(this);
    }
    if (event.equals(F_OPEN)) {
      File file = openFileDialog();
      if (file != null) {
        openFile(file.getAbsolutePath());
      }
    }
    if (event.equals(F_SAVE)) {
      if (hasSelectedDocument()) {
        fc.setFileFilter(pclFilter);
        File file = saveFileDialog();
        if (file != null) {
          saveFile(file.getAbsolutePath());
        }
      }
      else {
        saveFile(null);
      }
    }
    if (event.equals(E_IMAGE)) {
      saveHeatmap();
    }
    if (event.equals(E_POSTSCRIPT)) {
      showDialog("Not implemented", "Error");
    }
    if (event.equals(A_STEPMINER)) {
      runStepMiner();
    } 	
    if (event.equals(A_GOANA)) {
      runGOAnalysis();
    }
    if (event.equals(V_LOG)) {
      showLog();
    }
    if (event.equals(S_SETTINGS)) {
      showSettings();
    }
    if (event.equals(H_ABOUT)) {
      showAbout();
    }
  }

  public void itemStateChanged(ItemEvent e) {
    AbstractButton source = (AbstractButton)(e.getSource());
    String s = "Item event detected." + newline
      + "    Event source: " + source.getText()
      + " (an instance of " + getClassName(source) + ")" + newline
      + "    New state: "
      + ((e.getStateChange() == ItemEvent.SELECTED) ? "selected":"unselected");
    System.out.println(s);
  }
}
