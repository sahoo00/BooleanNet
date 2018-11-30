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

package tools;

import tools.microarray.Data;
import tools.microarray.FileReader.PCLFileReader;
import tools.microarray.FileWriter.PCLFileWriter;
import tools.microarray.StepMiner.StepMiner;
import tools.microarray.ArrayOrder;
import tools.microarray.Impute;
import java.io.*;
import java.util.*;

import jargs.gnu.CmdLineParser;
import tools.microarray.GeneNameScheme;

public class Analyze {

  private static void printUsage() {                                          
    String out = 
      "Usage: tools.Analyze [-d/--debug] <filename>\n" +
      "                     [-r/--reduceLog] <filename>\n" +
      "                     --gui [<filename>]\n" +
      "                     [options] <filename>\n" +
      "Options:                                 \n" +
      "       [-o/--outfile <file.ano>] [-o/--outfile <file.pcl>]\n" +
      "       [-o/--outfile [<tag>:]<file.ann>] [-o/--outfile <file.exp>]\n" +
      "       [-o/--outfile <file.cdt>] [-o/--outfile <file.gmx>]\n" +
      "       [-o/--outfile [<tag>:]<file.ps>] [-o/--outfile <file.png>]\n" +
      "       Example tags: All, Step, Up, Down, UpDown, DownUp, Rest\n" +
      "       [-t/--type <type>]\n" +
      "             <type> : OneStep, OneStepFdr, Fdr, \n" +
      "                      Order, Subset, ZeroCenter, MeanCenter, \n" +
      "                      None, ListGenes, Normalize,\n" +
      "                      KNNimpute, LLSimpute\n" +
      "       [--annFile <Gene Annotation File : 15 columns format>]\n" +
      "       [--onnFile <Ontology File : OBO format>]\n" +
      "       [--org <Organism: Mm/Hs/Sgd/Pombie>]\n" +
      "       [--geneIndex <arrayIndex of gene description>]\n" +
      "       [--splitString <Splitting regexp of the gene str>]\n" +
      "       [--splitIndex <Index of gene after splitting>]\n" +
      "       [--goPvalue <pvalue threshold of GOAnalysis>]\n" +
      "       [--range <ex: 4:17 Range of array indices for analysis>]\n" +
      "       [-p/--pvalue <pvalue threshold>]\n" +
      "       [--numMissing <Number of missing timepoints>]\n" +
      "       [--Intersect <file>]\n" +
      "       [--Union <file>]\n" +
      "       [--Select <file>] select ids with original order\n" +
      "       [--SelectOrder <file>] select ids with given order\n" +
      "       [--Diff <file>]\n" +
      "       [--SelectNames <file>] select names with original order\n" +
      "                                         \n";
    System.err.print(out);
  }  

  public static void printCurrentDirectory() {
    File f = new File(".");
    String[] list = f.list();
    for (int i=0; i<list.length; i++) {
      System.out.println(list[i]);
    }
  }

  public static void test() {
    GeneNameScheme ns = new GeneNameScheme();
    String a =  " || proteosome (prosome, macropain) subunit, beta type 9 (large || 2210417M10 || Proteosome (prosome, macropain) subunit, beta type 9 (large multifunctional peptidase 2) || Psmb9 || AV081184 || Mm.390983 ||  || 197923";
    String gene = ns.getGene(a);
    System.out.println(":"+gene+":");
  }

  public static void runStepMiner(String file, String type, Double pvalue,
    Vector outFiles, GeneNameScheme ns, Double goPvalue, String range) {
    try {
      Data data = PCLFileReader.readFile(file);
      data.setGeneNameScheme(ns);
      data.setRange(range);
      data.convertDoubles();
      if (type.startsWith("Order") || type.equals("Subset") ||
          type.equals("ZeroCenter") || type.equals("MeanCenter") ||
          type.equals("None") || type.equals("ListGenes") ||
          type.equals("Normalize") ) {
        // Ordering of Arrays
        ArrayOrder arr = new ArrayOrder(data);
        if (type.equals("Order")) {
          arr.calculateOrderOld();
          arr.performRandomExpt();
        }
        if (type.equals("OrderTree")) {
          arr.calculateOrder();
        }
        if (type.equals("OrderKL")) {
          arr.calculateOrderKL();
        }
        if (type.equals("MeanCenter")) {
          arr.meanCenter();
        }
        if (type.equals("ZeroCenter")) {
          arr.zeroCenter();
        }
        if (type.equals("Normalize")) {
          arr.normalize();
        }
        Data dataOut = data;
        if (type.startsWith("Order")) {
          dataOut = arr.getOrderedDataSorted();
        }
        if (type.equals("Subset") || type.equals("MeanCenter") ||
            type.equals("ZeroCenter") ||  type.equals("Normalize") ) {
          dataOut = arr.getOrderedData();
        }
        if (type.equals("ListGenes")) {
          dataOut = arr.listGenes();
        }
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String f = (String) e.nextElement();
          if (f.endsWith(".pcl")) {
            PCLFileWriter.writeFile(dataOut, f, null);
          }
          if (f.endsWith(".html")) {
            arr.performGOAnalysis(dataOut, f, goPvalue);
          }
          if (f.endsWith(".ps") || f.endsWith(".eps")) {
            arr.plotGenes(dataOut, f);
          }
        }
      }
      else if (type.equals("KNNimpute") || type.equals("LLSimpute")) {
        Impute im = new Impute(data, type);
        im.impute(10);
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String f = (String) e.nextElement();
          if (f.endsWith(".pcl")) {
            PCLFileWriter.writeFile(data, f, null);
          }
        }
      }
      else { // StepMiner Analysis
        StepMiner sm = new StepMiner(data);
        if (type.startsWith("OneStep")) {
          sm.setOneStepAnalysis();
        }
        if (type.startsWith("OneStepFdr")) {
          sm.setFdrAnalysis(true);
        }
        if (type.startsWith("Fdr")) {
          sm.setFdrAnalysis(true);
        }
        if (type.contains("NoCentering")) {
          sm.setStepCentering(false);
        }
        else {
          sm.setStepCentering(true);
        }
        sm.setDeleteGenesInMultipleGroups(false);
        sm.setPvalueThr(pvalue.doubleValue());
        sm.performAnalysis();
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String f = (String) e.nextElement();
          if (f.endsWith(".pcl")) {
            sm.writePCL(f);
          }
          if (f.endsWith(".ano")) {
            PCLFileWriter.writeFile(sm.getStepAnnotationData(), f);
          }
          if (f.endsWith(".ing")) {
            PCLFileWriter.writeFile(sm.getStepAnnotationByGenes(), f);
          }
          if (f.endsWith(".cdt")) {
            PCLFileWriter.writeFile(sm.getCdtAnnotationData(), f);
            String gtr = f.replaceFirst(".cdt$", ".gtr");
            PCLFileWriter.writeFile(sm.getGtrAnnotationData(), gtr);
          }
          if (f.endsWith(".html")) {
            sm.performGOAnalysis(f, goPvalue);
          }
          if (f.endsWith(".htm")) {
            sm.performGOBestAnalysis(f, goPvalue);
          }
          if (f.endsWith(".ps") || f.endsWith(".eps")) {
            sm.plotSteps(f);
          }
          if (f.endsWith(".ann")) {
            sm.writeAnnotations(f);
          }
          if (f.endsWith(".gmt") || f.endsWith(".gxa") || f.endsWith(".tab")) {
            sm.writeGeneSets(f);
          }
        }
      }
    }
    catch (Exception e) {
      e.printStackTrace();
    }
  }

  public static void main(String[] args) throws Exception {

    // Get current time
    long start = System.currentTimeMillis();
    
    CmdLineParser parser = new CmdLineParser();
    CmdLineParser.Option debugO = parser.addBooleanOption('d', "debug");
    CmdLineParser.Option logreduce0 = parser.addBooleanOption('r', "reduceLog");
    CmdLineParser.Option out0 = parser.addStringOption('o', "outfile");
    CmdLineParser.Option type0 = parser.addStringOption('t', "type");
    CmdLineParser.Option pvalue0 = parser.addDoubleOption('p', "pvalue");
    CmdLineParser.Option annFile0 = parser.addStringOption("annFile");
    CmdLineParser.Option onnFile0 = parser.addStringOption("onnFile");
    CmdLineParser.Option org0 = parser.addStringOption("org");
    CmdLineParser.Option geneIndex0 = parser.addIntegerOption("geneIndex");
    CmdLineParser.Option splitString0 = parser.addStringOption("splitString");
    CmdLineParser.Option splitIndex0 = parser.addIntegerOption("splitIndex");
    CmdLineParser.Option goPvalue0 = parser.addDoubleOption("goPvalue");
    CmdLineParser.Option range0 = parser.addStringOption("range");
    CmdLineParser.Option numMissing0 = parser.addIntegerOption("numMissing");
    CmdLineParser.Option intersect0 = parser.addStringOption("Intersect");
    CmdLineParser.Option union0 = parser.addStringOption("Union");
    CmdLineParser.Option select0 = parser.addStringOption("Select");
    CmdLineParser.Option selectOrder0 = parser.addStringOption("SelectOrder");
    CmdLineParser.Option diff0 = parser.addStringOption("Diff");
    CmdLineParser.Option selectNames0 = parser.addStringOption("SelectNames");
    CmdLineParser.Option gui0 = parser.addBooleanOption("gui");

    try {
      parser.parse(args);
    }
    catch (CmdLineParser.OptionException e) {
      System.err.println(e.getMessage());
      printUsage();
      System.exit(2);
    }

    Boolean debug = (Boolean)parser.getOptionValue(debugO, Boolean.FALSE);
    Boolean logreduce = (Boolean)parser.getOptionValue(logreduce0, Boolean.FALSE);
    String type = (String)parser.getOptionValue(type0, "BothStep");
    Double pvalue = (Double)parser.getOptionValue(pvalue0, new Double(0.05));
    String annFile = (String)parser.getOptionValue(annFile0, null);
    String onnFile = (String)parser.getOptionValue(onnFile0, null);
    String org = (String)parser.getOptionValue(org0, "Sgd");
    Integer geneIndex=(Integer)parser.getOptionValue(geneIndex0,new Integer(0));
    String splitString=(String)parser.getOptionValue(splitString0,"\\|\\|");
    Integer splitIndex=(Integer)parser.getOptionValue(splitIndex0,new Integer(-1));
    Double goPvalue = (Double)parser.getOptionValue(goPvalue0, new Double(0.05));
    Integer numMissing =(Integer)parser.getOptionValue(numMissing0,new Integer(0));
    String range=(String)parser.getOptionValue(range0,null);
    String intersectFile =(String)parser.getOptionValue(intersect0,null);
    String unionFile =(String)parser.getOptionValue(union0,null);
    String selectFile =(String)parser.getOptionValue(select0,null);
    String selectOrderFile =(String)parser.getOptionValue(selectOrder0,null);
    String diffFile =(String)parser.getOptionValue(diff0,null);
    String selectNamesFile =(String)parser.getOptionValue(selectNames0,null);
    Boolean gui =(Boolean)parser.getOptionValue(gui0, Boolean.FALSE);

    Vector outFiles = parser.getOptionValues(out0);

    String [] otherArgs = parser.getRemainingArgs();
    if (gui.booleanValue()) {
      SMGui.main(otherArgs);
      return;
    }
    if (otherArgs.length != 1) {
      printUsage();
      System.exit(2);
    }

    try {
      if (debug.booleanValue()) {
        System.out.println("Debugging is enabled.");
        test();
        System.exit(0);
      }

      if (logreduce.booleanValue()) {
        System.out.println("Log Reducing data ...");
        Data data = PCLFileReader.readFile(otherArgs[0]);
        data.reduceLog();
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String file = (String) e.nextElement();
          PCLFileWriter.writeFile(data, file);
          break;
        }
        System.exit(0);
      }

      if (intersectFile != null) {
        System.out.println("Intersecting data ...");
        Data data1 = PCLFileReader.readFile(otherArgs[0]);
        Data data2 = PCLFileReader.readFile(intersectFile);
        data1.setRange(range);
        data2.setRange(range);
        Data data = Data.mergeData(data1, data2, true, false);
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String file = (String) e.nextElement();
          PCLFileWriter.writeFile(data, file);
          break;
        }
        System.exit(0);
      }

      if (selectFile != null) {
        System.out.println("Selecting data ...");
        Data data1 = PCLFileReader.readFile(otherArgs[0]);
        Data data2 = PCLFileReader.readFile(selectFile);
        data1.setRange(range);
        data2.setRange(range);
        Data data = Data.mergeData(data1, data2, true, true);
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String file = (String) e.nextElement();
          PCLFileWriter.writeFile(data, file);
          break;
        }
        System.exit(0);
      }

      if (selectOrderFile != null) {
        System.out.println("Selecting data ...");
        Data data1 = PCLFileReader.readFile(otherArgs[0]);
        Data data2 = PCLFileReader.readFile(selectOrderFile);
        data1.setRange(range);
        data2.setRange(range);
        Data data = Data.selectOrder(data1, data2);
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String file = (String) e.nextElement();
          PCLFileWriter.writeFile(data, file);
          break;
        }
        System.exit(0);
      }

      if (diffFile != null) {
        System.out.println("Diff ...");
        Data data1 = PCLFileReader.readFile(otherArgs[0]);
        Data data2 = PCLFileReader.readFile(diffFile);
        data1.setRange(range);
        data2.setRange(range);
        Data data = Data.diffData(data1, data2);
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String file = (String) e.nextElement();
          PCLFileWriter.writeFile(data, file);
          break;
        }
        System.exit(0);
      }

      if (unionFile != null) {
        System.out.println("Union ...");
        Data data1 = PCLFileReader.readFile(otherArgs[0]);
        Data data2 = PCLFileReader.readFile(unionFile);
        data1.setRange(range);
        data2.setRange(range);
        Data data = Data.mergeData(data1, data2, false, false);
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String file = (String) e.nextElement();
          PCLFileWriter.writeFile(data, file);
          break;
        }
        System.exit(0);
      }

      GeneNameScheme ns = new GeneNameScheme(geneIndex.intValue(), 
          org, splitString, splitIndex.intValue());
      ns.setAnnotationFile(annFile);
      ns.setOntologyFile(onnFile);
      ns.setNumMissingPoints(numMissing.intValue());

      if (selectNamesFile != null) {
        System.out.println("Selecting Names ...");
        Data data1 = PCLFileReader.readFile(otherArgs[0]);
        Data data2 = PCLFileReader.readFile(selectNamesFile);
        data1.setRange(range);
        data2.setRange(range);
        data1.setGeneNameScheme(ns);
        data2.setGeneNameScheme(ns);
        Data data = Data.selectNames(data1, data2);
        Enumeration e = outFiles.elements();
        while (e.hasMoreElements()) {
          String file = (String) e.nextElement();
          PCLFileWriter.writeFile(data, file);
          break;
        }
        System.exit(0);
      }

      runStepMiner(otherArgs[0], type, pvalue, outFiles, ns, goPvalue, range);
    }
    catch(Exception e) {
        e.printStackTrace();
    }
    // Get elapsed time in milliseconds
    long elapsedTimeMillis = System.currentTimeMillis()-start;
    
    // Get elapsed time in seconds
    float elapsedTimeSec = elapsedTimeMillis/1000F;
    System.out.println("Elapsed Time : " + elapsedTimeSec + " sec");
  }

};

