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

package tools.graphs;

import java.io.*;
import java.awt.Color;
import java.util.*;
import java.text.MessageFormat;

public class PSPlot extends Plot {
  int rows_;
  int cols_;
  int currentIndex_;
  int maxIndex_;
  int currentPage_;
  int maxPage_;

  double paperWidth_;
  double paperHeight_;
  double topMargin_;
  double leftMargin_;
  double textWidth_;
  double textHeight_;
  double plotMargin_;

  double basex1_;
  double basey1_;
  double basex2_;
  double basey2_;
  double originx_;
  double originy_;
  double xlen_;
  double ylen_;

  double[] xrange_; // low, hi, tick
  double[] yrange_; // low, hi, tick

  public PSPlot(String filename) throws IOException {
    super(filename);
    rows_ = 1;
    cols_ = 1;
    currentIndex_ = -1;
    maxIndex_ = 0;
    currentPage_ = 1;
    maxPage_ = 1;
    xrange_ = new double[3];
    yrange_ = new double[3];
    xrange_[0] = 0; xrange_[1] = 1; xrange_[2] = 0.5;
    yrange_[0] = 0; yrange_[1] = 1; yrange_[2] = 0.5;
    setupLetterPaper();
  }

  public void open() throws IOException {
    String header = ""
      + "%!PS-Adobe-3.0\n"
      + "%%Creator: Java StepMiner - by Debashis Sahoo (Stanford University)\n"
      + "%%Title: StepMiner Plots \n"
      + "%%Pages: (atend)\n"
      + "%%PageOrder: Ascend\n"
      + "%%BoundingBox: 0 0 596 842\n"
      + "%%EndComments\n"
      + "%%BeginProlog\n"
      + "\n"
      + "/gs  { gsave } def\n"
      + "/gr  { grestore } def\n"
      + "/bp  { gs gs } def\n"
      + "/ep  { showpage gr gr } def\n"
      + "/in  {72 mul} def      % Convert inches->points (1/72 inch)\n"
      + "/rgb { setrgbcolor } def\n"
      + "/s   { scalefont setfont } def\n"
      + "/m   { moveto } def\n"
      + "/l   { rlineto } def\n"
      + "/o   { stroke } def\n"
      + "/c   { newpath 0 360 arc } def\n"
      + "/r   { 4 2 roll moveto 1 copy 3 -1 roll exch 0 exch rlineto 0 rlineto -1 mul 0 exch rlineto closepath } def\n"
      + "/t   { 6 -2 roll moveto gsave rotate\n"
      + "       ps mul neg 0 2 1 roll rmoveto\n"
      + "       1 index stringwidth pop\n"
      + "       mul neg 0 rmoveto show grestore } def\n"
      + "/cl  { grestore gsave newpath 3 index 3 index moveto 1 index\n"
      + "       4 -1 roll lineto  exch 1 index lineto lineto\n"
      + "       closepath clip newpath } def\n"
      + "\n"
      + "%%EndProlog\n"
      + "\n";
    print(header);
    print("%%Page: " + currentPage_ + " " + currentPage_ + "\nbp\n");
    setFont();
  }

  public void array(int row, int col) {
    rows_ = row;
    cols_ = col;
    maxIndex_ = row * col - 1;
    currentIndex_ = -1;
  }

  public void setupLetterPaper() {
    paperWidth_ = 8.5;
    paperHeight_ = 11.0;
    topMargin_ = 0.5;
    leftMargin_ = 0.5;
    textWidth_ = paperWidth_ - 2 * leftMargin_;
    textHeight_ = paperHeight_ - 2 * topMargin_;
    plotMargin_ = 0.25;
  }

  public void computeBase() throws IOException {
    double plotWidth = (textWidth_ - (cols_+1) * plotMargin_) / cols_;
    double plotHeight = (textHeight_ - (rows_+1) * plotMargin_) / rows_;
    int col = (currentIndex_ % cols_);
    int row = (currentIndex_ / cols_);
    basex1_ = leftMargin_ + plotMargin_ + col *
      (plotMargin_+plotWidth);
    basex2_ = basex1_ + plotWidth;
    basey1_ = paperHeight_ - topMargin_ - (row + 1) *
      (plotMargin_ + plotHeight);
    basey2_ = basey1_ + plotHeight;
    // Borders
    String res =  ""
      + basex1_ + " in " + basey1_ + " in m\n"
      + (basex2_-basex1_) + " in " + 0 + " in l\n"
      + 0 + " in " + (basey2_-basey1_) + " in l\n"
      + (basex1_-basex2_) + " in " + 0 + " in l\n"
      + 0 + " in " + (basey1_-basey2_) + " in l o\n";
    // print(res);
    originx_ = basex1_ + 0.35;
    originy_ = basey1_ + 0.35;
    xlen_ = plotWidth - 0.7;
    ylen_ = plotHeight - 0.7;
    // Axes :
    res =  ""
      + originx_ + " in " + (originy_ + ylen_) + " in m\n"
      + 0 + " in " + (-ylen_) + " in l\n"
      + xlen_ + " in " + 0 + " in l o\n";
    print(res);
  }

  public int computeOrd(double val) {
    if (val < 0) {
      val = - val;
    }
    if (val == 0) {
      return 0;
    }
    return (int)Math.floor(Math.log(val)/Math.log(10));
  }

  public int computeFactor(double val, boolean min) {
    int ord = computeOrd(val);
    int factor = 0;
    if (min) {
      factor = (int)Math.floor(val/Math.pow(10, ord));
    }
    else {
      factor = (int)Math.ceil(val/Math.pow(10, ord));
    }
    return factor;
  }

  public double computeRound(double val, boolean min) {
    int ord = computeOrd(val);
    int factor = 0;
    if (min) {
      factor = (int)Math.floor(val/Math.pow(10, ord));
    }
    else {
      factor = (int)Math.ceil(val/Math.pow(10, ord));
    }
    return factor * Math.pow(10, ord);
  }

  public double computeTick(double[] range) {
    int factor = computeFactor(range[1]-range[0], true);
    int ord = computeOrd(range[1]-range[0]);
    if (factor > 0 && factor < 5) {
      factor = 5;
    }
    else {
      factor = 1;
    }
    return Math.pow(10, ord)/factor;
  }

  public void computeRange(Vector<Double> v, double[] range) 
    throws IOException {
      Enumeration<Double> e = v.elements();
      double min = Double.MAX_VALUE;
      double max = Double.MIN_VALUE;
      int count = 0;
      while (e.hasMoreElements()) {
        Double val = (Double) e.nextElement();
        if (val != null && min > val.doubleValue()) {
          min = val.doubleValue();
        }
        if (val != null && max < val.doubleValue()) {
          max = val.doubleValue();
        }
        if (val != null) {
          count++;
        }
      }
      //print("% " + min + " " + convertString(max)  + "\n");
      if (min > max) {
        min = -1; max = 1;
      }
      if (min == max) {
        min -= 1; max += 1;
      }
      //print("% " + min + " " + convertString(max)  + "\n");
      range[0] = min;
      range[1] = max;
      double tick = computeTick(range);
      range[0] = computeRound(min, true);
      int numTicks = (int) Math.ceil((max - range[0]) / tick);
      range[1] = range[0] + numTicks * tick;
      range[0] = range[0] - tick;
      range[2] = tick;
    }

  public void plotXticks() throws IOException {
    int numTicks = (int) Math.round((xrange_[1] - xrange_[0]) / xrange_[2]);
    print("/ps 8 def\n");
    for (int i = 1; i <= numTicks; i++) {
      double t = xrange_[0] + i * xrange_[2];
      print( getXScale(t) + convertString(originy_) + " m\n");
      print( "0 -4 l o\n");
      print(getXScale(t) + convertString(originy_-0.15) + 
          " (" + getString(t) + ") .5 0 0 t\n");
    }
  }

  public void plotYticks() throws IOException {
    int numTicks = (int) Math.round((yrange_[1] - yrange_[0]) / yrange_[2]);
    print("/ps 8 def\n");
    int skip1 = numTicks / 6;
    if (skip1 == 0) {
        skip1 = 1;
    }
    int skip = skip1 * 3;
    for (int i = 0; i <= numTicks; i++) {
      double t = yrange_[0] + i * yrange_[2];
      if ((i % skip1) == 0) {
        print( convertString(originx_) + getYScale(t) + " m\n");
        print( "-4 0 l o\n");
      }
      if ( (i % skip) == 0) {
        print( convertString(originx_-0.1) + getYScale(t) + 
            " (" + getString(t) + ") .5 0 90 t\n");
      }
    }
  }

  public void plot(Double[] x, Double[] y) throws IOException {
    Vector<Double> vx = new Vector<Double>(Arrays.asList(x));
    Vector<Double> vy = new Vector<Double>(Arrays.asList(y));
    plot(vx, vy);
  }

  public void setupPlot(Double[] x, Double[] y) throws IOException {
    Vector<Double> vx = new Vector<Double>(Arrays.asList(x));
    Vector<Double> vy = new Vector<Double>(Arrays.asList(y));
    setupPlot(vx, vy);
  }

  public void setFont() throws IOException {
    // print("/Helvetica findfont 8 s\n");
    print("/Times-Roman findfont 8 s\n");
  }

  public void setFont(String font, int size) throws IOException {
    print("/"+font+" findfont "+ size +" s\n");
  }

  public void newPage() throws IOException {
    if (currentPage_ > 0) {
      print("ep\n");
    }
    currentPage_ ++;
    print("%%Page: " + currentPage_ + " " + currentPage_ + "\nbp\n");
    setFont();
    print("0 0 0 rgb\n");
    print("0.75 setlinewidth\n");
  }

  public void setupPlot(Vector<Double> x, Vector<Double> y) throws IOException {
    currentIndex_ ++;
    if (currentIndex_ > maxIndex_) {
      newPage();
      currentIndex_ = 0;
    }
    //print("% " + currentIndex_ + " " + currentPage_ + "\n");
    computeBase();
    computeRange(x, xrange_);
    computeRange(y, yrange_);
    //print("% X" + xrange_[0] + " " + xrange_[1] + " " + xrange_[2] + "\n");
    //print("% Y" + yrange_[0] + " " + yrange_[1] + " " + yrange_[2] + "\n");
    //print("% O" + originx_ + "-" + originy_ + " " + xlen_ + ":" + ylen_ + "\n");
    plotXticks();
    plotYticks();
  }

  public void plotPoints(Vector<Double> x, Vector<Double> y) throws IOException {
    Object[] xx = x.toArray();
    Object[] yy = y.toArray();
    for (int i = 0;i < xx.length; i++) {
        point((Double)xx[i], (Double)yy[i]);
    }
  }

  public void plot(Vector<Double> x, Vector<Double> y) throws IOException {
    setupPlot(x,y);
    plotPoints(x,y);
  }

  // Relative scale
  public double convertRefXScale(double x) {
    // range[0] -> originx_
    // range[1] -> originx_ + xlen_
    //  x -> xlen_ / (range[1] - range[0]) * x
    return  xlen_ / (xrange_[1] - xrange_[0]) * x;
  }

  public double convertRefYScale(double y) {
    return  ylen_ / (yrange_[1] - yrange_[0]) * y;
  }

  // Absolute scale
  public double convertXScale(double x) {
    return  originx_ + convertRefXScale(x - xrange_[0]);
  }

  public double convertYScale(double y) {
    return  originy_ + convertRefYScale(y - yrange_[0]);
  }

  public static String formatString(String format, double x) {
    MessageFormat mf = new MessageFormat(
        "{0,number," + format + "}");
    Object[] objs = {new Double(x)};
    return mf.format(objs);
  }

  public static String getString(double x) {
    String res = formatString("0.##", x);
    return res;
  }

  public static String convertString(double x) {
    String res = formatString("0.##", x);
    return res + " in ";
  }

  public String getXScale(double x) {
    return convertString(convertXScale(x));
  }

  public String getYScale(double y) {
    return convertString(convertYScale(y));
  }

  public String getRefXScale(double x) {
    return convertString(convertRefXScale(x));
  }

  public String getRefYScale(double y) {
    return convertString(convertRefYScale(y));
  }

  public void verticalLine(Double x1) throws IOException {
    line(x1, new Double(yrange_[0]), x1, new Double(yrange_[1]));
  }

  public void horizontalLine(Double y1) throws IOException {
    line(new Double(xrange_[0]), y1, new Double(xrange_[1]), y1);
  }

  public void line(Double x1, Double y1, Double x2, Double y2) 
    throws IOException {
      if (x1 == null || y1 == null || x2 == null || y2 == null) {
        return;
      }
      //print("% Line: " + x1 + " " + y1 + " " + x2 + " " + y2 + "\n");
      print(getXScale(x1.doubleValue())
          + getYScale(y1.doubleValue()) + " m\n" 
          + convertString(convertXScale(x2.doubleValue()) - convertXScale(x1.doubleValue()))
          + convertString(convertYScale(y2.doubleValue()) - convertYScale(y1.doubleValue()) )+ " l o\n");
    }

  public void point(Double x1, Double y1) throws IOException {
      if (x1 == null || y1 == null) {
        return;
      }
      //print("% Point: " + x1 + " " + y1 + "\n");
      print(getXScale(x1.doubleValue())
          + getYScale(y1.doubleValue()) + " 1.78 c o\n");
    }

  public void text(Double x, Double y, int rotate, String s) throws IOException {
    print("/ps 10 def " + getXScale(x.doubleValue())
        + getYScale(y.doubleValue()) + " (" + s + ") .5 0 "+rotate+" t\n");
  }
  public void title(String s) throws IOException {
    print("/ps 10 def " + convertString((basex1_ + basex2_)/2)
        + convertString(basey2_ - 0.1) + " (" + s + ") .5 0 0 t\n");
  }
  public void xlabel(String s) throws IOException {
    print("/ps 10 def " + convertString((basex1_ + basex2_)/2)
        + convertString(basey1_ + 0.1) + " (" + s + ") .5 0 0 t\n");
  }
  public void ylabel(String s) throws IOException {
    print("/ps 10 def " + convertString(basex1_ + 0.1)
        + convertString((basey1_+basey2_)/2) + " (" + s + ") .5 0 90 t\n");
  }

  public void setRGBcolor(double r, double g, double b) throws IOException {
    print(r + " " + g + " " + b + " rgb\n");
  }

  public void setRGBcolor(Color c) throws IOException {
    setRGBcolor(c.getRed()/255.0, c.getGreen()/255.0, c.getBlue()/255.0);
  }

  public void close() throws IOException {
    String end = "ep\n"
      + "%%Trailer\n"
      + "%%Pages: " + currentPage_ + "\n"
      + "%%EOF\n";
    print(end);
    super.close();
  }
}

