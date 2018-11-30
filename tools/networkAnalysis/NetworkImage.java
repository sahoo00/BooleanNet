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

package tools.networkAnalysis;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import java.awt.geom.*;
import java.awt.font.*;
import javax.imageio.*;

import java.util.*;
import java.io.*;

import tools.graphs.*;
import tools.graphs.parser.GML;

public class NetworkImage {

  Graph graph_;
  double left_, right_, bottom_, top_;
  // left_ < right_
  // top_ < bottom_
  double originx_, originy_;
  int imageWidth_;
  int imageHeight_;

  GraphicsConfiguration gc_;

  BufferedImage image_;

  public NetworkImage(Graph g) {
    graph_ = g;
    GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
    GraphicsDevice gs = ge.getDefaultScreenDevice();
    gc_ = gs.getDefaultConfiguration();
    imageWidth_ = 640;
    imageHeight_ = 480;
    init();
  }

  // 0 < ratio < 1 - zoom in
  // -1 < ratio < 0 - zoom out
  public void zoom(double ratio) {
    double xskew = ratio * (right_-left_);
    double yskew = ratio * (bottom_-top_);
    left_ -= xskew; right_ += xskew;
    top_ -= yskew; bottom_ += yskew;
  }

  public void adjustAspect() {
    double aspect = (right_-left_)/(bottom_-top_);
    if ((aspect * imageHeight_) < imageWidth_) {
        double shift = (imageWidth_ * (bottom_-top_)/imageHeight_) - (right_-left_);
        shift = shift / 2;
        left_ = left_ - shift;
        right_ = right_ + shift;
    }
    else {
        double shift = (imageHeight_ * (right_-left_)/imageWidth_) - (bottom_-top_);
        shift = shift / 2;
        bottom_ = bottom_ + shift;
        top_ = top_ - shift;
    }
  }

  public void findBoundaries() {
    left_ = 0; right_ = 0; top_ = 0; bottom_ = 0;
    // Traverse the graph
    class GraphDFSVisitor implements DFSVisitor {
        public boolean visitBefore(Graph g, Node v, Node parent) {
            return true;
        }
        public boolean visit(Graph g, Node v, Node parent) {
            String id = v.getID();
            HashMap map = (HashMap) v.getAttribute("graphics");
            if (map != null) {
              Double x = (Double) map.get("x");
              Double y = (Double) map.get("y");
              Double w = (Double) map.get("w");
              Double h = (Double) map.get("h");
              if (left_ > (x.doubleValue() - w.doubleValue())) {
                left_ = (x.doubleValue() - w.doubleValue());
              }
              if (right_ < (x.doubleValue() + w.doubleValue())) {
                right_ = (x.doubleValue() + w.doubleValue());
              }
              if (top_ > (y.doubleValue() - h.doubleValue())) {
                top_ = (y.doubleValue() - h.doubleValue());
              }
              if (bottom_ < (y.doubleValue() - h.doubleValue())) {
                bottom_ = (y.doubleValue() + h.doubleValue());
              }
            }
            return true;
        }
        public boolean visitAfter(Graph g, Node v, Node parent) {
            return true;
        }
    };
    graph_.traverseDFS(new GraphDFSVisitor());
    System.out.println("Size : " + left_ + "x" + top_+ "+" + 
        (int) (right_-left_) + "+" + (int) (bottom_-top_) );
  }

  public static String trimDoubleQuotes(String text) {
    int textLength = text.length();

    if (textLength >= 2 && text.charAt(0) == '"' && text.charAt(textLength - 1) == '"') {
      return text.substring(1, textLength - 1);
    }

    return text;
  }

  public static Color parseColor(String nm ) throws NumberFormatException {
    nm = trimDoubleQuotes(nm);
    //System.out.println("nm=" + nm );
    if ( nm.startsWith("#") ) {
      nm = nm.substring(1);
    }
    nm = nm.toLowerCase();
    if (nm.length() > 6) {
      throw new NumberFormatException("nm is not a 24 bit representation of the color, string too long"); 
    }
    //System.out.println("nm=" + nm );
    Color color = new Color( Integer.parseInt( nm , 16 ) );
    return color;
  }

  public Double getX(Double v) {
    double res = (v.doubleValue()-left_) * imageWidth_/(right_-left_);
    return new Double(res);
  }

  public Double getY(Double v) {
    double res = (v.doubleValue()-top_) * imageHeight_/(bottom_-top_);
    return new Double(res);
  }

  public Double getW(Double v) {
    double res = v.doubleValue() * imageWidth_/(right_-left_);
    return new Double(res);
  }

  public Double getH(Double v) {
    double res = v.doubleValue() * imageHeight_/(bottom_-top_);
    return new Double(res);
  }

  public void paintImage(Graphics2D g2) {
    // Traverse the graph
    class GraphDFSVisitor implements DFSVisitor {
        Graphics2D g2_;
        HashSet<Edge> hash_;
        public GraphDFSVisitor(Graphics2D g2) {
            g2_ = g2;
            hash_ = new HashSet<Edge>();
        }
        public boolean visitBefore(Graph g, Node v, Node parent) {
            Enumeration<Edge> e = v.getOutGoingEdges();
            while (e.hasMoreElements()) {
                Edge edge = (Edge) e.nextElement();
                if (!hash_.contains(edge)) {
                    hash_.add(edge);
                    visitEdge(g, v, edge);
                }
            }
            return true;
        }
        public void visitEdge(Graph g, Node v, Edge edge) {
          Node dst = edge.getDestination(v);
          HashMap map = (HashMap) v.getAttribute("graphics");
          Double srcx = getX((Double) map.get("x"));
          Double srcy = getY((Double) map.get("y"));
          Double srcw = getW((Double) map.get("w"));
          Double srch = getH((Double) map.get("h"));
          map = (HashMap) dst.getAttribute("graphics");
          Double dstx = getX((Double) map.get("x"));
          Double dsty = getY((Double) map.get("y"));
          Double dstw = getW((Double) map.get("w"));
          Double dsth = getH((Double) map.get("h"));
          Line2D.Double line = new Line2D.Double(srcx.doubleValue(), srcy.doubleValue(), dstx.doubleValue(), dsty.doubleValue());
          g2_.setPaint(Color.black);
          g2_.draw(line);
          //TODO: Draw arrow for the directed lines
        }
        public Ellipse2D.Double getEllipse(Double x, Double y, Double w, Double h) {
          double x1 = x.doubleValue() - w.doubleValue();
          double y1 = y.doubleValue() - h.doubleValue();
          double w1 = w.doubleValue() * 2;
          double h1 = h.doubleValue() * 2;
          Ellipse2D.Double ellipse = new Ellipse2D.Double(x1, y1, w1, h1);
            return ellipse;
        }
        public boolean visit(Graph g, Node v, Node parent) {
            String id = v.getID();
            String label = trimDoubleQuotes((String)v.getAttribute("label"));
            HashMap map = (HashMap) v.getAttribute("graphics");
            Font f = new Font("Courier New", Font.PLAIN,  12);
            FontRenderContext frc = g2_.getFontRenderContext();
            g2_.setFont(f);
            g2_.setRenderingHint(RenderingHints.KEY_ANTIALIASING, // Anti-alias!
                RenderingHints.VALUE_ANTIALIAS_ON);
            if (map != null) {
              Double x = getX((Double) map.get("x"));
              Double y = getY((Double) map.get("y"));
              Double w = getW((Double) map.get("w"));
              Double h = getH((Double) map.get("h"));
              String fill = (String) map.get("fill");
              String outline = (String) map.get("outline");
              Double outline_w = (Double) map.get("outline_width");
              Ellipse2D.Double ellipse = getEllipse(x, y, w, h);
              try {
                g2_.setPaint(parseColor(fill));
              }
              catch(Exception e) {
                g2_.setPaint(Color.cyan);
              }
              g2_.fill(ellipse);
              try {
                g2_.setPaint(parseColor(outline));
              }
              catch(Exception e) {
                g2_.setPaint(Color.red);
              }
              g2_.setStroke(new BasicStroke((int)outline_w.doubleValue()));
              g2_.draw(ellipse);
              g2_.drawString(label, x.intValue()-w.intValue(),
                  y.intValue()-h.intValue());
            }
            return true;
        }
        public boolean visitAfter(Graph g, Node v, Node parent) {
            return true;
        }
    };
    graph_.traverseDFS(new GraphDFSVisitor(g2));
  }

  public void init() {
    findBoundaries();
    // Increasing the window by 10%
    zoom(0.1);
    // Adjusting the aspect ratio
    adjustAspect();
  }

  public BufferedImage getImage() {
    BufferedImage image = gc_.createCompatibleImage(imageWidth_,
        imageHeight_, Transparency.OPAQUE);
    Graphics2D g2 = image.createGraphics();
    g2.setBackground(Color.white);
    g2.clearRect(0,0, imageWidth_, imageHeight_);
    paintImage(g2);
    return image;
  }

  public void prepareImage() {
    image_ = getImage();
  }

  public void drawString(String label, int x, int y) {
    Graphics2D g2 = image_.createGraphics();
    g2.setPaint(Color.blue);
    /*
    GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
    String[] names = ge.getAvailableFontFamilyNames();
    for ( int i=0; i<names.length; i++ ) {
      System.out.println( names[i] );
    }
    */
    Font f = new Font("Courier New", Font.BOLD,  22);
    FontRenderContext frc = g2.getFontRenderContext();
    g2.setFont(f);
    g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, // Anti-alias!
        RenderingHints.VALUE_ANTIALIAS_ON);
    g2.drawString(label, x, y);
  }

  public void writeImage(String file) throws IOException {
    ImageIO.write(image_, "jpeg", new File(file));
  }

  public static void main(String[] arg) throws Exception {
    String file = arg[0];
    Vector<Graph> graphs = GML.parseGMLFile(file);
    Graph graph = graphs.get(0);
    NetworkImage image = new NetworkImage(graph);
    image.prepareImage();
    image.writeImage(arg[1]);
  }

}
