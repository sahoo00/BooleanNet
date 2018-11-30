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
import java.awt.geom.*;
import java.awt.font.*;
import javax.imageio.*;

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

public class DocumentImage {

  Data data;
  GraphicsConfiguration gc;
  int headerSize;
  double xwidth;
  double ywidth;
  double border;
  double contrast;
  int imageWidth;
  int imageHeight;
  Color high, low, zero, missing;

  public DocumentImage(Data d) {
    GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
    GraphicsDevice gs = ge.getDefaultScreenDevice();
    gc = gs.getDefaultConfiguration();
    this.data = d;
    xwidth = 20; 
    ywidth = 10; 
    border = 5;
    contrast = 2;
    high = Color.red;
    low = Color.green;
    zero = Color.black;
    missing = Color.lightGray;
    int width = 20;
    int height = 20;
    BufferedImage image = gc.createCompatibleImage(width, height,
        Transparency.BITMASK);
    Graphics2D g2 = image.createGraphics();
    updateWidth(g2);
    updateHeight(g2);
    headerSize = getHeaderSize(g2);
  }

  public int getWidth() { return imageWidth; }
  public int getHeight() { return imageHeight; }

  public void update(SMSettings s) {
    Double xw = (Double) s.getOption("xwidth.fixed");
    Double yw = (Double) s.getOption("ywidth.fixed");
    Double bd = (Double) s.getOption("border");
    Double ct = (Double) s.getOption("contrast");
    xwidth = xw.doubleValue();
    ywidth = yw.doubleValue();
    border = bd.doubleValue();
    contrast = ct.doubleValue();
    high = (Color)  s.getOption("color.high");
    low = (Color)  s.getOption("color.low");
    zero = (Color)  s.getOption("color.zero");
    missing = (Color)  s.getOption("color.missing");
    int width = 20;
    int height = 20;
    BufferedImage image = gc.createCompatibleImage(width, height,
        Transparency.BITMASK);
    Graphics2D g2 = image.createGraphics();
    updateWidth(g2);
    updateHeight(g2);
    headerSize = getHeaderSize(g2);
    Boolean yfill = (Boolean) s.getOption("ywidth.fill");
    if (yfill.booleanValue()) {
       Integer h = (Integer) s.getOption("height"); 
       setMaximumHeight(h.intValue());
    }
    Boolean xfill = (Boolean) s.getOption("xwidth.fill");
    if (xfill.booleanValue()) {
       Integer w = (Integer) s.getOption("width"); 
       setMaximumWidth(w.intValue());
    }
  }

  public void updateWidth(Graphics2D g2) {
    Double[] times = data.getTimepoints();
    imageWidth = (int) (xwidth * times.length +  4 * border);
  }

  public void updateHeight(Graphics2D g2) {
    int width = (int) (headerSize +  4 * border);
    int gheight = (int) (data.getNumGenes() * ywidth + width + 3 * border);
    imageHeight = gheight;
  }

  public int getHeaderSize(Graphics2D g2) {
    Double[] times = data.getTimepoints();
    String[] header = new String[times.length];
    int maxSize = 0;
    Font f = new Font("Courior", 0, 8);
    FontRenderContext frc = g2.getFontRenderContext();
    for (int i =0; i < times.length; i++) {
      header[i] = "";
      if (times[i] != null) {
        header[i] = "" + times[i];
      }
      Rectangle2D rect = f.getStringBounds(header[i], frc);
      int length = (int) rect.getWidth();
      if (maxSize < length) {
        maxSize = length;
      }
    }
    return maxSize;
  }

  public void paintHeatmap(Graphics2D g2) {
    Double[] times = data.getTimepoints();
    String[] header = new String[times.length];
    for (int i =0; i < times.length; i++) {
      header[i] = "";
      if (times[i] != null) {
        header[i] = "" + times[i];
      }
    }
    Font f = new Font("Courior", 0, 8);
    FontRenderContext frc = g2.getFontRenderContext();

    // Printing Header
    int width = (int) (headerSize +  4 * border);
    int height = (int) (xwidth * times.length +  2 * border);
    System.out.println("Size : " + width + " x " + height);
    // Create an image that does not support transparency
    BufferedImage image = gc.createCompatibleImage(height, width,
        Transparency.BITMASK);
    Graphics2D g2d = image.createGraphics();
    g2d.translate(0, width-1);
    g2d.rotate(Math.toRadians(-90));
    g2d.setPaint(Color.blue);
    Rectangle2D rect = f.getStringBounds("1", frc);
    int h = (int) rect.getHeight();
    // Graphics gh = image.getGraphics();
    for (int i =0; i < times.length; i++) {
      int y = (int) (i * xwidth + border + xwidth/2 + h/2);
      g2d.drawString(header[i], (int) border, y);
    }
    g2d.draw(new Rectangle(0,0, width-1, height-1));
    g2.setPaint(Color.blue);
    g2.drawImage(image, (int) border, (int) border, null);

    // Printing Genes
    int gwidth = height;
    int gheight = (int) (data.getNumGenes() * ywidth);
    System.out.println("Size : " + gwidth + " x " + gheight);
    g2.draw(new Rectangle((int)border, (int) (width + 2 * border), gwidth-1, gheight-1));
    int startA = data.getNumArrayHeader(); // Array Index
    int endA = data.getNumColumns() - 1; // Array Index
    int head = data.getNumGeneHeader(); // Column index
    int end = data.getNumRows(); // Column index
    for (int i = head ; i < end; i++) {
      GeneData gene = data.getGeneData(i);
      Object[] geneData = gene.getData();
      int y = (int) ((i-head) * ywidth + width + 2 * border);
      int yn = (int) ((i+1-head) * ywidth + width + 2 * border);
      int h1 = yn - y;
      if (h1 <= 0) h1 = 1;
      for (int j = startA; j <= endA; j++) {
        Double entry = (Double)geneData[j];
        Color fill = getColor(entry);
        g2.setColor(fill);
        int x = (int) ((j-startA) * xwidth + 2 * border);
        int xn = (int) ((j+1-startA) * xwidth + 2 * border);
        int w1 = xn - x;
        if (w1 <= 0) w1 = 1;
        // System.out.println(x + "x" + y + "+" + w1 + "+" + h1);
        g2.fillRect(x, y, w1, h1);
      }
    }
  }

  public Color getColor(Color z, Color h, float val) {
    int red = (int) (h.getRed() * val + z.getRed() * ( 1 - val));
    int blue = (int) (h.getBlue() * val + z.getBlue() * ( 1 - val));
    int green = (int) (h.getGreen() * val + z.getGreen() * ( 1 - val));
    return new Color(red, green, blue);
  }

  public Color getColor(Double entry) {
    if (entry == null) {
        return missing;
    }
    float val = (float) entry.doubleValue();
    val = val / (float) contrast;
    if (val > 1) {
        val = 1;
    }
    if (val < -1) {
        val = -1;
    }
    if (val > 0) {
      return getColor(zero, high, val);
    }
    else {
      return getColor(zero, low, -val);
    }
  }

  public void setMaximumHeight(int height) {
    if (imageHeight > height) {
      int width = (int) (headerSize +  4 * border);
      ywidth = (height - width - 3 * border) / data.getNumGenes();
      imageHeight = height;
    }
  }

  public void setMaximumWidth(int width) {
    if (imageWidth > width) {
      Double[] times = data.getTimepoints();
      xwidth = (width - 4 * border) / times.length;
      imageWidth = (int) (xwidth * times.length +  4 * border);
    }
  }

  public BufferedImage getImage() {
    BufferedImage image = gc.createCompatibleImage(imageWidth, imageHeight,
        Transparency.BITMASK);
    Graphics2D g2 = image.createGraphics();
    paintHeatmap(g2);
    return image;
  }

  public void save(String file) throws IOException {
    BufferedImage img = getImage();
    ImageIO.write(img, "png", new File(file));
  }

  public static void main(String[] args) throws Exception {
    System.out.println("DocumentImage : ");
    Data data = PCLFileReader.readFile(args[0]);
    data.convertDoubles();
    DocumentImage im = new DocumentImage(data);
    im.setMaximumHeight(800);
    im.save("test.png");
  }

}


