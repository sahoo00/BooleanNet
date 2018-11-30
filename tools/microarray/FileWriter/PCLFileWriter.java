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

package tools.microarray.FileWriter;

import tools.microarray.Data;
import tools.microarray.GeneData;
import java.io.*;
import java.util.LinkedList;
import java.util.ListIterator;

public class PCLFileWriter {

  String filename_;
  BufferedWriter out_;

  public PCLFileWriter(String file) throws IOException {
    filename_ = file;
    out_ = new BufferedWriter(new FileWriter(filename_));
  }

  public void writeData(GeneData gene) throws IOException {
    if (gene != null) {
      out_.write(gene.toString());
    }
  }

  public void close() throws IOException {
    out_.close();
  }

  public static void writeFile(Data data, String filename) throws IOException {
    writeFile(data, filename, null);
  }

  public static void writeFile(Data data, String filename, int[] order) throws IOException {
    System.out.println("Writing file " + filename);
    BufferedWriter out = new BufferedWriter(new FileWriter(filename));
    int max = data.getNumRows();
    if (order != null) { // Specific gene order is given
      max = order.length;
      for (int i =0; i < data.getNumGeneHeader(); i++) {
        GeneData gene = data.getGeneData(i);
        out.write(gene.toString());
      }
    }

    for (int i =0; i < max; i++) {
      GeneData gene;
      if (order != null) {
        if (order[i] < data.getNumRows()) {
          gene = data.getGeneData(order[i]);
          if (gene == null) {
            System.err.println("Error at index" + i);
            System.exit(0);
          }
          out.write(gene.toString());
        }
      }
      else {
        gene = data.getGeneData(i);
        if (gene == null) {
            System.err.println("Error at index" + i);
            System.exit(0);
        }
        gene = gene.subset(0, data.getNumColumns()-1);
        out.write(gene.toString());
      }
    }
    out.close();
    System.out.println("Done");
  }
};

