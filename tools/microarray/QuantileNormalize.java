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

package tools.microarray;

import java.util.*;
import java.io.*;

import tools.io.*;
import tools.microarray.FileReader.*;
import tools.microarray.FileWriter.*;

public class QuantileNormalize {

    LinkedList<String> pclfiles_;
    String tmpfile_;
    String tmpfile_mean_;
    String outfile_;

    public static int BLOCK = 50;
    public static final short COOKIE = 0x11;
    public static final short PTRSPACE = 10;

    RandomAccessFile tmpfile_handle_;    
    RandomAccessFile tmpfile_mean_handle_;    
    BinaryFile tmpfile_obj_;    
    BinaryFile tmpfile_mean_obj_;    

    GeneData header_;
    LinkedList<GeneData> platform_;
    HashMap<String, Integer> idHash_; 
    long master_pointer_;
    int current_index_;
    double[] total_sum_;

    Vector<GeneData> buffer_;
    int numArrays_;

    public QuantileNormalize() {
        pclfiles_ = null;
        tmpfile_ = "tmpint";
        tmpfile_mean_ = "tmpint.mean";
        outfile_ = null;
    }

    public void setPclFiles(LinkedList<String> l) { pclfiles_ = l; }
    public void setTmpFile(String l) { tmpfile_ = l; tmpfile_mean_ = l+".mean";}
    public void setOutFile(String l) { outfile_ = l; }

    void init_io_() throws IOException {
       tmpfile_handle_ = new RandomAccessFile(tmpfile_,"rw"); 
       tmpfile_mean_handle_ = new RandomAccessFile(tmpfile_mean_,"rw"); 
       tmpfile_obj_ = new BinaryFile(tmpfile_handle_);
       tmpfile_mean_obj_ = new BinaryFile(tmpfile_mean_handle_);
    }

    void close() throws IOException {
       tmpfile_handle_.close();
       tmpfile_mean_handle_.close(); 
    }

    GeneData getHeader() throws IOException {
        GeneData header = null;
        if (pclfiles_ == null || pclfiles_.size() <= 0) {
          return header;
        }
        ListIterator<String> itr = pclfiles_.listIterator();
        while (itr.hasNext()) {
          String file = itr.next();
          PCLFileReader data = new PCLFileReader(file);
          data.begin();
          int start = data.getNumArrayHeader();
          int end = data.getNumColumns()-1;
          GeneData head = data.getHeader();
          if (header == null) {
            header = head.subset(0, start-1);
          }
          header = GeneData.merge(header, head.subset(start, end));
          data.close();
        }
        return header;
    }

    LinkedList<GeneData> getPlatform() throws IOException {
      LinkedList<GeneData> list = null;
      if (pclfiles_ == null || pclfiles_.size() <= 0) {
        return list;
      }
      String file = pclfiles_.getFirst();
      PCLFileReader data = new PCLFileReader(file);
      data.begin();
      int start = data.getNumArrayHeader();
      int end = data.getNumColumns()-1;
      GeneData head = data.getHeader();
      list = new LinkedList<GeneData>();
      list.add(head.subset(0, start-1));
      while (data.hasNext()) {
        GeneData gene = data.getData();
        // gene.print();
        if (gene == null) {
          break;
        }
        list.add(gene.subset(0, start-1));
      }
      data.close();
      return list;
    }

    HashMap<String, Integer> getIDHash(LinkedList<GeneData> platform) {
      HashMap<String, Integer> map = new HashMap<String, Integer>();
      if (platform == null) {
        return map;
      }
      ListIterator<GeneData> itr = platform.listIterator();
      int index = 0;
      while (itr.hasNext()) {
        GeneData gene = itr.next();
        String id = (String) gene.getDataAt(0);
        //System.out.println(id);
        map.put(id, new Integer(index));
        index++;
      }
      return map; 
    }

    void checkTmpHeader() throws Exception {
      tmpfile_handle_.seek(0);
      short cookie = 0;
      try {
        cookie = tmpfile_obj_.readByte();
      }
      catch(EOFException e) {
      }
      if (cookie == COOKIE) { // tmpint file present
        long numColumns = tmpfile_obj_.readDWord();
        if (numColumns != header_.size()) {
            throw new Exception("Number of columns is wrong in the intermediate file");
        }
        long numRows = tmpfile_obj_.readDWord();
        if (numRows != platform_.size()) {
            throw new Exception("Number of columns is wrong in the intermediate file");
        }
        short psize = tmpfile_obj_.readByte();
        if(psize != platform_.getFirst().size()) {
            throw new Exception("psize is wrong in the intermediate file");
        }
        master_pointer_ = tmpfile_handle_.getFilePointer();
        for (int i =0; i < PTRSPACE; i++) {
          tmpfile_obj_.readDWord(); // pointers
        }
        for (int i =0; i < numColumns; i++) {
          String item1 = tmpfile_obj_.readLengthPrefixString();
          String item2 = (String)  header_.getDataAt(i);
          if (!item1.equals(item2)) {
            throw new Exception("Header[" + i + "] - " + item1 + " != " + item2);
          }
        }
        for (int j =0; j < psize; j++) {
          short tmp = tmpfile_obj_.readByte();
        }
        for (int i =psize; i < numColumns; i++) {
          short tmp = tmpfile_obj_.readByte();
          if (tmp == 1) {
            System.out.println("Found (" +  header_.getDataAt(i) + ")\n");
            current_index_ = i - psize + 1;
          }
        }
        for (int i =0; i < numRows; i++) {
          for (int j =0; j < psize; j++) {
            String item1 = tmpfile_obj_.readLengthPrefixString();
          }
        }
      }
      else { // tmpint file not present create it
        tmpfile_handle_.seek(0);
        tmpfile_obj_.writeByte(COOKIE);
        tmpfile_obj_.writeDWord(header_.size()); // Num columns
        tmpfile_obj_.writeDWord(platform_.size()); // Num rows
        tmpfile_obj_.writeByte((short)platform_.getFirst().size()); // psize
        master_pointer_ = tmpfile_handle_.getFilePointer();
        for (int i =0; i < PTRSPACE; i++) {
          tmpfile_obj_.writeDWord(0); // pointers
        }
        writeHeaderAddress1();
        // Write header strings
        for (int i =0; i < header_.size(); i++) {
          String item = (String)  header_.getDataAt(i);
          tmpfile_obj_.writeLengthPrefixString(item);
        }
        writeHeaderAddress2();
        // Write header strings (completed)
        for (int i =0; i < header_.size(); i++) {
          tmpfile_obj_.writeByte((short)0x00);
        }
        writePlatformAddress();
        ListIterator<GeneData> itr = platform_.listIterator(); 
        while (itr.hasNext()) {
          GeneData gene = itr.next();
          for (int i =0; i < gene.size(); i++) {
            tmpfile_obj_.writeLengthPrefixString((String)gene.getDataAt(i));
          }
        }
        writeDataAddress();
        current_index_ = 0;
      }
    }

    void writeHeaderAddress1() throws IOException {
      long pos = tmpfile_handle_.getFilePointer();
      tmpfile_handle_.seek(master_pointer_ + 0 * 4);
      tmpfile_obj_.writeDWord(pos);
      tmpfile_handle_.seek(pos);
    }

    void writeHeaderAddress2() throws IOException {
      long pos = tmpfile_handle_.getFilePointer();
      tmpfile_handle_.seek(master_pointer_ + 1 * 4);
      tmpfile_obj_.writeDWord(pos);
      tmpfile_handle_.seek(pos);
    }

    void writePlatformAddress() throws IOException {
      long pos = tmpfile_handle_.getFilePointer();
      tmpfile_handle_.seek(master_pointer_ + 2 * 4);
      tmpfile_obj_.writeDWord(pos);
      tmpfile_handle_.seek(pos);
    }

    void writeDataAddress() throws IOException {
      long pos = tmpfile_handle_.getFilePointer();
      tmpfile_handle_.seek(master_pointer_ + 3 * 4);
      tmpfile_obj_.writeDWord(pos);
      tmpfile_handle_.seek(pos);
    }

    void loadMeanArray() throws IOException {
      total_sum_ = new double[platform_.size()];
      for(int i =0; i < total_sum_.length; i++) {
        total_sum_[i] = 0.0;
      }
      if (current_index_ <= 0) {
        return;
      }
      tmpfile_mean_handle_.seek(0);
      for (int i =0; i < total_sum_.length; i++) {
        long val = tmpfile_mean_obj_.read64Word();
        total_sum_[i] = Double.longBitsToDouble(val);
        // System.out.println(total_sum_[i] + " " + val);
      }
    }

    void writeMeanArray() throws IOException {
      tmpfile_mean_handle_.seek(0);
      for (int i =0; i < total_sum_.length; i++) {
        long val = Double.doubleToLongBits(total_sum_[i]);
        tmpfile_mean_obj_.write64Word(val);
        // System.out.println(total_sum_[i] + " " + val);
      }
    }

    void initBuffer() {
      Vector<GeneData> res = new Vector<GeneData>();
      for (int i =0; i < platform_.size(); i++) {
        Object[] o = new Object[BLOCK];
        res.add(new GeneData(o));
      }
      buffer_ = res;
      numArrays_ = 0;
    }

    boolean hasNextBlock() {
      return current_index_ < (header_.size() - platform_.getFirst().size());
    }

    void fillNextBlock() throws Exception {
      int numArrays = (header_.size() - platform_.getFirst().size()) - current_index_ ;
      if (numArrays > BLOCK) {
        numArrays = BLOCK;
      }
      numArrays_ = numArrays;
      int last = current_index_ + numArrays;
      int currentNum = 0;
      ListIterator<String> itr = pclfiles_.listIterator();
      while (itr.hasNext()) {
        String file = itr.next();
        PCLFileReader data = new PCLFileReader(file);
        data.begin();
        int num = data.getNumArrays();
        // System.out.println(current_index_+" "+currentNum+" "+last+" "+num);
        if (current_index_ <= (currentNum+num) && current_index_ < last) {
          int nH = data.getNumArrayHeader();
          int start = nH + current_index_ - currentNum;
          int end = data.getNumColumns() - 1;
          if ((currentNum + end - nH) >= last) {
            end = last - currentNum + nH - 1;
          }
          // System.out.println(start+" "+end);
          GeneData head = data.getHeader();
          for (int i = start ; i <= end; i++) {
            // System.out.print(head.getDataAt(i) + " ");
            buffer_.get(0).setDataAt(
                i+currentNum-nH-current_index_, head.getDataAt(i));
            buffer_.get(1).setDataAt(
                i+currentNum-nH-current_index_, "1");
          }
          // System.out.println();
          while (data.hasNext()) {
            long lineno = data.getLineNumber();
            GeneData gene = data.getData();
            if (gene == null) {
              break;
            }
            // gene.print();
            String id = (String) gene.getDataAt(0);
            if (idHash_.containsKey(id)) {
              gene.convertDouble(start, end);
              Integer loc = idHash_.get(id);
              for (int i = start; i <= end; i++) {
                // System.out.print(gene.getDataAt(i) + " ");
                buffer_.get(loc.intValue()).setDataAt(
                    i+currentNum-nH-current_index_, gene.getDataAt(i));
              }
              // System.out.println();
            }
          }
        }
        currentNum += data.getNumArrays();
        data.close();
      }
      current_index_ = last;
    }

    void setDataAddress(int currentIndex) throws IOException {
      tmpfile_handle_.seek(master_pointer_ + 3 * 4);
      long pos = tmpfile_obj_.readDWord();
      tmpfile_handle_.seek(pos + currentIndex * platform_.size() * 4);
    }

    public void saveRanks(int currentIndex, Integer[] rank) throws IOException {
      tmpfile_obj_.writeDWord(currentIndex); // headers
      tmpfile_obj_.writeDWord(0);
      for (int j =0; j < rank.length ; j++) {
        tmpfile_obj_.writeDWord(rank[j].longValue());
      }
    }

    public void updateMean(Integer[] rank, int num) throws IOException {
      for (int i = 0; i < rank.length; i++) {
        Double c1 = (Double) buffer_.get(rank[i].intValue()+2).getDataAt(num);
        if (c1 == null) {
          continue;
        }
        total_sum_[i+2] += c1.doubleValue();
      }
      // writeMeanArray();
    }

    void setHeaderAddress2(int index) throws IOException {
      int psize = platform_.getFirst().size();
      tmpfile_handle_.seek(master_pointer_ + 1 * 4);
      long pos = tmpfile_obj_.readDWord();
      tmpfile_handle_.seek(pos + index + psize);
      tmpfile_obj_.writeByte((short)0x1);
    }

    public void computeRanks(int i, int currentIndex) throws IOException {
      System.out.println("num -> " + (currentIndex + i));
      Integer[] rank = new Integer[buffer_.size()-2];
      for (int j =0; j < rank.length ; j++) {
        rank[j] = new Integer(j);
      }
      class RankComparator implements Comparator<Integer> {
        int num_;
        Vector<GeneData> data_;

        public RankComparator(int n,  Vector<GeneData> d) {
            num_ = n;
            data_ = d;
        }
        public int compare(Integer s1, Integer s2) {
          //System.out.println(s1 + " " + num_ + " " + data_.get(s1.intValue() + 2).getDataAt(num_));
          Double c1 = (Double) data_.get(s1.intValue() + 2).getDataAt(num_);
          Double c2 = (Double) data_.get(s2.intValue() + 2).getDataAt(num_);
          if (c1 == null) {
            return -1;
          }
          if (c2 == null) {
            return 1;
          }
          if (c1.doubleValue() > c2.doubleValue()) {
            return 1;
          }
          return -1;
        }
      };
      Arrays.sort(rank, new RankComparator(i, buffer_));
      updateMean(rank, i);
      Integer[] rankOrder = new Integer[buffer_.size()-2];
      int index = 0;
      while (index < rankOrder.length) {
        int j = index;
        while (j < (rankOrder.length - 1)) {
            Double c1 = (Double) buffer_.get(rank[j]+2).getDataAt(i);
            Double c2 = (Double) buffer_.get(rank[j+1]+2).getDataAt(i);
            if (c1 == null && c2 != null) {
                break;
            }
            if (c1 != null && c2 == null) {
                break;
            }
            if (c1 != null && !c1.equals(c2)) {
              break;
            }
            j++;
        }
        // System.out.println(index + " " + j);
        if (index != j) {
            for (int k = index; k <= j; k++) {
                rankOrder[rank[k].intValue()] = new Integer((int)Math.floor((index + j + 2)/2.0));
            }
        }
        else {
          rankOrder[rank[index].intValue()] = new Integer(index + 1);
        }
        index = j + 1;
      }
      /*
      for (int j = 0; j < rank.length; j++) {
        System.out.println(rank[j]);
      }
      for (int j = 0; j < rankOrder.length; j++) {
        System.out.println(rankOrder[j]);
      }
      */
      setDataAddress(currentIndex+i);
      saveRanks(currentIndex+i, rankOrder);
      setHeaderAddress2(currentIndex+i);
    }

    public void quantileNormalize() throws Exception {
       if (pclfiles_ == null || pclfiles_.size() <= 0) {
          return;
       }
       init_io_();
       header_ = getHeader();
       platform_ = getPlatform();
       idHash_ = getIDHash(platform_);
       checkTmpHeader();
       loadMeanArray();
       System.out.println(header_.size());
       System.out.println(platform_.size());
       System.out.println(current_index_);
       while (hasNextBlock()) {
         int index = current_index_;
         initBuffer();
         fillNextBlock();
         for (int i =0; i < numArrays_; i++) {
            computeRanks(i, index);
         }
       }
       writeMeanArray();
    }

    public void writePCLfile() throws Exception {
      PCLFileWriter writer = new PCLFileWriter(outfile_);
      writer.writeData(header_);
      int numArrays = header_.size() - platform_.getFirst().size();
      Object[] weight = new Object[header_.size()];
      for (int i =0; i < weight.length; i++) {
        weight[i] = "1";
      }
      weight[0] = "EWEIGHT"; weight[1] = ""; weight[1] = "";
      writer.writeData(new GeneData(weight));
      for (int i = 2; i < platform_.size(); i++) {
        GeneData h = platform_.get(i);
        Object[] d = new Object[header_.size()];
        for (int j =0; j < h.size(); j++) {
          d[j] = h.getDataAt(j);
        }
        for (int j = h.size(); j < header_.size(); j++) {
          tmpfile_handle_.seek(master_pointer_ + 3 * 4);
          long pos = tmpfile_obj_.readDWord();
          tmpfile_handle_.seek(pos + (j-h.size()) * platform_.size() * 4 + i*4);
          int rank = (int) tmpfile_obj_.readDWord();
          d[j] = new Double(total_sum_[rank+1]/numArrays);
        }
        writer.writeData(new GeneData(d));
      }
      writer.close();
    }

    public static void main(String arg[]) throws Exception {
        QuantileNormalize q = new QuantileNormalize();
        LinkedList<String> list = new LinkedList<String>(Arrays.asList(arg));
        String tmpfile = list.removeFirst();
        String outfile = list.removeFirst();
        QuantileNormalize.BLOCK = Integer.parseInt(list.removeFirst());
        q.setPclFiles(list);
        q.setTmpFile(tmpfile);
        q.setOutFile(outfile);
        q.quantileNormalize();
        q.writePCLfile();
    }
}

