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

package tools.microarray;

import java.text.MessageFormat;

public class Matrix {
  protected int numRows;
  protected int numCols;

  public double[][] cells;

  public Matrix(int rows, int cols) {
    numRows= rows; numCols= cols;
    cells= new double[numRows][numCols];
    for(int i=0; i < numRows; i++) {
      for(int j=0; j < numCols; j++) {
        cells[i][j]= 0.0;
      }
    }
  }

  /**
   * Copy constructor
   */
  public Matrix(Matrix m) {
    numRows= m.numRows;
    numCols= m.numCols;
    cells = new double[numRows][numCols];
    for(int i=0; i < numRows; i++) {
      for(int j=0; j < numCols; j++) {
        this.cells[i][j]= m.cells[i][j];
      }
    }
  }

  public int getNumRows() { return numRows; }
  public int getNumCols() { return numCols; }
  public int getNumCells() { return numCols * numRows; }
  public void setCell(int i, int j, double n) {cells[i][j] = n;}
  public double getCell(int i, int j) { return cells[i][j];}
  public double[] getRow(int i) { return cells[i];}

  // Averages of the all the cells
  public double mean() {
    double m = 0;
    for(int i=0; i < numRows; i++) {
      for(int j=0; j < numCols; j++) {
        m = m + this.cells[i][j];
      }
    }
    m = m / numRows / numCols;
    return m;
  }

  // Mean square error of the all the cells
  public double mse() {
    double m = mean();
    double var = 0;
    for(int i=0; i < numRows; i++) {
      for(int j=0; j < numCols; j++) {
        var = var + (this.cells[i][j] - m) * (this.cells[i][j] - m);
      }
    }
    return var;
  }

  // Variance of the all the cells
  public double var() {
    double var = mse();
    var = var / (numRows * numCols - 1);
    return var;
  }

  // Standard deviation of the all the cells
  public double sd() {
    double v = var();
    return Math.sqrt(v);
  }

  /**
   * @param tx ty tz translations in the x,y, and z directions
   * @returns the associated translation transformation matrix
   */
  public static Matrix translation(double tx, double ty, double tz) {
    Matrix T= new Matrix(4, 4);
    T.cells[0][0]= 1;
    T.cells[1][1]= 1;
    T.cells[2][2]= 1;
    T.cells[3][3]= 1;
    T.cells[3][0]= tx;
    T.cells[3][1]= ty;
    T.cells[3][2]= tz;
    return T;
  }

  /**
   * @param sx sy sz translations in the x,y, and z directions
   * @returns the associated scaling transformation matrix
   */
  public static Matrix scaling(double sx, double sy, double sz) {
    Matrix S= new Matrix(4,4);
    S.cells[0][0]= sx;
    S.cells[1][1]= sy;
    S.cells[2][2]= sz;
    S.cells[3][3]= 1;
    return S;
  }

  /**
   * @param theta an angle in radians
   * @return the associated x-axis rotation transformation matrix
   */
  public static Matrix xRotation(double theta) {
    Matrix R= new Matrix(4, 4);
    double c= Math.cos(theta);
    double s= Math.sin(theta);
    R.cells[0][0]= 1;
    R.cells[1][1]= c;
    R.cells[2][2]= c;
    R.cells[3][3]= 1;
    R.cells[1][2]= s;
    R.cells[2][1]= -s;
    return R;
  }

  /**
   * @param theta an angle in radians
   * @return the associated y-axis rotation transformation matrix
   */
  public static Matrix yRotation(double theta) {
    Matrix R= new Matrix(4, 4);
    double c= Math.cos(theta);
    double s= Math.sin(theta);
    R.cells[0][0]= c;
    R.cells[1][1]= 1;
    R.cells[2][2]= c;
    R.cells[3][3]= 1;
    R.cells[2][0]= s;
    R.cells[0][2]= -s;
    return R;
  }

  /**
   * @param theta an angle in radians
   * @return the associated z-axis rotation transformation matrix
   */
  public static Matrix zRotation(double theta) {
    Matrix R= new Matrix(4, 4);
    double c= Math.cos(theta);
    double s= Math.sin(theta);
    R.cells[0][0]= c;
    R.cells[1][1]= c;
    R.cells[2][2]= 1;
    R.cells[3][3]= 1;
    R.cells[0][1]= s;
    R.cells[1][0]= -s;
    return R;
  }

  /**
   * @param m a Matrix
   * @return a new matrix which is equal to the sum of this + m
   */
  public Matrix addTo (Matrix m) {
    if(numRows!=m.numRows || numCols!=m.numCols) {
      System.out.println("addTo: dimensions are different");
      System.exit(1);
    }
    Matrix ret= new Matrix(numRows, numCols);

    for(int i=0; i < numRows; i++) {
      for(int j=0; j < numCols; j++) {
        ret.cells[i][j]= this.cells[i][j] + m.cells[i][j];
      }
    }

    return ret;
  }

  /**
   * @param m a Matrix
   * @return a new matrix which is equal to the difference of this - m
   */
  public Matrix subtractFrom (Matrix m) {
    if(numRows!=m.numRows || numCols!=m.numCols) {
      System.out.println("subtractFrom: dimensions are different");
      System.exit(1);
    }
    Matrix ret= new Matrix(numRows, numCols);

    for(int i=0; i < numRows; i++) {
      for(int j=0; j < numCols; j++) {
        ret.cells[i][j]= this.cells[i][j] - m.cells[i][j];
      }
    }

    return ret;
  }

  /**
   * @param m a Matrix
   * @return a new matrix which is equal to the product of this*m
   */
  public Matrix multiply (Matrix m) {
    if(numCols!=m.numRows) {
      System.out.println("dimensions bad in multiply()");
      System.exit(1);
    }

    Matrix ret= new Matrix(numRows, m.numCols);

    for(int i=0; i < ret.numRows; i++) {
      for(int j=0; j < ret.numCols; j++) {
	for(int k=0; k < this.numCols; k++) {
          ret.cells[i][j] += this.cells[i][k] * m.cells[k][j];
        }
      }
    }

    return ret;
  }

  /**
   * @param m a Matrix
   * @return a new matrix which is equal to the weighted Average
   */
  public Matrix weightedAverage(Matrix m) {
    if(numRows != m.numRows && m.numCols != 1) {
      System.out.println("dimensions bad in weightedAverage");
      System.exit(1);
    }

    Matrix up = new Matrix(1, numCols);
    Matrix down = new Matrix(1, numCols);

    boolean foundZero = false;
    for(int i=0; i < numRows; i++) {
        if (m.cells[i][0] == 0) {
            foundZero = true;
            break;
        }
    }
 
    for(int i=0; i < numRows; i++) {
      for(int j=0; j < numCols; j++) {
        if (foundZero) {
          up.cells[0][j] += this.cells[i][j];
          down.cells[0][j] += 1;
        }
        else {
          up.cells[0][j] += this.cells[i][j] / m.cells[i][0];
          down.cells[0][j] += 1 / m.cells[i][0];
        }
      }
    }

    Matrix ret = new Matrix(numCols, 1);
    for(int i=0; i < numCols; i++) {
      ret.cells[i][0] = up.cells[0][i] / down.cells[0][i];
    }

    return ret;
  }

  /**
   * Scalar multiplication- multiplies each element by a scalar
   * @param s a scalar
   * @return a new matrix which is equal to the product of s*this
   */
  public Matrix multiply (double s) {
    Matrix ret= new Matrix(numRows, numCols);
    for(int i=0; i<numRows; i++)
      for(int j=0; j<numCols; j++)
	ret.cells[i][j]= cells[i][j]*s;
    return ret;
  }

  /**
   * @return the transposed matrix with dimensions numCols x numRows
   */
  public Matrix transpose() {
    Matrix ret= new Matrix(numCols, numRows);

    for(int r=0; r<numRows; r++)
      for(int c=0;c<numCols; c++)
	ret.cells[c][r]= cells[r][c];
    return ret;
  }

  /**
   * @param val a scalar
   * @returns true if and only if all elements of the matrix equal val
   */
  public boolean equals(double val) {
    for(int i=0; i<numRows; i++)
      for(int j=0; j<numCols; j++)
	if(Math.abs(cells[i][j]-val) > .0001) return false;
    return true;
  }

  /**
   * Computes the dot product (or scalar product) of two matrices by
   *  multiplying corresponding elements and summing all the products.
   * @param m A Matrix with the same dimensions
   * @returns the dot product (scalar product)
   */
  public double dot(Matrix m) {
    if(numRows!=m.numRows || numCols!=m.numCols) {
      System.out.println("dot: dimensions different"); 
      System.exit(1);
    }
    double sum= 0;

    for(int r=0; r<numRows; r++)
      for(int c=0; c<numCols; c++)
        sum += this.cells[r][c] * m.cells[r][c];

    return sum;
  }

  /**
   * Calculates the matrix's Moore-Penrose pseudoinverse
   * @return an MxN matrix which is the matrix's pseudoinverse.
   */
  public Matrix pseudoInverse() {

    int r,c;

    int k=1;
    Matrix ak= new Matrix(numRows, 1);
    Matrix dk, ck, bk;

    Matrix R_plus;

    for(r=0; r<numRows; r++)
      ak.cells[r][0]= this.cells[r][0];

    if(!ak.equals(0.0)) {
      R_plus= ak.transpose().multiply( 1.0/( ak.dot(ak) ) );
    }
    else {
      R_plus= new Matrix(1, numCols);
    }

    while(k< this.numCols) {

      for(r=0; r<numRows; r++)
        ak.cells[r][0]= this.cells[r][k];

      dk= R_plus.multiply(ak);
      Matrix T= new Matrix(numRows, k);
      for(r=0; r<numRows; r++)
        for(c=0; c<k; c++)
          T.cells[r][c]= this.cells[r][c];

      ck= ak.subtractFrom( T.multiply(dk) );

      if( !ck.equals(0.0) ) {
        bk= ck.transpose().multiply( 1.0/(ck.dot(ck)) );
      }
      else {
        bk= dk.transpose().multiply( 1.0/( 1.0 + dk.dot(dk) ) ).multiply(R_plus);
      }

      Matrix N= R_plus.subtractFrom( dk.multiply(bk) );
      R_plus= new Matrix(N.numRows+1, N.numCols);

      for(r=0; r< N.numRows; r++)
        for(c=0; c< N.numCols; c++)
          R_plus.cells[r][c]= N.cells[r][c];
      for(c=0; c<N.numCols; c++)
        R_plus.cells[R_plus.numRows-1][c]= bk.cells[0][c];

      k++;
    }
    return R_plus;
  }

  public static String formatString(String format, double x) {
    MessageFormat mf = new MessageFormat(
        "{0,number," + format + "}");
    Object[] objs = {new Double(x)};
    return mf.format(objs);
  }

  /**
   * @return a String representation of the matrix
   */
  public String toString() {
    StringBuffer buf= new StringBuffer();
    buf.append("[ \n");
    for(int r=0; r<numRows; r++) {
      buf.append("[ ");
      for(int c=0; c<numCols; c++) {
        buf.append(formatString("0.##", cells[r][c]));
        buf.append(" ");
      }
      buf.append("] \n");
    }
    buf.append("]\n");
    return buf.toString();
  }

  public void print() {
    System.out.println(toString());
  }
}
