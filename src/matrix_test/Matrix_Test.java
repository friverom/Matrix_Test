
package matrix_test;

/**
 * Matrix Library Test
 * @author Federico
 */
public class Matrix_Test {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        
        //Create a Matrix 3x3 of zeros
        Matrix m1=new Matrix(3,3);
        System.out.println("Create a 3x3 Matrix of zeros");
        m1.show();
        System.out.println("\n");
        
        
        //Create Matrix based on 2D array
        double [][] data={{1,2,3},{4,5,6},{7,8,9}};
        Matrix m2=new Matrix(data);
        System.out.println("Create a 3x3 matrix based on 2d Array");
        m2.show();
        System.out.println("\n");
        
        //Create a copy of a Matrix
        Matrix m3=new Matrix(m2);
        System.out.println("Create a copy of a Matrix");
        m3.show();
        System.out.println("\n");
        
        //Create a Matrix 3x3 of random numbers between -0.5 and 0.5
        Matrix m4=Matrix.random(3, 3);
        System.out.println("Create a Matrix 3x3 of random numbers between -0.5 and 0.5");
        m4.show();
        System.out.println("\n");
        
        // Create and return a random M-by-N matrix with values between -1 and 1
        //normally distributed with mean=0 and sd=1
        Matrix m5=Matrix.random_gaussian(3, 3);
        System.out.println(" Create and return a random M-by-N matrix with values between -1 and 1 \n normally distributed with mean=0 and sd=1");
        m5.show();
        System.out.println("\n");
        
        // Create and return a random M-by-N Xavier matrix with values between -1 and 1
        //normally distributed with mean=0 and sd=sqrt(1/M)
        Matrix m6=Matrix.xavier(6, 4);
        System.out.println(" Create and return a random M-by-N Xavier matrix with values between -1 and 1 \n normally distributed with mean=0 and sd=sqrt(1/M");
        m6.show();
        System.out.println("\n");

        // Create and return a random M-by-N HE matrix with values between -1 and 1
        //normally distributed with mean=0 and sd=sqrt(2/M)
        Matrix m7=Matrix.he(6, 4);
        System.out.println(" Create and return a random M-by-N HE matrix with values between -1 and 1 \n normally distributed with mean=0 and sd=sqrt(2/M");
        m7.show();
        System.out.println("\n");
        
        // Create and return an Identity matrix of NxN
        Matrix m8=Matrix.identity(6);
        System.out.println("Create and return an Identity matrix of N x N");
        m8.show();
        System.out.println("\n");
        
        //Create and returns the Transpose of Matrix
        //Transpose of matrix m3
        System.out.println("Transpose of matrix");
        System.out.println("Matrix m3");
        m3.show();
        System.out.println("\n Transpose matrix m3");
        Matrix m9=m3.transpose();
        m9.show();
        System.out.println("\n");
        
        //Sum 2 matrices
        //Sum matrix m3 + m3 transpose
        System.out.println("Sum of 2 matrices");
        System.out.println("Matrix m3");
        m3.show();
        System.out.println("\n Transpose matrix m3");
        m9.show();
        Matrix sum=m3.plus(m9);
        System.out.println("m3 + m9");
        sum.show();
        System.out.println("\n");
        
        //Substraction of 2 matrices
        //m3-m3 transpose
        System.out.println("Substraction of 2 matrices");
        System.out.println("Matrix m3");
        m3.show();
        System.out.println("\n Transpose matrix m3");
        m9.show();
        Matrix sub=m3.minus(m9);
        System.out.println("m3 - m9");
        sub.show();
        System.out.println("\n");
        
        //Matrix equality.
        Matrix m10=new Matrix(m3);
        System.out.println("Matrix m10 is equal to m3? "+m10.eq(m3));
        
        //Hadamar product.(element wise multiplication
        Matrix m11=new Matrix(new double[]{1,2,3});
        Matrix m12=new Matrix(new double[]{3,2,1});
        Matrix had=m11.times_hadamard(m12);
        System.out.println("Hadamar Product");
        System.out.println("Vector v1");
        m11.show();
        System.out.println("\n");
        System.out.println("Vector v2");
        m12.show();
        System.out.println("\n");
        System.out.println("Hadamar Product");
        had.show();
        System.out.println("\n");
        
        //Matrix multiplication
        Matrix m13=new Matrix(new double[][]{{1,2,3},{4,5,6},{7,8,9}});
        Matrix m14=new Matrix(new double[]{10,11,12});
        Matrix mul=m13.times(m14);
        System.out.println("Matrix Multiplication");
        System.out.println("Matrix A");
        m13.show();
        System.out.println("Matrix B");
        m14.show();
        System.out.println("A x B");
        mul.show();
        
        //Scalar product
        System.out.println("Matrix scaling");
        System.out.println("Matrix");
        m13.show();
        System.out.println("Scale by 2");
        Matrix sca=m13.scale(2);
        sca.show();
        
        //System of equations solver
        Matrix m15=new Matrix(new double[][]{{2,1},{-1,4}});
        Matrix m16=new Matrix(new double[]{4,7});
        System.out.println("m15 x = m16. Solve for x");
        System.out.println("Matrix m15");
        m15.show();
        System.out.println("Matrix m16");
        m16.show();
        Matrix m17=m15.solve(m16);
        System.out.println("Results:");
        m17.show();
        
        //Dot product of vectors
      //  Matrix a=new Matrix(new double[][]{{1,2,3},{4,5,6},{7,8,9}});
      //Matrix b=a.transpose();  
        Matrix a=new Matrix(new double[]{1,2,3});
        Matrix b=new Matrix(new double[]{1,2,3});
        System.out.println("Dot product of vectors");
        System.out.println("Matrix A");
        a.show();
        System.out.println("Matrix B");
        b.show();
        System.out.println("A dot B");
        Matrix c=a.dotMult(b);
        c.show();
        
        
    }
    
}
