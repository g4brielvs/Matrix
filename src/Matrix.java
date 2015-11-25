
/**
 * A bare-bones immutable data type for M-by-N matrices
 *
 * @author gabrielus
 */
final public class Matrix {
	
  /*
   * ------------------------ Class variables ------------------------
   */
	
	/**
	 * Epsilon for numerical methods
	 */
	private static final double EPSILON = 1e-10;
	
	/**
	 * Row and column dimensions.
	 */
    private final int M;  // number of rows           
    private final int N;  // number of columns        
    
    /**
     * Array for internal storage of elements.
     */
    private final double[][] data; // M-by-N array
    
    /*
     * ------------------------ Constructors ------------------------
     */
    
    /**
     * Constructor of a M-by-N matrix of 0's
     * @param M cardinality of rows
     * @param N cardinality of columns
     */
    public Matrix(int M, int N) {
        this.M = M;
        this.N = N;
        data = new double[M][N];
    }
        
    /**
     * Constructor of a matrix based on 2D array
     * @param data Two-dimensional array of doubles.
     */
    public Matrix(double[][] data) {
        M = data.length;
        N = data[0].length;
        this.data = new double[M][N];
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                    this.data[i][j] = data[i][j];
    }
    
    private Matrix(Matrix A) { this(A.data); }
    
    /*
     * ------------------------ Public Methods ------------------------
     */

    /**
     * Creates and returns a random M-by-N matrix with values between 0 and 1
     * @param M cardinality of rows
     * @param N cardinality of columns
     * @return a random M-by-N matrix with values between 0 and 1
     */
    public static Matrix random(int M, int N) {
        Matrix A = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[i][j] = Math.random();
        return A;
    }
    
    /**
     * Creates and returns the N-by-N identity matrix
     * @param order N
     * @return the N-by-N identity matrix
     */
    public static Matrix identity(int N) {
        Matrix I = new Matrix(N, N);
        for (int i = 0; i < N; i++)
            I.data[i][i] = 1;
        return I;
    }
    
    /**
     * Creates and returns the N-by-N random lower triangular matrix
     * @param order N
     * @return a N-by-N random symmetric matrix with values between 0 and 1
     */
    public static Matrix triangular(int N) {
        Matrix I = new Matrix(N, N);
        for(int i = 0; i < N; i++) {
        	for (int j = 0; j <= i; j++) {
        		I.data[i][j] = Math.random();
        	}
        }
        return I;
    }
    
    /**
     * Creates and returns the N-by-N random symmetric matrix
     * @param order N
     * @return a N-by-N random symmetric matrix with values between 0 and 1
     */
    public static Matrix symmetric(int N) {
        Matrix I = triangular(N);
        for(int i = 0; i < N; i++) {
        	for (int j = i + 1; j < N; j++) {
        		I.data[i][j] = I.data[j][i];
        	}
        }
        return I;
    }
        
    /**
     * Creates and returns the transpose of the invoking matrix
     * @return the transpose of the invoking matrix
     */
    public Matrix transpose() {
        Matrix A = new Matrix(N, M);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                A.data[j][i] = this.data[i][j];
        return A;
    }
    
    /**
     * Creates and returns A + B
     * @param B
     * @return A + B
     */
    public Matrix plus(Matrix B) {
        Matrix A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] + B.data[i][j];
        return C;
    }
    
    /**
     * Creates and returns A - B
     * @param B
     * @return A - B
     */
    public Matrix minus(Matrix B) {
        Matrix A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(M, N);
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                C.data[i][j] = A.data[i][j] - B.data[i][j];
        return C;
    }
    
    /**
     * Creates and returns A * B
     * @param B 
     * @return A * B
     */
    public Matrix times(Matrix B) {
        Matrix A = this;
        if (A.N != B.M) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(A.M, B.N);
        for (int i = 0; i < C.M; i++)
            for (int j = 0; j < C.N; j++)
                for (int k = 0; k < A.N; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }
    
    /**
     * Checks whether A equals B
     * @param B
     * @return whether A equals B
     */
    public boolean equals(Matrix B) {
        Matrix A = this;
        if (B.M != A.M || B.N != A.N) throw new RuntimeException("Illegal matrix dimensions.");
        for (int i = 0; i < M; i++)
            for (int j = 0; j < N; j++)
                if (A.data[i][j] != B.data[i][j]) return false;
        return true;
    }
    
    
    /**
     * 
     * @return whether matrix is square
     */
    public boolean isSquare() {
    	return M==N;
    }
    
    /**
     * 
     * @return whether matrix is triangular
     */
    public boolean isTriangular() {
        for (int i = 0; i < N; i++) {
            for (int j = i + 1; j < N; j++) {
                if (data[i][j] != 0.0) return false;
            }
        }
        return true;
    }
    
    /**
     * 
     * @return whether matrix is symmetric
     */
    public boolean isSymmetric() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < i; j++) {
                if (data[i][j] != data[j][i]) return false;
            }
        }
        return true;
    }
        
    /**
     * Swaps rows i and j
     * @param i i-th row
     * @param j j-ith row
     */
    private void swap(int i, int j) {
        double[] temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }
    
    /**
     * Returns x = A^-1 b, assuming A is square and has full rank
     * @param rhs Right Hand Side
     * @return x = A^-1 b
     */
    public Matrix solve(Matrix rhs) {
    	
        if (M != N || rhs.M != N || rhs.N != 1)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this);
        Matrix b = new Matrix(rhs);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < N; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < N; j++)
                if (Math.abs(A.data[j][i]) > Math.abs(A.data[max][i]))
                    max = j;
            A.swap(i, max);
            b.swap(i, max);

            // singular or nearly singular
           if (Math.abs(A.data[i][i]) <= EPSILON) {
                throw new RuntimeException("Matrix is singular or nearly singular");
            }

            // pivot within b
            for (int j = i + 1; j < N; j++)
                b.data[j][0] -= b.data[i][0] * A.data[j][i] / A.data[i][i];

            // pivot within A
            for (int j = i + 1; j < N; j++) {
                double m = A.data[j][i] / A.data[i][i];
                for (int k = i+1; k < N; k++) {
                    A.data[j][k] -= A.data[i][k] * m;
                }
                A.data[j][i] = 0.0;
            }
        }

        // back substitution
        Matrix x = new Matrix(N, 1);
        for (int j = N - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < N; k++)
                t += A.data[j][k] * x.data[k][0];
            x.data[j][0] = (b.data[j][0] - t) / A.data[j][j];
        }
        return x;
   
    }
    
    /**
     * Returns Cholesky factor L of positive-semidefinite matrix A = L L^T
     * @return Cholesky factor L of positive-semidefinite matrix A = L L^T
     * @exception RuntimeException Matrix must be positive-semidefinite
     */
    public Matrix cholesky() {
        if (!isSquare()) {
            throw new RuntimeException("Matrix is not square");
        }
        if (!isSymmetric()) {
            throw new RuntimeException("Matrix is not symmetric");
        }

        Matrix L = new Matrix(new double[N][N]);
        		
        for (int i = 0; i < N; i++)  {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += L.data[i][k] * L.data[j][k];
                }
                if (i == j) L.data[i][i] = Math.sqrt(this.data[i][i] - sum);
                else        L.data[i][j] = 1.0 / L.data[j][j] * (this.data[i][j] - sum);
            }
            if (L.data[i][i] <= 0) {
                throw new RuntimeException("Matrix not positive definite");
            }
        }
        return L;
    }
    
    /**
     * Computes the determinant of the matrix using Laplace's Formula
     * @param matrix
     * @param n
     * @return determinant of the matrix
     */
    private double determinant(double matrix[][], int n) { 
    	if (!isSquare()) {
            throw new RuntimeException("Matrix is not square");
        }
    	
    	double det = 0.0;
    	int sign = 1, p = 0, q = 0;
    	
    	if(n==1){
    		det = matrix[0][0];
    	}
    	else{
    		double factor[][] = new double[n-1][n-1];
    		for(int x = 0 ; x < n ; x++) {
    			p = 0;
    			q = 0;
    			for(int i = 1; i < n; i++) {
    				for(int j = 0; j < n; j++) {
    					if(j != x) {
    						factor[p][q++] = matrix[i][j];
    						if(q % (n - 1) == 0) {
    							p++;
    							q = 0;
    						}
    					}
    				}
    			}
    			det += sign * matrix[0][x] * determinant(factor, n-1);
    			sign = -sign;
    		}
    	}
    	return det;
    }
    
    /**
     * Computes the determinant of the matrix using Laplace's Formula
     * @return determinant of the matrix
     */
    public double determinant() {
    	return determinant(data, N);	
    }
    
	/**
     * Print matrix to standard output
     */
    public void show() {
    	
    	try {
    		int rows = data.length;
            int columns = data[0].length;
    		String str = "|\t";
    		
    		for(int i = 0; i < rows; i++) {
               for(int j = 0; j < columns; j++) {
            	   str += String.format("%9.4f", data[i][j]) + "\t";
               }
               System.out.println(str + "|");
               str = "|\t"; 
            }

       } catch(Exception e){
    	   System.out.println("Matrix is empty!!");
       }
    }
        
    //@Override
	//public String toString() {
	//	return "Matrix [M=" + M + ", N=" + N + ", data=" + Arrays.toString(data) + "]";
	//}
    
}
