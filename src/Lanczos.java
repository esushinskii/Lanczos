import org.apache.commons.math3.linear.*;


public class Lanczos {
    public static void main(String[] args) {
        double[][] a = new double[][] {{1, 2, 3}, {1, 4, 7}, {2, 1, 3}};
        int m = 3;

        double[][] t = getTridiagonalMatrix(a, m);
        printMatrix(t);
    }


    public static double[][] getTridiagonalMatrix(double[][] a, int m) {
        double[] alpha = new double[m + 1];
        double[] betta = new double[m + 1];

        int n = a.length;

        RealVector v0 = new ArrayRealVector(n);
        RealVector v1 = new ArrayRealVector(n, 1.0);

        RealVector wx = new ArrayRealVector(n);
        RealVector w = new ArrayRealVector(n);

        RealMatrix t = new Array2DRowRealMatrix(a);

        for (int i = 1; i < m; i++) {
            wx = t.operate(v1);
            alpha[i] = wx.dotProduct(v1);
            w = wx.subtract(v1.mapMultiply(alpha[i])).subtract(v0.mapMultiply(betta[i]));
            betta[i + 1] = w.getNorm();
            v0 = v1.copy();
            v1 = w.mapMultiply(1 / betta[i + 1]);
        }

        return makeTridiagonalMatrix(alpha, betta);
    }


    public static double[][] makeTridiagonalMatrix(double[] alpha, double[] betta) {
        int m = alpha.length - 1;
        double[][] T = new double[m][m];

        T[0][0] = alpha[1];
        for (int i = 1; i < alpha.length - 1; i++) {
            T[i][i] = alpha[i + 1];
            T[i - 1][i] = betta[i + 1];
            T[i][i - 1] = betta[i + 1];
        }

        return T;
    }


    public static void printMatrix(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.printf("%8.2f ", matrix[i][j]);
            }
            System.out.println();
        }
    }
}
