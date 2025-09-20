import java.math.BigInteger;

public class PolynomialSolver {

    public static BigInteger baseToDecimal(String value, int base) {
        String digits = "0123456789abcdef";
        BigInteger result = BigInteger.ZERO;
        BigInteger b = BigInteger.valueOf(base);

        value = value.toLowerCase();
        for (int i = 0; i < value.length(); i++) {
            int digitValue = digits.indexOf(value.charAt(i));
            if (digitValue < 0 || digitValue >= base) {
                throw new IllegalArgumentException("Invalid digit " + value.charAt(i) + " for base " + base);
            }
            result = result.multiply(b).add(BigInteger.valueOf(digitValue));
        }
        return result;
    }

    public static double[] powers(double x, int m) {
        double[] pows = new double[m + 1];
        pows[0] = 1.0;
        for (int i = 1; i <= m; i++) {
            pows[i] = pows[i - 1] * x;
        }
        return pows;
    }

   
    public static double[] solveSystem(double[][] V, int k) {
        int n = V.length;
        int m = k - 1;

        double[][] A = new double[m][m];
        double[] b = new double[m];

        for (int i = 0; i < m; i++) {
           
            for (int j = 0; j < m; j++) {
                A[i][j] = V[i][j];
            }
            
            b[i] = -V[i][m];
        }

        
        return gaussianElimination(A, b);
    }

   
    public static double[] gaussianElimination(double[][] A, double[] b) {
        int n = b.length;

        for (int p = 0; p < n; p++) {
            
            int max = p;
            for (int i = p + 1; i < n; i++) {
                if (Math.abs(A[i][p]) > Math.abs(A[max][p])) {
                    max = i;
                }
            }

           
            double[] tempRow = A[p];
            A[p] = A[max];
            A[max] = tempRow;

            double tempVal = b[p];
            b[p] = b[max];
            b[max] = tempVal;

           
            if (Math.abs(A[p][p]) <= 1e-12) {
                throw new ArithmeticException("Matrix is singular or nearly singular");
            }

            for (int i = p + 1; i < n; i++) {
                double alpha = A[i][p] / A[p][p];
                b[i] -= alpha * b[p];
                for (int j = p; j < n; j++) {
                    A[i][j] -= alpha * A[p][j];
                }
            }
        }

        
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += A[i][j] * x[j];
            }
            x[i] = (b[i] - sum) / A[i][i];
        }

        return x;
    }

    public static void main(String[] args) {
        // Given roots data (base, value)
        String[][] rootsData = {
                {"6", "13444211440455345511"},
                {"15", "aed7015a346d635"},
                {"15", "6aeeb69631c227c"},
                {"16", "e1b5e05623d881f"},
                {"8", "316034514573652620673"},
                {"3", "2122212201122002221120200210011020220200"},
                {"3", "20120221122211000100210021102001201112121"},
                {"6", "20220554335330240002224253"},
                {"12", "45153788322a1255483"},
                {"7", "1101613130313526312514143"}
        };

        int n = 10;  
        int k = 7;   

        double[] roots = new double[n];
        for (int i = 0; i < n; i++) {
            int base = Integer.parseInt(rootsData[i][0]);
            BigInteger val = baseToDecimal(rootsData[i][1], base);
            roots[i] = val.doubleValue();  
        }

     
        int m = k - 1;
        double[][] V = new double[m][k];
        for (int i = 0; i < m; i++) {
            double[] pows = powers(roots[i], m);
            System.arraycopy(pows, 0, V[i], 0, k);
        }

        
        double[] coeffsWithoutLast = solveSystem(V, k);

        
        System.out.println("Polynomial coefficients:");
        for (int i = 0; i < m; i++) {
            System.out.printf("a%d = %.6e\n", i, coeffsWithoutLast[i]);
        }
        System.out.printf("a%d = 1.0 (assumed)\n", m);
    }



Ouput:
    TestCase 1:
    Polynomial coefficients:
a0 = 28.000000
a1 = -11.000000
a2 = 1.0 (assumed)

    Testcase 2:
Polynomial coefficients:
a0 = 2.555330e+105
a1 = -2.592353e+90
a2 = 2.459174e+73
a3 = -7.010125e+55
a4 = 6.289601e+37
a5 = -1.603559e+19
a6 = 1.0 (assumed)
