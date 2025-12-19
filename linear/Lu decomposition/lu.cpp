#include <bits/stdc++.h>
using namespace std;

const double EPSILON = 1e-9;

bool isZero(double val) {
    return fabs(val) < EPSILON;
}

bool luDecompose(vector<vector<double>>& A, vector<int>& P, int n) {
    for (int i = 0; i < n; i++) {
        P[i] = i;
    }
    
    for (int k = 0; k < n; k++) {
        double maxVal = 0;
        int pivotRow = k;
        for (int i = k; i < n; i++) {
            if (fabs(A[i][k]) > maxVal) {
                maxVal = fabs(A[i][k]);
                pivotRow = i;
            }
        }

        if (pivotRow != k) {
            swap(A[k], A[pivotRow]);
            swap(P[k], P[pivotRow]);
        }
        
        if (isZero(A[k][k])) {
            return false;
        }
        
        for (int i = k + 1; i < n; i++) {
            A[i][k] = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }
    return true;
}


void forwardSubstitution(const vector<vector<double>>& LU, const vector<int>& P, 
                         const vector<double>& b, vector<double>& y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = b[P[i]];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j] * y[j];
        }
    }
}

void backwardSubstitution(const vector<vector<double>>& LU, 
                          const vector<double>& y, vector<double>& x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= LU[i][j] * x[j];
        }
        x[i] = x[i] / LU[i][i];
    }
}

int checkSystemType(vector<vector<double>> A, vector<double> b, int n) {
    vector<int> P(n);
    vector<vector<double>> augmented = A;
    
    for (int i = 0; i < n; i++) {
        augmented[i].push_back(b[i]);
    }
    
    int rank_A = 0, rank_Ab = 0;
    
    for (int col = 0; col < n; col++) {
        int pivotRow = -1;
        for (int row = rank_A; row < n; row++) {
            if (! isZero(A[row][col])) {
                pivotRow = row;
                break;
            }
        }
        
        if (pivotRow == -1) continue;
        
        if (pivotRow != rank_A) {
            swap(A[rank_A], A[pivotRow]);
            swap(augmented[rank_A], augmented[pivotRow]);
        }

        for (int row = rank_A + 1; row < n; row++) {
            if (! isZero(A[row][col])) {
                double factor = A[row][col] / A[rank_A][col];
                for (int j = 0; j < n; j++) {
                    A[row][j] -= factor * A[rank_A][j];
                }
                for (int j = 0; j <= n; j++) {
                    augmented[row][j] -= factor * augmented[rank_A][j];
                }
            }
        }
        rank_A++;
    }
    
    rank_Ab = rank_A;
    for (int i = rank_A; i < n; i++) {
        if (!isZero(augmented[i][n])) {
            return 2; 
        }
    }
    
    if (rank_A < n) {
        return 1;
    }
    
    return 0;
}

void solveSystem(ifstream& input, ofstream& output, int testCase) {
    int n;
    input >> n;
    
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n), x(n), y(n);
    vector<int> P(n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            input >> A[i][j];
        }
    }
    
    for (int i = 0; i < n; i++) {
        input >> b[i];
    }
    
    output << "========== TEST CASE " << testCase << " ==========\n";
    output << "System size: " << n << "x" << n << "\n\n";
    
    int systemType = checkSystemType(A, b, n);
    
    if (systemType == 2) {
        output << "System Type: NO SOLUTION\n";
        output << "The system is inconsistent.\n\n";
        return;
    } else if (systemType == 1) {
        output << "System Type: INFINITE SOLUTIONS\n";
        output << "The system has infinitely many solutions.\n\n";
        return;
    }
    
    output << "System Type:  UNIQUE SOLUTION\n\n";
    vector<vector<double>> LU = A;
    if (!luDecompose(LU, P, n)) {
        output << "Error: Matrix is singular\n\n";
        return;
    }
    
    forwardSubstitution(LU, P, b, y, n);
    backwardSubstitution(LU, y, x, n);
    
    output << "Solution:\n";
    for (int i = 0; i < n; i++) {
        output << "x[" << i << "] = " << fixed << setprecision(6) << x[i] << "\n";
    }
    output << "\n";
}

int main() {
    ifstream input("input.txt");
    ofstream output("output.txt");
    
    if (!input.is_open() || !output.is_open()) {
        cerr << "Error opening files!\n";
        return 1;
    }
    
    int numTests;
    input >> numTests;
    
    output << "LU DECOMPOSITION - LINEAR EQUATION SOLVER\n";
    output << "==========================================\n\n";
    
    for (int i = 1; i <= numTests; i++) {
        solveSystem(input, output, i);
    }
    
    input.close();
    output.close();
    
    cout << "LU decomposition completed.  Check lu_output.txt\n";
    return 0;
}