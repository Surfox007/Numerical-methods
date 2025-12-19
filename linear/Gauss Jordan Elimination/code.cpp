#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

int matrixRank(vector<vector<double>> M) {
    int n = M.size();
    int m = M[0].size();
    int rank = 0;

    for (int col = 0, row = 0; col < m && row < n; col++) {
        int sel = row;
        for (int i = row; i < n; i++) {
            if (fabs(M[i][col]) > fabs(M[sel][col])) sel = i;
        }
        if (fabs(M[sel][col]) < EPS) continue;

        swap(M[sel], M[row]);
        double div = M[row][col];
        for (int j = col; j < m; j++) M[row][j] /= div;

        for (int i = 0; i < n; i++) {
            if (i != row && fabs(M[i][col]) > EPS) {
                double factor = M[i][col];
                for (int j = col; j < m; j++) {
                    M[i][j] -= factor * M[row][j];
                }
            }
        }
        row++;
        rank++;
    }
    return rank;
}

int checkSystemType(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> Aug(n, vector<double>(n + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) Aug[i][j] = A[i][j];
        Aug[i][n] = b[i];
    }

    int rankA = matrixRank(A);
    int rankAug = matrixRank(Aug);

    if (rankA == rankAug) {
        if (rankA == n) return 0;
        else return 1;        
    } else {
        return 2;              
    }
}

vector<double> gaussJordan(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        if (fabs(A[i][i]) < EPS) {
            throw runtime_error("Singular matrix or no unique solution.");
        }

        double pivot = A[i][i];
        for (int j = 0; j < n; j++) {
            A[i][j] /= pivot;
        }
        b[i] /= pivot;

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; j++) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    }
    return b;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T; 
    fin >> T;

    for (int t = 1; t <= T; t++) {
        int n; 
        fin >> n;

        vector<vector<double>> A(n, vector<double>(n));
        vector<double> b(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fin >> A[i][j];
            }
        }
        for (int i = 0; i < n; i++) {
            fin >> b[i];
        }

        fout << "Case #" << t << ":\n";

        int systemType = checkSystemType(A, b);
        if (systemType == 2) {
            fout << "System Type: NO SOLUTION\n\n";
            continue;
        } else if (systemType == 1) {
            fout << "System Type: INFINITE SOLUTIONS\n\n";
            continue;
        } else {
            fout << "System Type: UNIQUE SOLUTION\n";
            try {
                vector<double> solution = gaussJordan(A, b);
                for (int i = 0; i < n; i++) {
                    fout << "x" << i + 1 << " = " 
                         << fixed << setprecision(6) << solution[i] << "\n";
                }
            } catch (runtime_error &e) {
                fout << "Error: " << e.what() << "\n";
            }
            fout << "\n";
        }
    }

    fin.close();
    fout.close();
    return 0;
}
