#include <bits/stdc++.h>
#define int long long
const double EPSILON = 1e-9;
using namespace std;

void get_submatrix(vector<vector<double>> mat, vector<vector<double>> &temp, int p, int q, int n)
{
    int r = 0;
    for (int i = 0; i < n; i++){
        if (i != p){
            int s = 0;
            for (int j = 0; j < n; j++){
                if (j != q){
                    temp[r][s++] = mat[i][j];
                }
            }

            r++;
        }
    }
}

double get_determinant(vector<vector<double>> mat, int n)
{
    if (n == 1){
        return mat[0][0];
    }

    if (n == 2){
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }

    int D = 0;
    int sign = 1;
    vector<vector<double>> temp(n - 1, vector<double>(n - 1, 0));
    for (int c = 0; c < n; c++){
        get_submatrix(mat, temp, 0, c, n);

        D += sign * mat[0][c] * get_determinant(temp, n - 1);
        sign *= -1;
    }

    return D;
}

bool get_inverse(vector<vector<double>> mat, vector<vector<double>> &inverse, int n)
{
    double A = get_determinant(mat, n);
    if (abs(A) < EPSILON){
        return false;
    }

    if (n == 1){
        inverse[0][0] = 1.0 / A;
        return true;
    }
    double inv_det = 1.0 / A;
    vector<vector<double>> temp(n - 1, vector<double>(n - 1, 0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            get_submatrix(mat, temp, i, j, n);

            double minor = get_determinant(temp, n - 1);

            double sign = ((i + j) % 2 == 0 ? 1.0 : -1.0);

            double cofactor = sign * minor;

            inverse[j][i] = cofactor * inv_det;
        }
    }

    return true;
}

void print_matrix(vector<vector<double>> &m, int n){
    cout << fixed << setprecision(3);
    for (int i = 0; i < n; i++){
        cout << "|";
        for (int j = 0; j < n; j++)
        {
            cout << setw(4) << m[i][j] << " ";
        }
        cout << "|\n";
    }
}

int32_t main(){

    ifstream inputFile("input.txt");
    if (!inputFile.is_open()){
        cerr << "Error while opening input file" << endl;
    }

    ofstream outputFile("output.txt", ios::trunc);
    if (!inputFile.is_open()){
        cerr << "Error while opening input file" << endl;
    }

    int n;
    // cin>>n;
    if (!(inputFile >> n)){
        cerr << "Error in reading matrix dimension N" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }

    vector<vector<double>> matrix(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (!(inputFile >> matrix[i][j])){
                cerr << "Error in reading matrix element at (" << i << "," << j << ")" << endl;
                inputFile.close();
                outputFile.close();
                return 1;
            }
        }
    }
    vector<vector<double>> inverse(n, vector<double>(n, 0));
    get_inverse(matrix, inverse, n);

    outputFile << fixed << setprecision(3);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            outputFile << setw(6) << inverse[i][j] << "\t";
        }
        outputFile << "\n";
    }
    return 0;
}
