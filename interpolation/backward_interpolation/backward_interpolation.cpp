#include <bits/stdc++.h>
#define int long long
using namespace std;

double u_cal(double u, int n)
{
    double temp = u;
    for (int i = 1; i < n; i++)
    {
        temp *= (u + i);
    }
    return temp;
}

double fact(int n)
{
    if (n == 1)
        return 1;
    return n * fact(n - 1);
}

double backward_interpolation(vector<double> &x, vector<vector<double>> &y, int n, int a)
{
    for (int i = 1; i < n; i++)
    {
        for (int j = n - 1; j <= i; j++)
        {
            y[j][i] = y[j][i - 1] - y[j - 1][i - 1];
        }
    }

    double u = (a - x[n - 1]) / (x[1] - x[0]);
    double sum = y[n - 1][0];
    for (int i = 1; i < n; i++)
    {
        sum += (u_cal(u, i) * y[n - 1][i]) / fact(i);
    }
    return sum;
}

int32_t main()
{
    ifstream inputFile("input.txt");
    ofstream outputFile("output.txt", ios::trunc);
    if (!inputFile.is_open())
    {
        cerr << "Error while opening input file" << endl;
        return 1;
    }

    if (!outputFile.is_open())
    {
        cerr << "Error while opening output file" << endl;
        inputFile.close();
        return 1;
    }

    int n;
    // cin>>n;
    if (!(inputFile >> n))
    {
        cerr << "Error reading number of data points" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }
    vector<double> x(n, 0);
    vector<vector<double>> y(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++)
    {
        if (!(inputFile >> x[i] >> y[i][0]))
        {
            cerr << "Error reading data point " << i + 1 << endl;
            inputFile.close();
            outputFile.close();
            return 1;
        }
        // cin >> x[i] >> y[i][0];
    }
    double a;
    // cout << "Interpole at: ";
    //  cin >> a;
    if (!(inputFile >> a))
    {
        cerr << "Error reading interpolation point" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }
    double ans = backward_interpolation(x, y, n, a);

    outputFile << fixed << setprecision(4);
    // printing backward difference table
    for (int i = 0; i < n; i++)
    {
        // cout<<setw(4)<<x[i]<<"\t";
        if (!(outputFile << setw(4) << x[i] << "\t"))
        {
            cerr << "Error writing output" << endl;
            inputFile.close();
            outputFile.close();
            return 1;
        }
        for (int j = 0; j < n; j++)
        {
            // cout<<setw(4)<<y[i][j]<<"\t";
            if (!(outputFile << setw(4) << y[i][j] << "\t"))
            {
                cerr << "Error writing output" << endl;
                inputFile.close();
                outputFile.close();
                return 1;
            }
        }
        outputFile << "\n";
    }
    // cout << ans << endl;
    if (!(outputFile << "Interpole value at " << a << " = " << ans << endl))
    {
        cerr << "Error writing output" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }
    return 0;
}
