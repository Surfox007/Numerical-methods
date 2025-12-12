#include <bits/stdc++.h>
#define int long long
using namespace std;

double u_cal(double u, int n)
{
    double temp = u;
    for (int i = 1; i < n; i++)
    {
        temp *= (u - i);
    }
    return temp;
}

int fact(int n)
{
    if (n == 1)
        return 1;
    return n * fact(n - 1);
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
    if (!(inputFile >> n))
    {
        cerr << "Error reading number of data points" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }
    // cin>>n;
    vector<vector<double>> y(n, vector<double>(n, 0));
    vector<double> x(n, 0);
    for (int i = 0; i < n; i++)
    {
        // cin>>x[i]>> y[i][0];
        if (!(inputFile >> x[i] >> y[i][0]))
        {
            cerr << "Error reading data point " << i + 1 << endl;
            inputFile.close();
            outputFile.close();
            return 1;
        }
    }

    // forward difference
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < n - i; j++)
        {
            y[j][i] = y[j + 1][i - 1] - y[j][i - 1];
        }
    }

    outputFile << fixed << setprecision(4);

    // forward difference table
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

    double a;
    // cout<<"Interpolate at: ";
    // cin>> a;
    if (!(inputFile >> a))
    {
        cerr << "Error reading interpolation point" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }

    double sum = y[0][0];
    double u = (a - x[0]) / (x[1] - x[0]);
    for (int i = 1; i < n; i++)
    {
        sum += u_cal(u, i) * (1 / fact(i)) * y[0][1];
    }

    // cout<<"Value at "<<a<<": "<<sum<<endl;
    if (!(outputFile << "Value at " << a << " = " << sum << endl))
    {
        cerr << "Error writing output" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }

    return 0;
}
/* Testcase
4
45 0.7071
50 0.766
55 0.8192
60 0.866
52
*/
