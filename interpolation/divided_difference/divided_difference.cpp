#include <bits/stdc++.h>
#define int long long
using namespace std;

double proterm(vector<double> &x, int n, double a)
{
    double temp = 1;
    for (int i = 0; i < n; i++)
    {
        temp *= (a - x[i]);
    }
    return temp;
}

int32_t main()
{
    ifstream inputFile("input.txt");
    
    if (!inputFile.is_open())
    {
        cerr << "Error while opening input file" << endl;
        return 1;
    }
    
    int n;
    // cin>>n;
    if (!(inputFile >> n))
    {
        cerr << "Error reading number of data points" << endl;
        inputFile.close();
        return 1;
    }
    vector<vector<double>> y(n, vector<double>(n, 0));
    vector<double> x(n, 0);
    for (int i = 0; i < n; i++)
    {
        // cin>>x[i]>> y[i][0];
        if (!(inputFile >> x[i] >> y[i][0]))
        {
            cerr << "Error reading data point " << i + 1 << endl;
            inputFile.close();
        
            return 1;
        }
    }

    // divided difference table
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < n - i; j++)
        {
            y[j][i] = (y[j][i - 1] - y[j + 1][i - 1]) / (x[j] - x[j + i]);
        }
    }

    //cout << fixed << setprecision(3);
    ofstream outputFile ("output.txt", ios::trunc);
    if (!outputFile.is_open())
    {
        cerr << "Error while opening output file" << endl;
        inputFile.close();
        return 1;
    }
    outputFile << fixed << setprecision(6);
    // printing divided difference table
    for (int i = 0; i < n; i++)
    {
        //cout << setw(7) << x[i] << "\t";
        if (!(outputFile << setw(7) << x[i] << "\t"))
        {
            cerr << "Error writing to output file" << endl;
            inputFile.close();
            outputFile.close();
            return 1;
        }
        for (int j = 0; j < n; j++)
        {
            //cout << setw(7) << y[i][j] << "\t";
            if (!(outputFile << setw(7) << y[i][j] << "\t"))
            {
                cerr << "Error writing to output file" << endl;
                inputFile.close();
                outputFile.close();
                return 1;
            }
        }
        //cout << endl;
        if (!(outputFile << endl))
        {
            cerr << "Error writing to output file" << endl;
            inputFile.close();
            outputFile.close();
            return 1;
        }
    }

    double a;
    //cout << "Value  at: ";

    //cin >> a;
    if (!(inputFile >> a))
    {
        cerr << "Error reading interpolation point" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }

    

    double sum = y[0][0];
    for (int i = 1; i < n; i++)
    {
        sum += proterm(x, i, a) * y[0][i];
    }

    //cout << "Value at " << a << ": " << sum << endl;
    outputFile << "Value at " << a << " : "<< sum << endl;
    return 0;
}
/* Testcase
4
300	2.4771
304	2.4829
305	2.4843
307	2.4871
301
*/

/*
Testcase 2
4
5 12
6 13
9 14
11 16
7
*/
