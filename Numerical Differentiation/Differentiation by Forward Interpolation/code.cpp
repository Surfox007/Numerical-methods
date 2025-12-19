#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

ll fact(ll n)
{
    if (n <= 1)
        return 1;
    return n * fact(n - 1);
}

vector<ld> coef;

// f(x)
ld f(ld x)
{
    ld sum = 0;
    for (ll i = 0; i < coef.size(); i++)
        sum += coef[i] * pow(x, i);
    return sum;
}

// f'(x)
ld f1(ld x)
{
    ld sum = 0;
    for (ll i = 1; i < coef.size(); i++)
        sum += coef[i] * i * pow(x, i - 1);
    return sum;
}

// f''(x)
ld f2(ld x)
{
    ld sum = 0;
    for (ll i = 2; i < coef.size(); i++)
        sum += coef[i] * i * (i - 1) * pow(x, i - 2);
    return sum;
}

// Difference table
vector<vector<ld>> diff_table(vector<ld> &y)
{
    ll n = y.size();
    vector<vector<ld>> diff(n, vector<ld>(n, 0));

    for (ll i = 0; i < n; i++)
        diff[i][0] = y[i];

    for (ll j = 1; j < n; j++)
        for (ll i = 0; i < n - j; i++)
            diff[i][j] = diff[i + 1][j - 1] - diff[i][j - 1];

    return diff;
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
    ll tc = 1;
    ll n, deg;

    while (cin >> n >> deg)
    {
        coef.assign(deg + 1, 0);
        for (ll i = 0; i <= deg; i++)
            cin >> coef[i];

        ld a, b, X;
        cin >> a >> b >> X;

        ld h = (b - a) / n;

        vector<ld> x(n), y(n);
        for (ll i = 0; i < n; i++)
        {
            x[i] = a + i * h;
            y[i] = f(x[i]);
        }

        vector<vector<ld>> diff = diff_table(y);
        ld u = (X - x[0]) / h;

        ld y1 = (diff[0][1] + (2 * u - 1) * diff[0][2] / fact(2) + (3 * u * u - 6 * u + 2) * diff[0][3] / fact(3)) / h;

        ld y2 = (diff[0][2] + (u - 1) * diff[0][3]) / (h * h);

        ld exact_y1 = f1(X);
        ld exact_y2 = f2(X);

        ld err1 = fabsl((exact_y1 - y1) / exact_y1) * 100;
        ld err2 = fabsl((exact_y2 - y2) / exact_y2) * 100;

        cout << "\nTEST CASE #" << tc++ << "\n";
        cout << "n = " << n << ", degree = " << deg
             << ", a = " << a << ", b = " << b << ", X = " << X << "\n";

        cout << "Difference Table:\n";
        for (ll i = 0; i < n; i++)
        {
            for (ll j = 0; j < n; j++)
                cout << setw(12) << diff[i][j] << " ";
            cout << "\n";
        }

        cout << fixed << setprecision(10);
        cout << "Numerical y'  = " << y1 << "\n";
        cout << "Numerical y'' = " << y2 << "\n";
        cout << "Exact y'      = " << exact_y1 << "\n";
        cout << "Exact y''     = " << exact_y2 << "\n";
        cout << "First derivative error  = " << err1 << "%\n";
        cout << "Second derivative error = " << err2 << "%\n";
        cout << "---------------\n";
    }

    return 0;
}
