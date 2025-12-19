#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

ld polyValue(const vector<ld>& coeff, ld x)
{
    ld val = 0, p = 1;
    for (ld c : coeff)
    {
        val += c * p;
        p *= x;
    }
    return val;
}

ld trueDerivative(const vector<ld>& coeff, ld x)
{
    ld val = 0;
    for (ll i = 1; i < coeff.size(); i++)
        val += i * coeff[i] * pow(x, i - 1);
    return val;
}

vector<vector<ld>> buildBackwardDiffTable(const vector<ld>& y)
{
    ll n = y.size();
    vector<vector<ld>> diff(n, vector<ld>(n, 0));
    for (ll i = 0; i < n; i++) diff[i][0] = y[i];

    for (ll j = 1; j < n; j++)
        for (ll i = n - 1; i >= j; i--)
            diff[i][j] = diff[i][j - 1] - diff[i - 1][j - 1];

    return diff;
}

ld newtonBackwardDerivative5(const vector<vector<ld>>& diff, ld xn, ld h, ld xp, ll n)
{
    ld s = (xp - xn) / h;
    ld der = 0;

    if (n >= 2) der += diff[n - 1][1];
    if (n >= 3) der += ((2 * s + 1) / 2.0) * diff[n - 1][2];
    if (n >= 4) der += ((3 * s * s + 6 * s + 2) / 6.0) * diff[n - 1][3];
    if (n >= 5) der += ((4 * s * s * s + 18 * s * s + 22 * s + 6) / 24.0) * diff[n - 1][4];
    if (n >= 6) der += ((5 * s * s * s * s + 40 * s * s * s + 105 * s * s + 100 * s + 24) / 120.0) * diff[n - 1][5];

    return der / h;
}

string polyString(const vector<ld>& c)
{
    stringstream ss;
    bool first = true;

    for (ll i = 0; i < c.size(); i++)
    {
        if (fabs(c[i]) < 1e-12) continue;
        if (!first) ss << (c[i] >= 0 ? " + " : " - ");
        else if (c[i] < 0) ss << "-";

        ss << fixed << setprecision(4) << fabs(c[i]);
        if (i == 1) ss << "x";
        else if (i > 1) ss << "x^" << i;

        first = false;
    }

    return ss.str();
}

int main()
{
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    ll n, deg;
    ll tc = 1;

    // EOF loop
    while (cin >> n >> deg)
    {
        vector<ld> coeff(deg + 1);
        for (ll i = 0; i <= deg; i++) cin >> coeff[i];

        vector<ld> x(n), y(n);
        for (ll i = 0; i < n; i++) cin >> x[i];

        for (ll i = 0; i < n; i++) y[i] = polyValue(coeff, x[i]);

        ld xp;
        cin >> xp;

        ld h = x[1] - x[0];
        auto diff = buildBackwardDiffTable(y);

        ld approx = newtonBackwardDerivative5(diff, x[n - 1], h, xp, n);
        ld trueVal = trueDerivative(coeff, xp);
        ld err = fabs((trueVal - approx) / trueVal) * 100.0;

        string poly = polyString(coeff);

        cout << "TEST CASE #" << tc++ << "\n";
        cout << "Polynomial: f(x) = " << poly << "\n";
        cout << "Number of points: " << n << "\n";

        cout << "x-values: ";
        for (ld v : x) cout << v << " ";
        cout << "\n";

        cout << "y-values: ";
        for (ld v : y) cout << fixed << setprecision(6) << v << " ";
        cout << "\n";

        cout << fixed << setprecision(6);
        cout << "Step size (h): " << h << "\n";
        cout << "Differentiation point: " << xp << "\n";

        cout << "Backward Difference Table:\n";
        for (ll i = 0; i < n; i++)
        {
            cout << "Row " << i << ": ";
            for (ll j = 0; j <= i; j++) cout << setw(12) << diff[i][j] << " ";
            cout << "\n";
        }

        cout << "Approximate derivative: " << approx << "\n";
        cout << "True derivative: " << trueVal << "\n";
        cout << "Percentage error: " << err << " %\n";
        cout << "---------------\n";
    }

    return 0;
}
