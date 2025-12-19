#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

const ld eps = 1e-9;
// f(x)
ld f(ld x, vector<pair<ld, ll>> &poly)
{
    ld sum = 0;
    for (auto [coeff, power] : poly)
    {
        sum += coeff * pow(x, power);
    }
    return sum;
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    ll n;
    while (cin >> n)
    {
        vector<pair<ld, ll>> poly(n + 1);
        ld an = 0;
        for (ll i = 0; i <= n; i++)
        {
            ld c;
            ll p;
            cin >> c >> p;
            poly[i] = {c, p};
            if (p == n)
                an = c;
        }

        ld xmax = 0;
        for (auto [c, p] : poly)
        {
            if (p != n)
                xmax = max(xmax, fabsl(c / an));
        }
        xmax += 1;
        ld step;
        cin >> step;
        vector<pair<ld, ld>> interval;
        ld l = -xmax;
        while (l <= xmax)
        {
            while (l + step <= xmax && f(l, poly) * f(l + step, poly) > 0.0)
            {
                l += step;
            }
            if (l + step <= xmax)
                interval.push_back({l, l + step});
            l += step;
        }
        for (auto [u, v] : interval)
        {
            ld x0 = u;
            ld x1 = v;
            ll iter = 0;

            while (1)
            {
                iter++;
                ld f0 = f(x0, poly);
                ld f1 = f(x1, poly);

                if (fabsl(f1 - f0) < eps) // avoid division by zero
                    break;

                ld x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

                if (fabsl(x2 - x1) < eps)
                {
                    cout << fixed << setprecision(10);
                    cout << "Root = " << x2 << "  Iterations = " << iter << "\n";
                    break;
                }

                x0 = x1;
                x1 = x2;
            }
        }
        if (interval.size() != n)
        {
            cout << "All other roots are imaginary\n";
        }
        cout << "---------------\n";
    }
    return 0;
}