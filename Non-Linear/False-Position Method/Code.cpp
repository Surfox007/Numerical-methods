#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

const ld eps = 1e-9;

ld f(ld x, const vector<pair<ld, ll>> &poly)
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
            ll iter = 0;
            while (u <= v)
            {
                iter++;
                ld fu = f(u, poly);
                ld fv = f(v, poly);

                ld mid = (u * fv - v * fu) / (fv - fu);
                ld fmid = f(mid, poly);

                if (fabsl(fmid) < eps)
                {
                    cout << fixed << setprecision(10);
                    cout << "Root = " << mid << "  Iterations = " << iter << "\n";
                    break;
                }

                if (f(u, poly) * fmid < 0)
                    v = mid;
                else
                    u = mid;
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