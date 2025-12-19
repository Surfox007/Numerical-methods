#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    ll n;
    while(cin >> n)
    {
        vector<ld> x(n), y(n);
        for(ll i=0;i<n;i++)
            cin >> x[i] >> y[i];

        ld sumX=0, sumY=0, sumXY=0, sumX2=0;
        for(ll i=0;i<n;i++)
        {
            sumX += x[i];
            sumY += y[i];
            sumXY += x[i]*y[i];
            sumX2 += x[i]*x[i];
        }

        ld a = (n*sumXY - sumX*sumY)/(n*sumX2 - sumX*sumX);
        ld b = (sumY - a*sumX)/n;

        cout << fixed << setprecision(10);
        cout << "Linear Regression: y = " << a << " * x + " << b << "\n";
        cout << "---------------\n";
    }
    return 0;
}
