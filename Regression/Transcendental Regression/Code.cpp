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

        vector<ld> lnY(n);
        bool valid=true;
        for(ll i=0;i<n;i++)
        {
            if(y[i]<=0) { valid=false; break; }
            lnY[i] = log(y[i]);
        }

        if(valid)
        {
            ld sumX=0, sumY=0, sumXY=0, sumX2=0;
            for(ll i=0;i<n;i++)
            {
                sumX += x[i];
                sumY += lnY[i];
                sumXY += x[i]*lnY[i];
                sumX2 += x[i]*x[i];
            }

            ld b = (n*sumXY - sumX*sumY)/(n*sumX2 - sumX*sumX);
            ld ln_a = (sumY - b*sumX)/n;
            ld a = exp(ln_a);

            cout << fixed << setprecision(10);
            cout << "Exponential Regression: y = " << a << " * e^(" << b << "*x)\n";
        }
        else
        {
            cout << "Exponential Regression: Not valid (y <= 0)\n";
        }

        cout << "---------------\n";
    }
    return 0;
}
