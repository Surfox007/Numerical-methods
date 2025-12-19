#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

vector<ld> gauss(vector<vector<ld>> A, vector<ld> B)
{
    ll n = B.size();
    for(ll i=0;i<n;i++)
    {
        ll pivot=i;
        for(ll j=i+1;j<n;j++)
            if(fabsl(A[j][i])>fabsl(A[pivot][i])) pivot=j;
        swap(A[i], A[pivot]);
        swap(B[i], B[pivot]);

        for(ll j=i+1;j<n;j++)
        {
            ld factor = A[j][i]/A[i][i];
            for(ll k=i;k<n;k++) A[j][k]-=factor*A[i][k];
            B[j]-=factor*B[i];
        }
    }

    vector<ld> X(n);
    for(ll i=n-1;i>=0;i--)
    {
        X[i]=B[i];
        for(ll j=i+1;j<n;j++) X[i]-=A[i][j]*X[j];
        X[i]/=A[i][i];
    }
    return X;
}

int main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    ll n;
    while(cin >> n)
    {
        ll degree;
        cin >> degree;
        vector<ld> x(n), y(n);
        for(ll i=0;i<n;i++)
            cin >> x[i] >> y[i];

        vector<vector<ld>> A(degree+1, vector<ld>(degree+1,0));
        vector<ld> B(degree+1,0);

        for(ll i=0;i<=degree;i++)
        {
            for(ll j=0;j<=degree;j++)
                for(ll k=0;k<n;k++)
                    A[i][j] += pow(x[k], i+j);

            for(ll k=0;k<n;k++)
                B[i] += y[k]*pow(x[k], i);
        }

        vector<ld> coeff = gauss(A,B);

        cout << fixed << setprecision(10);
        cout << "Polynomial Regression: y = ";
        for(ll i=0;i<=degree;i++)
        {
            if(i>0) cout << " + ";
            cout << coeff[i] << "*x^" << i;
        }
        cout << "\n";
        cout << "---------------\n";
    }
    return 0;
}
