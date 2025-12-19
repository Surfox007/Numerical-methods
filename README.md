# Numerical Methods Implementation

This repository contains implementations of various numerical methods for solving mathematical problems, including linear and non-linear equations, interpolation, regression, integration, and ordinary differential equations.

# Table of Contents

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)
  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  - [Newton Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)

- [Solution of Linear Equations](#solution-of-linear-equations)
  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion Method](#matrix-inversion-method)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Interpolation](#interpolation)
  - [Newton's Forward Interpolation](#newtons-forward-interpolation)
    - [Theory](#forward-interpolation-theory)
    - [Code](#forward-interpolation-code)
    - [Input](#forward-interpolation-input)
    - [Output](#forward-interpolation-output)
  - [Newton's Backward Interpolation](#newtons-backward-interpolation)
    - [Theory](#backward-interpolation-theory)
    - [Code](#backward-interpolation-code)
    - [Input](#backward-interpolation-input)
    - [Output](#backward-interpolation-output)
  - [Divided Difference Interpolation](#divided-difference-interpolation)
    - [Theory](#divided-difference-theory)
    - [Code](#divided-difference-code)
    - [Input](#divided-difference-input)
    - [Output](#divided-difference-output)
- [Solution of Numerical Differentiation](#solution-of-numerical-differentiation)
  - [Differentiation by Forward Interpolation Method](#differentiation-by-forward-interpolation-method)
    - [Theory](#differentiation-by-forward-interpolation-theory)
    - [Code](#differentiation-by-forward-interpolation-code)
    - [Input](#differentiation-by-forward-interpolation-input)
    - [Output](#differentiation-by-forward-interpolation-output)
  - [Differentiation by Backward Interpolation Method](#differentiation-by-backward-interpolation-method)
    - [Theory](#differentiation-by-backward-interpolation-theory)
    - [Code](#differentiation-by-backward-interpolation-position-code)
    - [Input](#differentiation-by-backward-interpolation-input)
    - [Output](#differentiation-by-backward-interpolation-output)
- [Regression (Curve Fitting)](#regression-curve-fitting)
  - [Linear Regression](#linear-regression)
    - [Theory](#linear-regression-theory)
    - [Code](#linear-regression-code)
    - [Input](#linear-regression-input)
    - [Output](#linear-regression-output)
  - [Transcendental Equation Regression](#transcendental-equation-regression)
    - [Theory](#transcendental-regression-theory)
    - [Code](#transcendental-regression-code)
    - [Input](#transcendental-regression-input)
    - [Output](#transcendental-regression-output)
  - [Polynomial Regression](#polynomial-regression)
    - [Theory](#polynomial-regression-theory)
    - [Code](#polynomial-regression-code)
    - [Input](#polynomial-regression-input)
    - [Output](#polynomial-regression-output)

- [Numerical Integration](#numerical-integration)
  - [Simpson's 1/3rd Rule](#simpsons-13rd-rule)
    - [Theory](#simpsons-13rd-theory)
    - [Code](#simpsons-13rd-code)
    - [Input](#simpsons-13rd-input)
    - [Output](#simpsons-13rd-output)
  - [Simpson's 3/8th Rule](#simpsons-38th-rule)
    - [Theory](#simpsons-38th-theory)
    - [Code](#simpsons-38th-code)
    - [Input](#simpsons-38th-input)
    - [Output](#simpsons-38th-output)

- [Solution of Ordinary Differential Equations (ODE)](#solution-of-ordinary-differential-equations-ode)
  - [Interpolation Method (ODE)](#interpolation-method-ode)
    - [Theory](#ode-interpolation-theory)
    - [Code](#ode-interpolation-code)
    - [Input](#ode-interpolation-input)
    - [Output](#ode-interpolation-output)
  - [Runge-Kutta (RK) Method](#runge-kutta-rk-method)
    - [Theory](#rk-method-theory)
    - [Code](#rk-method-code)
    - [Input](#rk-method-input)
    - [Output](#rk-method-output)

---

# Solution of Non-Linear Equations

## Bisection Method

### Bisection Theory


### Bisection Code
```cpp
#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

const ld eps = 1e-9;

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
            ll iter = 0;
            while (u <= v)
            {
                iter++;
                ld mid = (u + v) / 2.0;
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
        cout<<"---------------\n";
    }
    return 0;
}

```

### Bisection Input

```
2
1 2
0 1
-4 0
0.01
3
1 3
0 2
-1 1
-2 0
0.01

```

### Bisection Output

```
Root = -1.9999999999  Iterations = 26
Root = 2.0000000001  Iterations = 26
---------------
Root = 1.5213797069  Iterations = 22
All other roots are imaginary
---------------

```

---

## False Position Method

### False Position Theory

[Add your theory content here]

### False Position Code

```cpp
#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

const ld eps = 1e-9;

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

```

### False Position Input

```
2
1 2
0 1
-4 0
0.01
3
1 3
0 2
-1 1
-2 0
0.01

```

### False Position Output

```
Root = -2.0000000000  Iterations = 1
Root = 2.0000000000  Iterations = 1
---------------
Root = 1.5213797068  Iterations = 4
All other roots are imaginary
---------------

```

---

## Secant Method

### Secant Theory

[Add your theory content here]

### Secant Code

```cpp
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

```

### Secant Input

```
2
1 2
0 1
-4 0
0.01
3
1 3
0 2
-1 1
-2 0
0.01

```

### Secant Output

```
Root = -2.0000000000  Iterations = 2
Root = 2.0000000000  Iterations = 2
---------------
Root = 1.5213797068  Iterations = 4
All other roots are imaginary
---------------

```

---

## Newton Raphson Method

### Newton Raphson Theory

[Add your theory content here]

### Newton Raphson Code

```cpp
#include <bits/stdc++.h>
using namespace std;
#define ll long long
#define ld long double

const ld eps = 1e-9;
//f(x)
ld f(ld x, vector<pair<ld, ll>> &poly)
{
    ld sum = 0;
    for (auto [coeff, power] : poly)
    {
        sum += coeff * pow(x, power);
    }
    return sum;
}

// f'(x)
ld df(ld x, vector<pair<ld, ll>> &poly)
{
    ld sum = 0;
    for (auto [c, p] : poly)
    {
        if (p > 0)
            sum += c * p * pow(x, p - 1);
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
            ld x = u;
            ll iter = 0;

            while (1) 
            {
                iter++;
                ld fx = f(x, poly);
                ld dfx = df(x, poly);

                ld x1 = x - fx / dfx;

                if (fabsl(x1 - x) < eps)
                {
                    cout << fixed << setprecision(10);
                    cout << "Root = " << x1
                         << "  Iterations = " << iter << "\n";
                    break;
                }

                x = x1;
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

```

### Newton Raphson Input

```
2
1 2
0 1
-4 0
0.01
3
1 3
0 2
-1 1
-2 0
0.01

```

### Newton Raphson Output

```
Root = -2.0000000000  Iterations = 1
Root = 2.0000000000  Iterations = 1
---------------
Root = 1.5213797068  Iterations = 3
All other roots are imaginary
---------------

```

---

# Solution of Linear Equations

## Gauss Elimination Method

### Gauss Elimination Theory

Gauss Elimination is a method to solve systems of linear equations by transforming the coefficient matrix into an upper triangular form. Once in triangular form, the unknowns can be solved using back-substitution.



#### Algorithm:



**Step 1:** Write equations in matrix form A × X = B

**Example:**
2x + y - z = 8
-3x - y + 2z = -11
-2x + y + 2z = -3



Matrix form:
A = \[\[2, 1, -1], \[-3, -1, 2], \[-2, 1, 2]]
X = \[\[x], \[y], \[z]]
B = \[\[8], \[-11], \[-3]]



**Step 2:** Transform A into upper triangular form

Use row operations to eliminate variables below the diagonal.



After elimination, the matrix looks like:

\[\[2, 1, -1],
\[0, -0.5, 0.5],
\[0, 0, 1]]



**Step 3:** Back-substitution

Start from the last equation and solve upwards.

Solve for z, then substitute into the second equation to find y, then substitute into the first to find x.



**Final Answer:**
x = 2, y = 3, z = -1



### Gauss Elimination Code
```cpp
#include <bits/stdc++.h>
using namespace std;

int matrixRank(vector<vector<double>> M) {
    int n = M.size();
    int m = M[0].size();
    int rank = 0;
    const double EPS = 1e-9;

    for (int col = 0, row = 0; col < m && row < n; col++) {
        int sel = row;
        for (int i = row; i < n; i++) {
            if (fabs(M[i][col]) > fabs(M[sel][col])) sel = i;
        }
        if (fabs(M[sel][col]) < EPS) continue;

        swap(M[sel], M[row]);
        double div = M[row][col];
        for (int j = col; j < m; j++) M[row][j] /= div;

        for (int i = 0; i < n; i++) {
            if (i != row && fabs(M[i][col]) > EPS) {
                double factor = M[i][col];
                for (int j = col; j < m; j++) {
                    M[i][j] -= factor * M[row][j];
                }
            }
        }
        row++;
        rank++;
    }
    return rank;
}

int checkSystemType(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> Aug(n, vector<double>(n + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) Aug[i][j] = A[i][j];
        Aug[i][n] = b[i];
    }

    int rankA = matrixRank(A);
    int rankAug = matrixRank(Aug);

    if (rankA == rankAug) {
        if (rankA == n) return 0;
        else return 1;            
    } else {
        return 2;                
    }
}

vector<double> gaussElimination(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        for (int k = i + 1; k < n; k++) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    vector<double> x(n);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
    return x;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T; 
    fin >> T;

    for (int t = 1; t <= T; t++) {
        int n; 
        fin >> n;

        vector<vector<double>> A(n, vector<double>(n));
        vector<double> b(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fin >> A[i][j];
            }
        }
        for (int i = 0; i < n; i++) {
            fin >> b[i];
        }

        fout << "Case #" << t << ":\n";

        int systemType = checkSystemType(A, b);
        if (systemType == 2) {
            fout << "System Type: NO SOLUTION\n\n";
            continue;
        } else if (systemType == 1) {
            fout << "System Type: INFINITE SOLUTIONS\n\n";
            continue;
        } else {
            fout << "System Type: UNIQUE SOLUTION\n";
            vector<double> solution = gaussElimination(A, b);
            for (int i = 0; i < n; i++) {
                fout << "x" << i + 1 << " = " << fixed << setprecision(6) << solution[i] << "\n";
            }
            fout << "\n";
        }
    }

    fin.close();
    fout.close();
    return 0;
}

```

### Gauss Elimination Input
```
3
3
2 1 -1
-3 -1 2
-2 1 2
8 -11 -3
2
1 2
2 4
5 10
2
1 2
2 4
5 11
```

### Gauss Elimination Output
```
Case #1:
System Type: UNIQUE SOLUTION
x1 = 2.000000
x2 = 3.000000
x3 = -1.000000

Case #2:
System Type: INFINITE SOLUTIONS

Case #3:
System Type: NO SOLUTION
```

---

## Gauss Jordan Elimination Method

### Gauss Jordan Theory

Gauss-Jordan Elimination is an extension of Gauss Elimination. Instead of stopping at upper triangular form, it continues until the matrix is in reduced row echelon form. This makes the solution appear directly without back-substitution.



#### Algorithm:

**Step 1:** Write equations in matrix form A × X = B

**Example:**
Same system as before:
2x + y - z = 8
-3x - y + 2z = -11
-2x + y + 2z = -3



**Step 2:** Form the augmented matrix \[A|B]

\[\[2, 1, -1 | 8],
\[-3, -1, 2 | -11],
\[-2, 1, 2 | -3]]



**Step 3:** Apply row operations

Make pivot elements = 1

Eliminate all other entries in the pivot column (both above and below)

Final matrix:

\[\[1, 0, 0 | 2],
\[0, 1, 0 | 3],
\[0, 0, 1 | -1]]



**Step 4:** Read off the solution directly
x = 2, y = 3, z = -1



### Gauss Jordan Code

```cpp
#include <bits/stdc++.h>
using namespace std;

const double EPS = 1e-9;

int matrixRank(vector<vector<double>> M) {
    int n = M.size();
    int m = M[0].size();
    int rank = 0;

    for (int col = 0, row = 0; col < m && row < n; col++) {
        int sel = row;
        for (int i = row; i < n; i++) {
            if (fabs(M[i][col]) > fabs(M[sel][col])) sel = i;
        }
        if (fabs(M[sel][col]) < EPS) continue;

        swap(M[sel], M[row]);
        double div = M[row][col];
        for (int j = col; j < m; j++) M[row][j] /= div;

        for (int i = 0; i < n; i++) {
            if (i != row && fabs(M[i][col]) > EPS) {
                double factor = M[i][col];
                for (int j = col; j < m; j++) {
                    M[i][j] -= factor * M[row][j];
                }
            }
        }
        row++;
        rank++;
    }
    return rank;
}

int checkSystemType(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    vector<vector<double>> Aug(n, vector<double>(n + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) Aug[i][j] = A[i][j];
        Aug[i][n] = b[i];
    }

    int rankA = matrixRank(A);
    int rankAug = matrixRank(Aug);

    if (rankA == rankAug) {
        if (rankA == n) return 0;
        else return 1;        
    } else {
        return 2;              
    }
}

vector<double> gaussJordan(vector<vector<double>> A, vector<double> b) {
    int n = A.size();

    for (int i = 0; i < n; i++) {
        int maxRow = i;
        for (int k = i + 1; k < n; k++) {
            if (fabs(A[k][i]) > fabs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        swap(A[i], A[maxRow]);
        swap(b[i], b[maxRow]);

        if (fabs(A[i][i]) < EPS) {
            throw runtime_error("Singular matrix or no unique solution.");
        }

        double pivot = A[i][i];
        for (int j = 0; j < n; j++) {
            A[i][j] /= pivot;
        }
        b[i] /= pivot;

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; j++) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    }
    return b;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    ifstream fin("input.txt");
    ofstream fout("output.txt");

    int T; 
    fin >> T;

    for (int t = 1; t <= T; t++) {
        int n; 
        fin >> n;

        vector<vector<double>> A(n, vector<double>(n));
        vector<double> b(n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fin >> A[i][j];
            }
        }
        for (int i = 0; i < n; i++) {
            fin >> b[i];
        }

        fout << "Case #" << t << ":\n";

        int systemType = checkSystemType(A, b);
        if (systemType == 2) {
            fout << "System Type: NO SOLUTION\n\n";
            continue;
        } else if (systemType == 1) {
            fout << "System Type: INFINITE SOLUTIONS\n\n";
            continue;
        } else {
            fout << "System Type: UNIQUE SOLUTION\n";
            try {
                vector<double> solution = gaussJordan(A, b);
                for (int i = 0; i < n; i++) {
                    fout << "x" << i + 1 << " = " 
                         << fixed << setprecision(6) << solution[i] << "\n";
                }
            } catch (runtime_error &e) {
                fout << "Error: " << e.what() << "\n";
            }
            fout << "\n";
        }
    }

    fin.close();
    fout.close();
    return 0;
}

```

### Gauss Jordan Input

```
3
3
2 1 -1
-3 -1 2
-2 1 2
8 -11 -3
2
1 2
2 4
5 10
2
1 2
2 4
5 11
```

### Gauss Jordan Output

```
Case #1:
System Type: UNIQUE SOLUTION
x1 = 2.000000
x2 = 3.000000
x3 = -1.000000

Case #2:
System Type: INFINITE SOLUTIONS

Case #3:
System Type: NO SOLUTION

```

---

## LU Decomposition Method

### LU Decomposition Theory

LU Decomposition breaks a matrix A into the product of two matrices:
A = L × U
**Where:**

L = Lower triangular matrix

U = Upper triangular matrix

This makes solving systems of equations easier and faster, especially for multiple right-hand sides.

#### Algorithm:

**Step 1:** Write equations in matrix form A × X = B

**Example:**
2x + y - z = 8
-3x - y + 2z = -11
-2x + y + 2z = -3

Matrix form:
A = \[\[2, 1, -1], \[-3, -1, 2], \[-2, 1, 2]]
X = \[\[x], \[y], \[z]]
B = \[\[8], \[-11], \[-3]]

**Step 2:** Decompose A into L and U

L is lower triangular (entries below diagonal)

U is upper triangular (entries above diagonal)

For this example:
L = \[\[1, 0, 0], \[-1.5, 1, 0], \[-1, -1, 1]]
U = \[\[2, 1, -1], \[0, 0.5, 0.5], \[0, 0, 1]]

**Step 3:** Solve L × Y = B (forward substitution)

Find intermediate vector Y.

**Step 4:** Solve U × X = Y (back substitution)

Find final solution vector X.

Final Answer:
x = 2, y = 3, z = -1



### LU Decomposition Code

```cpp
#include <bits/stdc++.h>
using namespace std;

const double EPSILON = 1e-9;

bool isZero(double val) {
    return fabs(val) < EPSILON;
}

bool luDecompose(vector<vector<double>>& A, vector<int>& P, int n) {
    for (int i = 0; i < n; i++) {
        P[i] = i;
    }
    
    for (int k = 0; k < n; k++) {
        double maxVal = 0;
        int pivotRow = k;
        for (int i = k; i < n; i++) {
            if (fabs(A[i][k]) > maxVal) {
                maxVal = fabs(A[i][k]);
                pivotRow = i;
            }
        }

        if (pivotRow != k) {
            swap(A[k], A[pivotRow]);
            swap(P[k], P[pivotRow]);
        }
        
        if (isZero(A[k][k])) {
            return false;
        }
        
        for (int i = k + 1; i < n; i++) {
            A[i][k] = A[i][k] / A[k][k];
            for (int j = k + 1; j < n; j++) {
                A[i][j] = A[i][j] - A[i][k] * A[k][j];
            }
        }
    }
    return true;
}


void forwardSubstitution(const vector<vector<double>>& LU, const vector<int>& P, 
                         const vector<double>& b, vector<double>& y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = b[P[i]];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j] * y[j];
        }
    }
}

void backwardSubstitution(const vector<vector<double>>& LU, 
                          const vector<double>& y, vector<double>& x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= LU[i][j] * x[j];
        }
        x[i] = x[i] / LU[i][i];
    }
}

int checkSystemType(vector<vector<double>> A, vector<double> b, int n) {
    vector<int> P(n);
    vector<vector<double>> augmented = A;
    
    for (int i = 0; i < n; i++) {
        augmented[i].push_back(b[i]);
    }
    
    int rank_A = 0, rank_Ab = 0;
    
    for (int col = 0; col < n; col++) {
        int pivotRow = -1;
        for (int row = rank_A; row < n; row++) {
            if (! isZero(A[row][col])) {
                pivotRow = row;
                break;
            }
        }
        
        if (pivotRow == -1) continue;
        
        if (pivotRow != rank_A) {
            swap(A[rank_A], A[pivotRow]);
            swap(augmented[rank_A], augmented[pivotRow]);
        }

        for (int row = rank_A + 1; row < n; row++) {
            if (! isZero(A[row][col])) {
                double factor = A[row][col] / A[rank_A][col];
                for (int j = 0; j < n; j++) {
                    A[row][j] -= factor * A[rank_A][j];
                }
                for (int j = 0; j <= n; j++) {
                    augmented[row][j] -= factor * augmented[rank_A][j];
                }
            }
        }
        rank_A++;
    }
    
    rank_Ab = rank_A;
    for (int i = rank_A; i < n; i++) {
        if (!isZero(augmented[i][n])) {
            return 2; 
        }
    }
    
    if (rank_A < n) {
        return 1;
    }
    
    return 0;
}

void solveSystem(ifstream& input, ofstream& output, int testCase) {
    int n;
    input >> n;
    
    vector<vector<double>> A(n, vector<double>(n));
    vector<double> b(n), x(n), y(n);
    vector<int> P(n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            input >> A[i][j];
        }
    }
    
    for (int i = 0; i < n; i++) {
        input >> b[i];
    }
    
    output << "========== TEST CASE " << testCase << " ==========\n";
    output << "System size: " << n << "x" << n << "\n\n";
    
    int systemType = checkSystemType(A, b, n);
    
    if (systemType == 2) {
        output << "System Type: NO SOLUTION\n";
        output << "The system is inconsistent.\n\n";
        return;
    } else if (systemType == 1) {
        output << "System Type: INFINITE SOLUTIONS\n";
        output << "The system has infinitely many solutions.\n\n";
        return;
    }
    
    output << "System Type:  UNIQUE SOLUTION\n\n";
    vector<vector<double>> LU = A;
    if (!luDecompose(LU, P, n)) {
        output << "Error: Matrix is singular\n\n";
        return;
    }
    
    forwardSubstitution(LU, P, b, y, n);
    backwardSubstitution(LU, y, x, n);
    
    output << "Solution:\n";
    for (int i = 0; i < n; i++) {
        output << "x[" << i << "] = " << fixed << setprecision(6) << x[i] << "\n";
    }
    output << "\n";
}

int main() {
    ifstream input("input.txt");
    ofstream output("output.txt");
    
    if (!input.is_open() || !output.is_open()) {
        cerr << "Error opening files!\n";
        return 1;
    }
    
    int numTests;
    input >> numTests;
    
    output << "LU DECOMPOSITION - LINEAR EQUATION SOLVER\n";
    output << "==========================================\n\n";
    
    for (int i = 1; i <= numTests; i++) {
        solveSystem(input, output, i);
    }
    
    input.close();
    output.close();
    
    cout << "LU decomposition completed.  Check lu_output.txt\n";
    return 0;
}

```

### LU Decomposition Input

```
4
3
2 1 -1
-3 -1 2
-2 1 2
8
-11
-3
2
1 2
3 4
5
6
3
2 4 6
1 2 3
3 6 9
1
2
3
3
1 0 0
0 1 0
0 0 0
1
2
3

```

### LU Decomposition Output

```
LU DECOMPOSITION - LINEAR EQUATION SOLVER
==========================================

========== TEST CASE 1 ==========
System size: 3x3

System Type:  UNIQUE SOLUTION

Solution:
x[0] = 2.000000
x[1] = 3.000000
x[2] = -1.000000

========== TEST CASE 2 ==========
System size: 2x2

System Type:  UNIQUE SOLUTION

Solution:
x[0] = -4.000000
x[1] = 4.500000

========== TEST CASE 3 ==========
System size: 3x3

System Type: NO SOLUTION
The system is inconsistent.

========== TEST CASE 4 ==========
System size: 3x3

System Type: NO SOLUTION
The system is inconsistent.



```

---
## Matrix Inversion Method

### Matrix Inversion Theory


For solving a linear equation we can implement: If A × X = B, then X = A⁻¹ × B



Here:



A = coefficient matrix



X = vector of unknowns



B = constants



A⁻¹ = inverse of matrix A



### **Algorithm:**



**Step 1:** Write equations in matrix form A × X = B



**Example:**

2x + y = 5

x + 3y = 10



Matrix form:

A = \[\[2, 1], \[1, 3]]

X = \[\[x], \[y]]

B = \[\[5], \[10]]



**Step 2:** Find the inverse of A (A⁻¹)



For a 2×2 matrix:

A⁻¹ = (1 / (ad - bc)) × \[\[d, -b], \[-c, a]]



**Here:**

A = \[\[2, 1], \[1, 3]]

Determinant = (2×3) - (1×1) = 6 - 1 = 5

So A⁻¹ = (1/5) × \[\[3, -1], \[-1, 2]]



**Step 3:** Multiply both sides by A⁻¹



X = A⁻¹ × B



X = (1/5) × \[\[3, -1], \[-1, 2]] × \[\[5], \[10]]



**Step 4:** Do the multiplication



X = (1/5) × \[\[(3×5) + (-1×10)], \[(-1×5) + (2×10)]]

X = (1/5) × \[\[15 - 10], \[-5 + 20]]

X = (1/5) × \[\[5], \[15]]

X = \[\[1], \[3]]



Final Answer:

x = 1

y = 3



### Matrix Inversion Code

```cpp
#include <bits/stdc++.h>
#define int long long
const double EPSILON = 1e-9;
using namespace std;

void get_submatrix(vector<vector<double>> mat, vector<vector<double>> &temp, int p, int q, int n)
{
    int r = 0;
    for (int i = 0; i < n; i++){
        if (i != p){
            int s = 0;
            for (int j = 0; j < n; j++){
                if (j != q){
                    temp[r][s++] = mat[i][j];
                }
            }

            r++;
        }
    }
}

double get_determinant(vector<vector<double>> mat, int n)
{
    if (n == 1){
        return mat[0][0];
    }

    if (n == 2){
        return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
    }

    int D = 0;
    int sign = 1;
    vector<vector<double>> temp(n - 1, vector<double>(n - 1, 0));
    for (int c = 0; c < n; c++){
        get_submatrix(mat, temp, 0, c, n);

        D += sign * mat[0][c] * get_determinant(temp, n - 1);
        sign *= -1;
    }

    return D;
}

bool get_inverse(vector<vector<double>> mat, vector<vector<double>> &inverse, int n)
{
    double A = get_determinant(mat, n);
    if (abs(A) < EPSILON){
        return false;
    }

    if (n == 1){
        inverse[0][0] = 1.0 / A;
        return true;
    }
    double inv_det = 1.0 / A;
    vector<vector<double>> temp(n - 1, vector<double>(n - 1, 0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            get_submatrix(mat, temp, i, j, n);

            double minor = get_determinant(temp, n - 1);

            double sign = ((i + j) % 2 == 0 ? 1.0 : -1.0);

            double cofactor = sign * minor;

            inverse[j][i] = cofactor * inv_det;
        }
    }

    return true;
}

void print_matrix(vector<vector<double>> &m, int n){
    cout << fixed << setprecision(3);
    for (int i = 0; i < n; i++){
        cout << "|";
        for (int j = 0; j < n; j++)
        {
            cout << setw(4) << m[i][j] << " ";
        }
        cout << "|\n";
    }
}

int32_t main(){

    ifstream inputFile("input.txt");
    if (!inputFile.is_open()){
        cerr << "Error while opening input file" << endl;
    }

    ofstream outputFile("output.txt", ios::trunc);
    if (!inputFile.is_open()){
        cerr << "Error while opening input file" << endl;
    }

    int n;
    // cin>>n;
    if (!(inputFile >> n)){
        cerr << "Error in reading matrix dimension N" << endl;
        inputFile.close();
        outputFile.close();
        return 1;
    }

    vector<vector<double>> matrix(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            if (!(inputFile >> matrix[i][j])){
                cerr << "Error in reading matrix element at (" << i << "," << j << ")" << endl;
                inputFile.close();
                outputFile.close();
                return 1;
            }
        }
    }
    vector<vector<double>> inverse(n, vector<double>(n, 0));
    get_inverse(matrix, inverse, n);

    outputFile << fixed << setprecision(3);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            outputFile << setw(6) << inverse[i][j] << "\t";
        }
        outputFile << "\n";
    }
    return 0;
}

```

### Matrix Inversion Input

```
4
1 0 0 0
0 2 0 0
0 0 3 0
0 0 0 4

```

### Matrix Inversion Output

```
 1.000	-0.000	 0.000	-0.000	
-0.000	 0.500	-0.000	 0.000	
 0.000	-0.000	 0.333	-0.000	
-0.000	 0.000	-0.000	 0.250	

```

---

# Interpolation

## Newton's Forward Interpolation

### Forward Interpolation Theory

[Add your theory content here]

### Forward Interpolation Code

```cpp
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


```

### Forward Interpolation Input

```
4
45 0.7071
50 0.766
55 0.8192
60 0.866
52

```

### Forward Interpolation Output

```
45.0000	0.7071	0.0589	-0.0057	-0.0007	
50.0000	0.7660	0.0532	-0.0064	0.0000	
55.0000	0.8192	0.0468	0.0000	0.0000	
60.0000	0.8660	0.0000	0.0000	0.0000	
Value at 52.0000 = 0.7896

```

---

## Newton's Backward Interpolation

### Backward Interpolation Theory

[Add your theory content here]

### Backward Interpolation Code

```cpp
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

```

### Backward Interpolation Input

```
4
5 12
6 13
9 14
11 16
7

```

### Backward Interpolation Output

```
5.0000	12.0000	0.0000	0.0000	0.0000	
6.0000	13.0000	0.0000	0.0000	0.0000	
9.0000	14.0000	0.0000	0.0000	0.0000	
11.0000	16.0000	0.0000	0.0000	0.0000	
Interpole value at 7.0000 = 16.0000

```

---

## Divided Difference Interpolation

### Divided Difference Theory

[Add your theory content here]

### Divided Difference Code

```cpp
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

```

### Divided Difference Input

```
4
300	2.4771
304	2.4829
305	2.4843
307	2.4871
301

```

### Divided Difference Output

```
300.000000	2.477100	0.001450	-0.000010	0.000001	
304.000000	2.482900	0.001400	-0.000000	0.000000	
305.000000	2.484300	0.001400	0.000000	0.000000	
307.000000	2.487100	0.000000	0.000000	0.000000	
Value at 301.000000 : 2.478597

```

---
#Numerical Differentiation

###Differentiation by Forward Interpolation

###Differentiation Forward Theory

[Add your theory content here]

###Differentiation Forward Code
```cpp
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
```

###Differentiation Forward Input
```
5 2
1 0 3
0 2 1
4 2
1 0 3
0 2 1
```
###Differentiation Forward Output
```

TEST CASE #1
n = 5, degree = 2, a = 0, b = 2, X = 1
Difference Table:
           1         0.48         0.96            0 -8.88178e-16 
        1.48         1.44         0.96 -8.88178e-16            0 
        2.92          2.4         0.96            0            0 
        5.32         3.36            0            0            0 
        8.68            0            0            0            0 
Numerical y'  = 6.0000000000
Numerical y'' = 6.0000000000
Exact y'      = 6.0000000000
Exact y''     = 6.0000000000
First derivative error  = 0.0000000000%
Second derivative error = 0.0000000000%
---------------

TEST CASE #2
n = 4, degree = 2, a = 0.0000000000, b = 2.0000000000, X = 1.0000000000
Difference Table:
1.0000000000 0.7500000000 1.5000000000 0.0000000000 
1.7500000000 2.2500000000 1.5000000000 0.0000000000 
4.0000000000 3.7500000000 0.0000000000 0.0000000000 
7.7500000000 0.0000000000 0.0000000000 0.0000000000 
Numerical y'  = 6.0000000000
Numerical y'' = 6.0000000000
Exact y'      = 6.0000000000
Exact y''     = 6.0000000000
First derivative error  = 0.0000000000%
Second derivative error = 0.0000000000%
---------------
```
###Differentiation by Backward Interpolation Method

###Differentiation Backward Theory

[Add your theory content here]

###Differentiation Backward Code
```cpp
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
```

###Differentiation Backward Input
```
5 2
1 0 2
0 1 2 3 4
2
4 3
0 1 0 1
0 1 2 3 4
1.5
3 1
2 3 1
0 0.5 1.0
0.5

```
###Differentiation Backward Output
```
TEST CASE #1
Polynomial: f(x) = 1.0000 + 2.0000x^2
Number of points: 5
x-values: 0 1 2 3 4 
y-values: 1.000000 3.000000 9.000000 19.000000 33.000000 
Step size (h): 1.000000
Differentiation point: 2.000000
Backward Difference Table:
Row 0:     1.000000 
Row 1:     3.000000     2.000000 
Row 2:     9.000000     6.000000     4.000000 
Row 3:    19.000000    10.000000     4.000000     0.000000 
Row 4:    33.000000    14.000000     4.000000     0.000000     0.000000 
Approximate derivative: 8.000000
True derivative: 8.000000
Percentage error: 0.000000 %
---------------
TEST CASE #2
Polynomial: f(x) = 1.0000x + 1.0000x^3
Number of points: 4
x-values: 0.000000 1.000000 2.000000 3.000000 
y-values: 0.000000 2.000000 10.000000 30.000000 
Step size (h): 1.000000
Differentiation point: 4.000000
Backward Difference Table:
Row 0:     0.000000 
Row 1:     2.000000     2.000000 
Row 2:    10.000000     8.000000     6.000000 
Row 3:    30.000000    20.000000    12.000000     6.000000 
Approximate derivative: 49.000000
True derivative: 49.000000
Percentage error: 0.000000 %
---------------

```
# Regression (Curve Fitting)

## Linear Regression

### Linear Regression Theory

[Add your theory content here]

### Linear Regression Code

```cpp
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

```

### Linear Regression Input

```
5
1 2
2 3
3 5
4 4
5 6
3
1 2
2 4
3 6

```

### Linear Regression Output

```
Linear Regression: y = 0.9000000000 * x + 1.3000000000
---------------
Linear Regression: y = 2.0000000000 * x + 0.0000000000
---------------

```

---

## Transcendental Equation Regression

### Transcendental Regression Theory

[Add your theory content here]

### Transcendental Regression Code

```cpp
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

```

### Transcendental Regression Input

```
5
0 1
1 2.718
2 7.389
3 20.085
4 54.598
3
0 2
1 4
2 8

```

### Transcendental Regression Output

```
Exponential Regression: y = 0.9999575583 * e^(1.0000071456*x)
---------------
Exponential Regression: y = 2.0000000000 * e^(0.6931471806*x)
---------------

```

---

## Polynomial Regression

### Polynomial Regression Theory

[Add your theory content here]

### Polynomial Regression Code

```cpp
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

```

### Polynomial Regression Input

```
5 2
1 1
2 4
3 9
4 16
5 25
3 1
1 2
2 4
3 6

```

### Polynomial Regression Output

```
Polynomial Regression: y = 0.0000000000*x^0 + -0.0000000000*x^1 + 1.0000000000*x^2
---------------
Polynomial Regression: y = 0.0000000000*x^0 + 2.0000000000*x^1
---------------

```

---

# Numerical Integration

## Simpson's 1/3rd Rule

### Simpson's 1/3rd Theory

Simpson’s 1/3 Rule is a numerical integration method. It approximates the area under a curve by dividing the interval into an even number of subintervals through the points.

#### Algorithm:

**Step 1:** Divide the interval \[a, b] into n subintervals, where n must be even.
**Step 2:** Compute the width of each subinterval:
h = (b - a) / n
**Step 3:** Apply the formula:

Integral = (h/3) \* \[f(x0) + 4f(x1+x3+x5+....) + 2f(x2+x4+x6+....) + f(xn)]

**Step 4:** Add up the terms to get the approximate value of the integral.

### Simpson's 1/3rd Code

```cpp
#include <bits/stdc++.h>

using namespace std;


double f1(double x) {
    return x * x;  
}

double f2(double x) {
    return 1.0 / (1.0 + x * x);  
}

double f3(double x) {
    return sin(x);  
}

double f4(double x) {
    return exp(x);  
}

double f5(double x) {
    return sqrt(x);  
}


double simpsonsOneThird(double (*f)(double), double a, double b, int n) {
    if (n % 2 != 0) {
        n++;
    }
    
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; i += 2) {
        sum += 4 * f(a + i * h);
    }
    
    for (int i = 2; i < n; i += 2) {
        sum += 2 * f(a + i * h);
    }
    
    return (h / 3.0) * sum;
}

int main() {
    ifstream input("input.txt");
    ofstream output("output.txt");
    
    if (!input.is_open() || !output.is_open()) {
        cerr << "Error opening files!\n";
        return 1;
    }
    
    output << "SIMPSON'S 1/3 RULE INTEGRATION\n";
    output << "===============================\n\n";
    
    int numTests;
    input >> numTests;
    
    for (int t = 1; t <= numTests; t++) {
        int funcNum, n;
        double a, b;
        
        input >> funcNum >> a >> b >> n;
        
        output << "Test Case " << t << ":\n";
        output << "Function: f" << funcNum << "\n";
        output << "Limits: [" << a << ", " << b << "]\n";
        output << "Intervals: " << n << "\n";
        
        double result;
        switch(funcNum) {
            case 1: result = simpsonsOneThird(f1, a, b, n); 
                    output << "f(x) = x^2\n"; break;
            case 2: result = simpsonsOneThird(f2, a, b, n);
                    output << "f(x) = 1/(1+x^2)\n"; break;
            case 3: result = simpsonsOneThird(f3, a, b, n);
                    output << "f(x) = sin(x)\n"; break;
            case 4: result = simpsonsOneThird(f4, a, b, n);
                    output << "f(x) = e^x\n"; break;
            case 5: result = simpsonsOneThird(f5, a, b, n);
                    output << "f(x) = sqrt(x)\n"; break;
            default: 
                output << "Invalid function number\n\n";
                continue;
        }
        
        output << "Result: " << fixed << setprecision(8) << result << "\n\n";
    }
    
    input.close();
    output.close();
    
    cout << "Simpson's 1/3 rule completed. Check simpson_one_third_output.txt\n";
    return 0;
}

```

### Simpson's 1/3rd Input

```
6
1 0 1 6
2 0 1 9
3 0 3.14159265 9
4 0 1 12
5 1 4 9
6 0 2 6

```

### Simpson's 1/3rd Output

```SIMPSON'S 1/3 RULE INTEGRATION
===============================

Test Case 1:
Function: f1
Limits: [0, 1]
Intervals: 6
f(x) = x^2
Result: 0.33333333

Test Case 2:
Function: f2
Limits: [0.00000000, 1.00000000]
Intervals: 9
f(x) = 1/(1+x^2)
Result: 0.78539815

Test Case 3:
Function: f3
Limits: [0.00000000, 3.14159265]
Intervals: 9
f(x) = sin(x)
Result: 2.00010952

Test Case 4:
Function: f4
Limits: [0.00000000, 1.00000000]
Intervals: 12
f(x) = e^x
Result: 1.71828229

Test Case 5:
Function: f5
Limits: [1.00000000, 4.00000000]
Intervals: 9
f(x) = sqrt(x)
Result: 4.66665163

Test Case 6:
Function: f6
Limits: [0.00000000, 2.00000000]
Intervals: 6
Invalid function number

```

---

## Simpson's 3/8th Rule

### Simpson's 3/8th Theory

Simpson’s 3/8 Rule is another numerical integration method. It approximates the area under a curve by dividing the interval into subintervals that are multiples of 3 through the points.



#### Algorithm:

**Step 1:** Divide the interval \[a, b] into n subintervals, where n must be a multiple of 3.


**Step 2:** Compute the width of each subinterval:
h = (b - a) / n


**Step 3:** Apply the formula:

Integral = (3h/8) \* \[f(x0) + 3f(x1+x2+x4+x5+...) + 2f(x3+x6+x9+......) + f(xn)]



**Step 4:** Add up the terms to get the approximate value of the integral.



### Simpson's 3/8th Code

```cpp
#include <bits/stdc++.h>

using namespace std;

double f1(double x) {
    return x * x;  
}

double f2(double x) {
    return 1.0 / (1.0 + x * x);
}

double f3(double x) {
    return sin(x); 
}

double f4(double x) {
    return exp(x);  
}

double f5(double x) {
    return sqrt(x);
}

double f6(double x) {
    return x * x * x; 
}


double simpsonsThreeEighth(double (*f)(double), double a, double b, int n) {
    while (n % 3 != 0) {
        n++;
    }
    
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    
    for (int i = 1; i < n; i++) {
        if (i % 3 == 0) {
            sum += 2 * f(a + i * h);
        } else {
            sum += 3 * f(a + i * h);
        }
    }
    
    return (3.0 * h / 8.0) * sum;
}

int main() {
    ifstream input("input.txt");
    ofstream output("output.txt");
    
    if (!input.is_open() || !output.is_open()) {
        cerr << "Error opening files!\n";
        return 1;
    }
    
    output << "SIMPSON'S 3/8 RULE INTEGRATION\n";
    output << "===============================\n\n";
    
    int numTests;
    input >> numTests;
    
    for (int t = 1; t <= numTests; t++) {
        int funcNum, n;
        double a, b;
        
        input >> funcNum >> a >> b >> n;
        
        output << "Test Case " << t << ":\n";
        output << "Function: f" << funcNum << "\n";
        output << "Limits:  [" << a << ", " << b << "]\n";
        output << "Intervals: " << n << "\n";
        
        double result;
        switch(funcNum) {
            case 1: result = simpsonsThreeEighth(f1, a, b, n); 
                    output << "f(x) = x^2\n"; break;
            case 2: result = simpsonsThreeEighth(f2, a, b, n);
                    output << "f(x) = 1/(1+x^2)\n"; break;
            case 3: result = simpsonsThreeEighth(f3, a, b, n);
                    output << "f(x) = sin(x)\n"; break;
            case 4: result = simpsonsThreeEighth(f4, a, b, n);
                    output << "f(x) = e^x\n"; break;
            case 5: result = simpsonsThreeEighth(f5, a, b, n);
                    output << "f(x) = sqrt(x)\n"; break;
            case 6: result = simpsonsThreeEighth(f6, a, b, n);
                    output << "f(x) = x^3\n"; break;
            default: 
                output << "Invalid function number\n\n";
                continue;
        }
        
        output << "Result: " << fixed << setprecision(8) << result << "\n\n";
    }
    
    input.close();
    output.close();
    
    cout << "Simpson's 3/8 rule completed. Check simpson_three_eighth_output. txt\n";
    return 0;
}

```

### Simpson's 3/8th Input

```
6
1 0 1 6
2 0 1 9
3 0 3.14159265 9
4 0 1 12
5 1 4 9
6 0 2 6

```

### Simpson's 3/8th Output

```
SIMPSON'S 3/8 RULE INTEGRATION
===============================

Test Case 1:
Function: f1
Limits:  [0, 1]
Intervals: 6
f(x) = x^2
Result: 0.33333333

Test Case 2:
Function: f2
Limits:  [0.00000000, 1.00000000]
Intervals: 9
f(x) = 1/(1+x^2)
Result: 0.78539808

Test Case 3:
Function: f3
Limits:  [0.00000000, 3.14159265]
Intervals: 9
f(x) = sin(x)
Result: 2.00038224

Test Case 4:
Function: f4
Limits:  [0.00000000, 1.00000000]
Intervals: 12
f(x) = e^x
Result: 1.71828286

Test Case 5:
Function: f5
Limits:  [1.00000000, 4.00000000]
Intervals: 9
f(x) = sqrt(x)
Result: 4.66661964

Test Case 6:
Function: f6
Limits:  [0.00000000, 2.00000000]
Intervals: 6
f(x) = x^3
Result: 4.00000000



```

---

# Solution of Ordinary Differential Equations (ODE)

## Interpolation Method (ODE)

### ODE Interpolation Theory

[Add your theory content here]

### ODE Interpolation Code

```python
# Add your code here

```

### ODE Interpolation Input

```
[Add your input format here]

```

### ODE Interpolation Output

```
[Add your output format here]

```

---

## Runge-Kutta (RK) Method

### RK Method Theory

Runge-Kutta methods are used to solve ordinary differential equations (ODEs). The most common is the 4th-order Runge-Kutta method (RK4), which gives a very accurate solution.

#### 

#### Algorithm:

**Step 1:** Start with an initial condition (x0, y0).


**Step 2:** Choose a step size h.


**Step 3:** Compute four slopes:

k1 = h \* f(x0, y0)
k2 = h \* f(x0 + h/2, y0 + k1/2)
k3 = h \* f(x0 + h/2, y0 + k2/2)
k4 = h \* f(x0 + h, y0 + k3)



**Step 4:** Update the solution:
y1 = y0 + (1/6)(k1 + 2k2 + 2k3 + k4)



**Step 5:** Repeat for the next step until the desired x is reached.



### RK Method Code

```cpp
#include<bits/stdc++.h>
using namespace std;

double dydx(double x, double y){
    return (x-y)/2.0;
}

double rkmethod(double x0, double y0, double x, double h){

    int n = (x-x0)/h;
    double k1, k2, k3, k4;
    double y = y0;
    for(int i =0;i<n;i++){
        k1= h*dydx(x0,y);
        k2 = h*dydx(x0+0.5*h, y+0.5*k1);
        k3 = h*dydx(x0+0.5*h, y+0.5*k2);
        k4 = h*dydx(x0+h, y+k3);

        y = y + (1.0/6.0)*(k1+2*k2+2*k3+k4);

        x0 =x0+h;
    }

    return y;
}

int32_t main(){
    double x0 = 0.0, y0 = 1.0, x =2.0,  h = 0.2;
    ifstream inputFile("input.txt");
    ofstream outputFile("output.txt", ios::trunc);

    if(!inputFile.is_open()){
        cerr<<"Error while opening input file"<<endl;
        return 1;
    }
    if(!(inputFile>>x0>>y0>>x>>h)){
        cerr<<"Error reading input values"<<endl;
        inputFile.close();
        return 1;
    }

    inputFile.close();
    if(!outputFile.is_open()){
        cerr<<"Error while opening output file"<<endl;
        return 1;
    }
    outputFile<<fixed<<setprecision(6);
    outputFile<<"The value of y at "<<x<<": "<<rkmethod(x0,y0,x,h)<<endl;
    outputFile.close();
    
    return 0;
}
```

### RK Method Input

```
0.0 1.0 2.0 0.2

```

### RK Method Output

```
The value of y at 2.000000: 1.019710

```

```

```
