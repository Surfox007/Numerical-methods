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