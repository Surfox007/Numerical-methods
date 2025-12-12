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
