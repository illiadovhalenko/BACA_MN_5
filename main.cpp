//Illia Dovhalenko
#include <iostream>
#include<cmath>
#include <iomanip>
#include "funkcja.h"
using namespace std;
struct Jet{
public:
    double tab[6];
    Jet(){
        std::fill(tab, tab+6, 0);
    }
    Jet(double a){
        std::fill(tab, tab+6, 0);
        tab[0]=a;
    }
    Jet(const double x, const int poz){
        std::fill(tab, tab+6, 0);
        tab[0]=x;
        tab[poz]=1;
    }
    Jet(const Jet& a){
        std::copy(a.tab, a.tab+6, this->tab);
    }
    friend Jet operator+(const Jet&x, const Jet& a){
        Jet temp;
        for(int i=0;i<6; i++){
            temp.tab[i]=x.tab[i]+a.tab[i];
        }
        return temp;
    }
    Jet operator+(const double a){
        Jet temp(*this);
        temp.tab[0]+=a;
        return temp;
    }
    friend Jet operator-(const Jet&x, const Jet& a){
        Jet temp;
        for(int i=0;i<6; i++){
            temp.tab[i]=x.tab[i]-a.tab[i];
        }
        return temp;
    }
    Jet operator-(const double a){
        Jet temp(*this);
        temp.tab[0]-=a;
        return temp;
    }
    friend Jet operator-(const Jet&a){
        Jet temp;
        for(int i=0; i<6; i++)
            temp.tab[i]=-a.tab[i];
        return temp;
    }
    Jet operator*(double a){
        Jet temp;
        for(int i=0; i<6; i++){
            temp.tab[i]=tab[i]*a;
        }
        return temp;
    }
    friend Jet operator*(const Jet&x, const Jet& a){
        Jet temp;
        temp.tab[0]= x.tab[0]*a.tab[0];
        temp.tab[1]=x.tab[1]*a.tab[0]+a.tab[1]*x.tab[0];
        temp.tab[2]=x.tab[2]*a.tab[0]+a.tab[2]*x.tab[0];
        temp.tab[3]=x.tab[0]*a.tab[3]+2* x.tab[1]*a.tab[1]+x.tab[3]*a.tab[0];
        temp.tab[4]=x.tab[0]*a.tab[4]+a.tab[1]*x.tab[2]+a.tab[2]*x.tab[1]+x.tab[4]*a.tab[0];
        temp.tab[5]=x.tab[0]*a.tab[5]+2* x.tab[2]*a.tab[2]+x.tab[5]*a.tab[0];
        return temp;
    }
   friend Jet operator/(const Jet& x, const Jet& a){
        Jet temp;
        temp.tab[0]=x.tab[0]/a.tab[0];
        temp.tab[1]=x.tab[1]/a.tab[0]-(a.tab[1]*x.tab[0])/std::pow(a.tab[0], 2);
        temp.tab[2]=x.tab[2]/a.tab[0]-(a.tab[2]*x.tab[0])/std::pow(a.tab[0], 2);
        temp.tab[3]=(-a.tab[0]*(2*x.tab[1]*a.tab[1]+x.tab[0]*a.tab[3])+x.tab[3]*pow(a.tab[0], 2)+2*x.tab[0]*pow(a.tab[1], 2))/pow(a.tab[0], 3);
        temp.tab[4]=(-a.tab[0]*(x.tab[1]*a.tab[2]+x.tab[2]*a.tab[1]+x.tab[0]*a.tab[4])+x.tab[4]*pow(a.tab[0], 2)+2*x.tab[0]*a.tab[1]*a.tab[2])/pow(a.tab[0], 3);
        temp.tab[5]=(-a.tab[0]*(2*x.tab[2]*a.tab[2]+x.tab[0]*a.tab[5])+x.tab[5]*pow(a.tab[0], 2)+2*x.tab[0]*pow(a.tab[2], 2))/pow(a.tab[0], 3);
        return temp;
    }
   Jet operator/(const double a){
        Jet temp(*this);
        for(int i=0; i<6; i++){
            temp.tab[i]/=a;
        }
        return temp;
    }
    friend Jet sin(const Jet& a){
        Jet temp;
        temp.tab[0]= std::sin(a.tab[0]);
        temp.tab[1]=a.tab[1]*std::cos(a.tab[0]);
        temp.tab[2]=a.tab[2]*std::cos(a.tab[0]);
        temp.tab[3]=a.tab[3]*std::cos(a.tab[0])-pow(a.tab[1], 2)* std::sin(a.tab[0]);
        temp.tab[4]=a.tab[4]*std::cos(a.tab[0])-a.tab[1]*a.tab[2]*std::sin(a.tab[0]);
        temp.tab[5]=a.tab[5]*std::cos(a.tab[0])- pow(a.tab[2], 2)*std::sin(a.tab[0]);
        return temp;
    }
    friend Jet cos(const Jet& a){
        Jet temp;
        temp.tab[0]= std::cos(a.tab[0]);
        temp.tab[1]=-a.tab[1]*std::sin(a.tab[0]);
        temp.tab[2]=-a.tab[2]*std::sin(a.tab[0]);
        temp.tab[3]=-a.tab[3]*std::sin(a.tab[0])-pow(a.tab[1], 2)*std::cos(a.tab[0]);
        temp.tab[4]=-a.tab[4]*std::sin(a.tab[0])-a.tab[1]*a.tab[2]*std::cos(a.tab[0]);
        temp.tab[5]=-a.tab[5]*std::sin(a.tab[0])- pow(a.tab[2], 2)*std::cos(a.tab[0]);
        return temp;
    }
    friend Jet exp(const Jet& a){
        Jet temp;
        temp.tab[0]= std::exp(a.tab[0]);
        temp.tab[1]=a.tab[1]*std::exp(a.tab[0]);
        temp.tab[2]=a.tab[2]*std::exp(a.tab[0]);
        temp.tab[3]=a.tab[3]*std::exp(a.tab[0])+pow(a.tab[1], 2)*std::exp(a.tab[0]);
        temp.tab[4]=a.tab[4]*std::exp(a.tab[0])+a.tab[1]*a.tab[2]*std::exp(a.tab[0]);
        temp.tab[5]=a.tab[5]*std::exp(a.tab[0])+ pow(a.tab[2], 2)*std::exp(a.tab[0]);
        return temp;
    }
   friend std::ostream& operator<<(std::ostream& ost, const Jet& a){
        ost<< a.tab[0]<<" "<<a.tab[1]<<" "<<a.tab[2]<<" "<<a.tab[3]<<" "<<a.tab[4]<<" "<<a.tab[5];
        return ost;
    }
//    operator double()const{
//        return tab[0];
//    }
};

int main() {
    int n=0;
    std::cin>>n;
    double a, b;
    std::cout<<std::setprecision(15);
    for(int i=0; i<n;i++){
        std::cin>>a>>b;
        std::cout<<funkcja(Jet(a, 1),  Jet(b, 2))<<std::endl;
    }
}
