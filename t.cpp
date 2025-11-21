 #include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
#include <algorithm>
using namespace std;       
int main(){        
        vector<int> times ;
        vector<vector<double>>orderpara(10,vector<double>(5,0));
        for(int i=0;i<10;i++){
            for(int j=0;j<5;j++)
            {   if(j==0)times.push_back(i);
                orderpara[i][j]=(j+10);
            }    
        }
        ofstream order_file("order_parameter.csv");
        string a="";
        for (int i=0;i<orderpara[0].size();i++)a+=to_string(i)+",";      
        order_file<<a<<",t\n";
        for(int i=0;i<orderpara.size();i++){
            for(int j=0;j<orderpara[0].size();j++){
                order_file<<orderpara[i][j]<<",";
            }
            order_file<<times[i]<<"\n";
        }
        order_file.close();
}