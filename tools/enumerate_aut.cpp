// Exhaustive projective monomial stabilizer of the explicitly switched witness.
// Input is 15x15 ternary, last row zero, canonical Gram as checked in Python.
#include <algorithm>
#include <array>
#include <iostream>
#include <unordered_map>
#include <vector>
using namespace std;

int main(){
    int A[15][15];
    for(auto &r:A)for(int &x:r)if(!(cin>>x))return 2;
    unordered_map<int,int> codes;
    for(int j=0;j<15;j++){
        int z=0;
        for(int i=0;i<14;i++)z=3*z+A[i][j];
        if(codes.count(z))return 3;
        codes[z]=j;
    }
    array<int,7> pair={0,1,2,3,4,5,6};
    int count=0,seen=0;
    cout<<"[\n";
    do{
        for(int mask=0;mask<128;mask++){
            seen++;
            array<int,15> p,c;
            p[14]=14;
            for(int i=0;i<7;i++){
                p[2*i]=2*pair[i]+((mask>>i)&1);
                p[2*i+1]=2*pair[i]+(1-((mask>>i)&1));
            }
            bool ok=true;
            for(int j=0;j<15;j++){
                int z=0;
                for(int i=0;i<14;i++)z=3*z+A[p[i]][j];
                auto it=codes.find(z);
                if(it==codes.end()){ok=false;break;}
                c[j]=it->second;
            }
            if(ok){
                if(count)cout<<",\n";
                cout<<"{\"rows\":[";
                for(int i=0;i<15;i++){if(i)cout<<',';cout<<p[i];}
                cout<<"],\"columns\":[";
                for(int i=0;i<15;i++){if(i)cout<<',';cout<<c[i];}
                cout<<"]}";
                count++;
            }
        }
    }while(next_permutation(pair.begin(),pair.end()));
    cout<<"\n]\n";
    cerr<<"searched="<<seen<<" projective_order="<<count<<" monomial_order="<<3*count<<'\n';
}
