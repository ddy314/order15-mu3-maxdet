// Matrix-level verification, independent of the Fano potential orbit calculation.
#include <algorithm>
#include <array>
#include <iostream>
#include <unordered_map>
using namespace std;
int main(){int a[9][15][15];for(auto&m:a)for(auto&r:m)for(int&v:r)if(!(cin>>v))return 2;
 unordered_map<int,int> target[9];for(int t=0;t<9;t++)for(int j=0;j<15;j++){int code=0;for(int i=0;i<14;i++)code=3*code+a[t][i][j];if(target[t].count(code))return 3;target[t][code]=j;}
 long seen=0;int counts[9][9]{};
 for(int t=0;t<9;t++){array<int,7> p={0,1,2,3,4,5,6};do{for(int mask=0;mask<128;mask++){seen++;int perm[15];perm[14]=14;for(int i=0;i<7;i++){perm[2*i]=2*p[i]+((mask>>i)&1);perm[2*i+1]=2*p[i]+1-((mask>>i)&1);}for(int u=0;u<9;u++){bool ok=true;for(int j=0;j<15;j++){int code=0;for(int i=0;i<14;i++)code=3*code+a[t][perm[i]][j];auto it=target[u].find(code);if(it==target[u].end()){ok=false;break;}for(int i=0;i<15;i++)if(a[t][perm[i]][j]!=a[u][i][it->second])return 4;}if(ok)counts[t][u]++;}}}while(next_permutation(p.begin(),p.end()));}
 cout<<"{\"row_transformations\":"<<seen<<",\"projective_equivalence_counts\":[";for(int i=0;i<9;i++){if(i)cout<<",";cout<<"[";for(int j=0;j<9;j++){if(j)cout<<",";cout<<counts[i][j];}cout<<"]";}cout<<"]}\n";}
