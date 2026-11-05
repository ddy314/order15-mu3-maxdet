// Independent exact checker: scaled inverse identity and incremental quadratic scan.
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;
struct Z {long long a=0,b=0;Z operator+(Z y)const{return{a+y.a,b+y.b};}Z operator-(Z y)const{return{a-y.a,b-y.b};}Z operator*(Z y)const{return{a*y.a-b*y.b,a*y.b+b*y.a-b*y.b};}Z bar()const{return{a-b,-b};}Z times(long long t)const{return{a*t,b*t};}bool operator==(Z y)const{return a==y.a&&b==y.b;}bool zero()const{return a==0&&b==0;}long long norm()const{return a*a-a*b+b*b;}};
Z roots[3]={{1,0},{0,1},{-1,-1}},units[4]={{1,0},{-1,0},{0,1},{0,-1}};
int main(int argc,char**argv){if(argc!=2)throw runtime_error("text certificate path required");ifstream in(argv[1]);int ns[4];for(int&i:ns)in>>i;Z eps[7];int pos=0,K=0;for(int t=0;t<4;t++)for(int j=0;j<ns[t];j++){eps[pos++]=units[t];K+=2*units[t].a-units[t].b;}int M=290+K;Z k[15],r[15][15]{},g[15][15]{};k[14]={-24,24};for(int i=0;i<7;i++){k[2*i]=Z{5,0}-eps[i].bar();k[2*i+1]=Z{5,0}-eps[i];g[2*i][2*i+1]=eps[i].times(3);g[2*i+1][2*i]=eps[i].bar().times(3);r[2*i][2*i]=r[2*i+1][2*i+1]={5*M,0};r[2*i][2*i+1]=eps[i].times(-M);r[2*i+1][2*i]=eps[i].bar().times(-M);}for(int i=0;i<15;i++){g[i][i]={15,0};if(i<14)g[i][14]={1,-1},g[14][i]={2,1};for(int j=0;j<15;j++)r[i][j]=r[i][j]+k[i].bar()*k[j];}
 for(int i=0;i<15;i++)for(int j=0;j<15;j++){Z s;for(int q=0;q<15;q++)s=s+g[i][q]*r[q][j];if(!(s==Z{i==j?72*M:0,0}))throw runtime_error("inverse identity");}
 int n;in>>n;vector<unsigned> codes(n);vector<int> colors(n);for(auto&c:codes)in>>c;for(auto&c:colors)in>>c;vector<array<Z,15>> cols(n),inv(n);vector<int> lookup(4782969,-1);for(int q=0;q<n;q++){lookup[codes[q]]=q;unsigned t=codes[q];for(int j=0;j<7;j++){int d=t%9;t/=9;cols[q][2*j]=roots[d/3];cols[q][2*j+1]=roots[d%3];}cols[q][14]=roots[0];for(int i=0;i<15;i++)for(int j=0;j<15;j++)inv[q][i]=inv[q][i]+r[i][j]*cols[q][j];}
 array<int,14> digits{};array<Z,15> c,rc{};for(auto&z:c)z=roots[0];for(int i=0;i<15;i++)for(int j=0;j<15;j++)rc[i]=rc[i]+r[i][j];long long q=0;for(auto z:rc)q+=z.a; // sum is real
 int found=0;for(int serial=0;serial<4782969;serial++){if(q==72*M){unsigned code=0,p=1;for(int j=0;j<7;j++){code+=(3*digits[2*j]+digits[2*j+1])*p;p*=9;}if(lookup[code]<0)throw runtime_error("omitted candidate column");found++;}for(int i=0;i<14;i++){int next=(digits[i]+1)%3;Z d=roots[next]-c[i];Z cross=d.bar()*rc[i];q+=2*cross.a-cross.b+d.norm()*r[i][i].a;for(int j=0;j<15;j++)rc[j]=rc[j]+r[j][i]*d;c[i]=roots[next];digits[i]=next;if(next)break;}}
 if(found!=n)throw runtime_error("extra column");long long pairs=0;auto edge=[&](int a,int b){Z s;for(int j=0;j<15;j++)s=s+cols[a][j].bar()*inv[b][j];return s.zero();};for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(colors[i]==colors[j]){pairs++;if(edge(i,j))throw runtime_error("improper global coloring");}
 vector<int> covered(n);int locals;in>>locals;for(int z=0;z<locals;z++){int v,len;in>>v>>len;covered[v]=1;vector<int> neigh(len),co(len);for(auto&x:neigh)in>>x;for(auto&x:co)in>>x;vector<int> seen(n);for(int i=0;i<len;i++){if(neigh[i]<0||neigh[i]>=n||seen[neigh[i]]++)throw runtime_error("duplicate neighbor");if(co[i]<0||co[i]>=13)throw runtime_error("local color bound");}for(int u=0;u<n;u++)if(u!=v && edge(v,u)!=(bool)seen[u])throw runtime_error("omitted actual neighbor");for(int i=0;i<len;i++)for(int j=0;j<i;j++)if(co[i]==co[j]&&edge(neigh[i],neigh[j]))throw runtime_error("improper local coloring");}
 for(int i=0;i<n;i++)if(colors[i]<0||(colors[i]>=14&&!covered[i]))throw runtime_error("uncovered high color");cout<<found<<" "<<pairs<<" "<<locals<<endl;
}
