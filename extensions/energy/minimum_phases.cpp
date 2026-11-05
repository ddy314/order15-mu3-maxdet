// Complete integer column scans and explicit proper-coloring certificates.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <vector>
using namespace std;
struct Z {int a=0,b=0;Z operator+(Z y)const{return{a+y.a,b+y.b};}Z operator-(Z y)const{return{a-y.a,b-y.b};}Z operator*(Z y)const{return{a*y.a-b*y.b,a*y.b+b*y.a-b*y.b};}Z times(int n)const{return{a*n,b*n};}Z bar()const{return{a-b,-b};}int norm()const{return a*a-a*b+b*b;}bool zero()const{return a==0&&b==0;}};
Z roots[3]={{1,0},{0,1},{-1,-1}}, units[4]={{1,0},{-1,0},{0,1},{0,-1}};
struct Col {uint32_t code;array<uint8_t,14> e;Z x;array<Z,14> inv;};
bool norm_ok(int n){for(int p=2;p*p<=n;p++){int k=0;while(n%p==0)n/=p,k++;if(p%3==2&&k%2)return false;}return n%3!=2;}
Z bilinear(const Col&c,const Col&d,int M){Z z;for(int i=0;i<14;i++)z=z+c.inv[i]*roots[d.e[i]];return z.times(M)+c.x.bar()*d.x;}
vector<int> color(const vector<vector<int>>&adj,const vector<int>&subset){int n=adj.size();vector<int> deg(n,-1),removed(n),order;priority_queue<pair<int,int>,vector<pair<int,int>>,greater<pair<int,int>>> pq;for(int v:subset)deg[v]=0;for(int v:subset){for(int u:adj[v])if(deg[u]>=0)deg[v]++;pq.push({deg[v],v});}while(!pq.empty()){auto[d,v]=pq.top();pq.pop();if(removed[v]||deg[v]!=d)continue;removed[v]=1;order.push_back(v);for(int u:adj[v])if(deg[u]>=0&&!removed[u])pq.push({--deg[u],u});}vector<int> co(n,-1);for(auto it=order.rbegin();it!=order.rend();it++){int v=*it;uint64_t used=0;for(int u:adj[v])if(co[u]>=0){if(co[u]>=63)throw runtime_error("color overflow");used|=1ULL<<co[u];}co[v]=__builtin_ctzll(~used);}return co;}
int main(int argc,char**argv){if(argc!=2)throw runtime_error("output directory required");string directory=argv[1];int kept=0;vector<int> counts;for(int n0=0;n0<=7;n0++)for(int n1=0;n1<=7-n0;n1++)for(int n2=0;n2<=7-n0-n1;n2++){int ns[4]={n0,n1,n2,7-n0-n1-n2};array<Z,7> eps;int j=0,K=0;for(int a=0;a<4;a++)for(int b=0;b<ns[a];b++){eps[j++]=units[a];K+=2*units[a].a-units[a].b;}int M=290+K;if(!norm_ok(M))continue;kept++;
 Z ls[7][9];int bs[7][9];array<array<Z,2>,7> unused;for(j=0;j<7;j++)for(int x=0;x<3;x++)for(int y=0;y<3;y++){ls[j][3*x+y]=(Z{5,0}-eps[j].bar())*roots[x]+(Z{5,0}-eps[j])*roots[y];Z b=roots[x].bar()*eps[j]*roots[y];bs[j][3*x+y]=2*b.a-b.b;}
 vector<Col> cols;array<int,7> digits{};Z L;int b=0;for(j=0;j<7;j++){L=L+ls[j][0];b+=bs[j][0];}const int LIMIT=4782969;
 for(int code=0;code<LIMIT;code++){Z X=Z{24,-24}-L;if(X.norm()==M*(2+b)){Col c;c.code=code;c.x=X;for(j=0;j<7;j++){int x=digits[j]/3,y=digits[j]%3;c.e[2*j]=x;c.e[2*j+1]=y;c.inv[2*j]=(roots[x].bar().times(5)-roots[y].bar()*eps[j].bar());c.inv[2*j+1]=(roots[y].bar().times(5)-roots[x].bar()*eps[j]);}cols.push_back(c);}for(j=0;j<7;j++){int old=digits[j],nw=(old+1)%9;digits[j]=nw;L=L-ls[j][old]+ls[j][nw];b+=bs[j][nw]-bs[j][old];if(nw)break;}}
 int n=cols.size();if(K==14){if(n!=77260)throw runtime_error("friendship count");cout<<"friendship "<<n<<endl;continue;}vector<vector<int>> adj(n);long edges=0;
 for(int i=0;i<n;i++)for(int k=0;k<i;k++)if(bilinear(cols[i],cols[k],M).zero()){adj[i].push_back(k);adj[k].push_back(i);edges++;}
 vector<int> all(n);iota(all.begin(),all.end(),0);auto co=color(adj,all);int colors=co.empty()?0:1+*max_element(co.begin(),co.end());vector<pair<int,vector<int>>> locals;
 if(colors>14){for(int v:all)if(co[v]>=14){auto cc=color(adj,adj[v]);int mx=0;for(int u:adj[v])mx=max(mx,cc[u]+1);if(mx>13)throw runtime_error("local coloring insufficient");locals.push_back({v,cc});}}
 // Independently validate every produced certificate against the exact form.
 for(int i=0;i<n;i++)for(int k=0;k<i;k++)if(co[i]==co[k]&&bilinear(cols[i],cols[k],M).zero())throw runtime_error("invalid coloring");;
 string path=directory+"/phase_"+to_string(ns[0])+"_"+to_string(ns[1])+"_"+to_string(ns[2])+".json";ofstream out(path);out<<"{\"multiplicities\":["<<ns[0]<<","<<ns[1]<<","<<ns[2]<<","<<ns[3]<<"],\"K\":"<<K<<",\"edges\":"<<edges<<",\"columns\":[";for(int i=0;i<n;i++){if(i)out<<",";out<<cols[i].code;}out<<"],\"coloring\":[";for(int i=0;i<n;i++){if(i)out<<",";out<<co[i];}out<<"],\"local_colorings\":[";for(int i=0;i<(int)locals.size();i++){if(i)out<<",";int v=locals[i].first;out<<"{\"vertex\":"<<v<<",\"neighbors\":[";for(int k=0;k<(int)adj[v].size();k++){if(k)out<<",";out<<adj[v][k];}out<<"],\"colors\":[";for(int k=0;k<(int)adj[v].size();k++){if(k)out<<",";out<<locals[i].second[adj[v][k]];}out<<"]}";}out<<"]}\n";cout<<ns[0]<<","<<ns[1]<<","<<ns[2]<<","<<ns[3]<<": "<<n<<" columns "<<edges<<" edges "<<colors<<" colors, local "<<locals.size()<<endl;
 }if(kept!=35)throw runtime_error("norm survivors");}
