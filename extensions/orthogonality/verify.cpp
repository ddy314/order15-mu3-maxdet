// Exact balanced-row graph generation and maximum clique replay.
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <vector>
using namespace std;
using Bits=vector<uint64_t>;
struct Row {uint16_t a=0,b=0;array<int,15> e{};};
Row row(array<int,15> e){Row r;r.e=e;for(int i=0;i<15;i++){if(e[i]==1)r.a|=1<<i;if(e[i]==2)r.b|=1<<i;}return r;}
bool orth(const Row&x,const Row&y){uint16_t x0=32767^(x.a|x.b),y0=32767^(y.a|y.b);return __builtin_popcount((x0&y0)|(x.a&y.a)|(x.b&y.b))==5 && __builtin_popcount((x.a&y0)|(x.b&y.a)|(x0&y.b))==5;}
bool any(const Bits&x){for(auto a:x)if(a)return true;return false;}
int first(const Bits&x){for(int i=0;i<(int)x.size();i++)if(x[i])return 64*i+__builtin_ctzll(x[i]);return -1;}
void erase(Bits&x,int i){x[i/64]&=~(1ULL<<(i%64));}
struct Clique {vector<Bits> adj;vector<int> chosen,best;uint64_t nodes=0;
 void expand(Bits p){nodes++;vector<int> vs,cs;Bits left=p;int color=0;
  while(any(left)){color++;Bits q=left;while(any(q)){int v=first(q);erase(left,v);erase(q,v);vs.push_back(v);cs.push_back(color);for(int j=0;j<(int)q.size();j++)q[j]&=~adj[v][j];}}
  for(int i=(int)vs.size()-1;i>=0;i--){if(chosen.size()+cs[i]<=best.size())return;int v=vs[i];chosen.push_back(v);Bits q=p;for(int j=0;j<(int)q.size();j++)q[j]&=adj[v][j];if(any(q))expand(q);else if(chosen.size()>best.size())best=chosen;chosen.pop_back();erase(p,v);}}
 vector<int> solve(){int n=adj.size();Bits p((n+63)/64,~0ULL);if(n%64)p.back()=(1ULL<<(n%64))-1;expand(p);return best;}
};
void check_small(){for(int n=0;n<=6;n++){int m=n*(n-1)/2;for(int mask=0;mask<(1<<m);mask++){Clique c;c.adj.assign(n,Bits((n+63)/64));int k=0;for(int i=0;i<n;i++)for(int j=0;j<i;j++,k++)if(mask>>k&1)c.adj[i][0]|=1ULL<<j,c.adj[j][0]|=1ULL<<i;int best=0;for(int s=0;s<(1<<n);s++){bool ok=true;for(int i=0;i<n;i++)if(s>>i&1)for(int j=0;j<i;j++)if((s>>j&1)&&!(c.adj[i][0]>>j&1))ok=false;if(ok)best=max(best,__builtin_popcount((unsigned)s));}if((int)c.solve().size()!=best)throw runtime_error("small-graph negative control");}}}
int main(int argc,char**argv){if(argc!=2)throw runtime_error("output JSON required");check_small();array<int,15> ee{};for(int i=0;i<15;i++)ee[i]=i/5;Row b=row(ee);vector<Row> common;int limit=4782969;
 for(int code=0;code<limit;code++){int t=code;array<int,15> e{};int cnt[3]={1,0,0};for(int i=1;i<15;i++){e[i]=t%3;t/=3;cnt[e[i]]++;}if(cnt[0]!=5||cnt[1]!=5)continue;Row x=row(e);if(orth(x,b))common.push_back(x);}
 if(common.size()!=37251)throw runtime_error("common count");
 array<array<int,3>,5> types={{{5,0,0},{4,1,0},{3,2,0},{3,1,1},{2,2,1}}};int expectedN[]={0,5184,5022,5130,5112};long expectedE[]={0,1741824,1604286,1655802,1667484};int expectedM[]={0,5,6,6,6};
 ofstream out(argv[1]);out<<"{\"small_graphs_checked\":33868,\"common_vertices\":37251,\"cases\":[";
 for(int typ=0;typ<5;typ++){auto a=types[typ];array<int,15> e{};int pos=0;int counts[3][3]={{a[0],a[1],a[2]},{a[1],a[2],a[0]},{a[2],a[0],a[1]}};for(int i=0;i<3;i++)for(int j=0;j<3;j++)for(int k=0;k<counts[i][j];k++)e[pos++]=j;Row c=row(e);vector<Row> vertices;for(auto& x:common)if(orth(x,c))vertices.push_back(x);int n=vertices.size();if(n!=expectedN[typ])throw runtime_error("candidate count");Clique search;search.adj.assign(n,Bits((n+63)/64));long edges=0;
 for(int i=0;i<n;i++)for(int j=0;j<i;j++)if(orth(vertices[i],vertices[j])){search.adj[i][j/64]|=1ULL<<(j%64);search.adj[j][i/64]|=1ULL<<(i%64);edges++;}auto best=search.solve();if(edges!=expectedE[typ]||(int)best.size()!=expectedM[typ])throw runtime_error("clique mismatch");
 if(typ)out<<",";out<<"{\"vertices\":"<<n<<",\"edges\":"<<edges<<",\"maximum\":"<<best.size()<<",\"nodes\":"<<search.nodes<<",\"witness\":[";vector<Row> witness={row({}),b,c};for(int v:best)witness.push_back(vertices[v]);for(int i=0;i<(int)witness.size();i++){if(i)out<<",";out<<"[";for(int j=0;j<15;j++){if(j)out<<",";out<<witness[i].e[j];}out<<"]";}out<<"]}";cout<<"type "<<typ<<": "<<n<<" vertices, "<<edges<<" edges, maximum "<<best.size()<<", nodes "<<search.nodes<<endl;
 }out<<"]}\n";
}
