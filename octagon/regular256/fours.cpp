#include <cstdio>
#include <vector>
#include <array>
#include <algorithm>
using namespace std;
vector<int> read(const char*s){FILE*f=fopen(s,"rb");fseek(f,0,SEEK_END);vector<int>v(ftell(f)/4);rewind(f);fread(v.data(),4,v.size(),f);fclose(f);return v;}
int main(int argc,char**argv){auto g=read("good.bin"),ids=read("ids.bin");vector<int>ix(1<<25,-1);vector<char>seen(20000);int n=g.size();
for(int i=0;i<n;i++)ix[g[i]]=i;vector<array<signed char,640>>L(n);signed char B[640];FILE*f=fopen("loads.bin","rb");fread(L.data(),640,n,f);fclose(f);f=fopen("B.bin","rb");fread(B,1,640,f);fclose(f);FILE*out=fopen("four.txt","w");
for(int i=0;i<n;i++){int a=g[i],h=ids[a];if(seen[h])continue;seen[h]=1;for(int j=0;j<n;j++){int b=g[j],c=a^b;if(b>=c||ix[c]<0||ids[b]<h||ids[c]<h)continue;bool ok=true;
for(int t=0;t<640;t++){int x=L[i][t],y=L[j][t],z=L[ix[c]][t];if(B[t]+max({x+y+z,x-y-z,-x+y-z,-x-y+z})>16){ok=false;break;}}
if(ok)fprintf(out,"0 %d %d\n",a,b);}}
for(int j=1;j<argc;j++){f=fopen(argv[j],"r");int a,n,v,b[25],s;fscanf(f,"%d%d%d",&a,&n,&v);for(int i=0;i<n;i++)fscanf(f,"%d",&b[i]);vector<array<int,9>>rows;
while(fscanf(f,"%d",&s)==1){array<int,9>r{};r[0]=s;for(int i=1;i<=s;i++)fscanf(f,"%d",&r[i]);rows.push_back(r);}fclose(f);
sort(rows.begin(),rows.end(),[&](auto&p,auto&q){auto sc=[&](auto&r){int z=0;for(int i=1;i<=r[0];i++)z|=r[i];return __builtin_popcount(z&((1<<n)-1));};return sc(p)<sc(q);});
for(int k=0;k<1<<n;k++){int m=k|(1<<n);bool ok=true;for(auto &r:rows){int c=0;for(int i=1;i<=r[0];i++)c+=__builtin_parity(m&r[i]);if(c>4||r[0]-c>4){ok=false;break;}}if(ok){int x=v;for(int i=0;i<n;i++)if(k>>i&1)x^=b[i];fprintf(out,"1 %d %d\n",a,x);}}}fclose(out);}
