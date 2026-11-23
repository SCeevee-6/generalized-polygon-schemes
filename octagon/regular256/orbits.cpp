#include <vector>
#include <cstdio>
using namespace std;
int tab[3][4][256],cols[3][25];
int main(){FILE*f=fopen("act.bin","rb");fread(cols,4,75,f);fclose(f);
for(int h=0;h<3;h++)for(int i=0;i<25;i++)for(int j=1<<(i%8);j<2<<(i%8);j++)tab[h][i/8][j]=cols[h][i]^tab[h][i/8][j^(1<<(i%8))];
vector<int> ids(1<<25,-1),q;int n=0;FILE*r=fopen("reps.txt","w");
for(int v=0;v<1<<25;v++)if(ids[v]<0){q={v};ids[v]=n;for(unsigned k=0;k<q.size();k++){int x=q[k];for(auto &t:tab){int y=t[0][x&255]^t[1][x>>8&255]^t[2][x>>16&255]^t[3][x>>24];if(ids[y]<0){ids[y]=n;q.push_back(y);}}}fprintf(r,"%d %d\n",v,(int)q.size());n++;}fclose(r);f=fopen("ids.bin","wb");fwrite(ids.data(),4,ids.size(),f);fclose(f);}
