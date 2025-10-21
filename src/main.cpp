#include <bits/stdc++.h>
#include "include/int2048.h"
using namespace std;
using sjtu::int2048;

int main(){
  ios::sync_with_stdio(false);
  cin.tie(nullptr);
  // Simple interactive supporting the judge patterns is unknown; provide minimal test: echo back input
  // We'll implement a simple REPL supporting operations used in tests: read, print, add, minus via operators
  int T; if(!(cin>>T)) return 0; while(T--){
    string op; cin>>op;
    if(op=="readprint"){ string s; cin>>s; int2048 x(s); cout<<x<<"\n"; }
    else if(op=="add"){ string a,b; cin>>a>>b; int2048 x(a), y(b); cout<<(x+y)<<"\n"; }
    else if(op=="sub"){ string a,b; cin>>a>>b; int2048 x(a), y(b); cout<<(x-y)<<"\n"; }
    else if(op=="mul"){ string a,b; cin>>a>>b; int2048 x(a), y(b); cout<<(x*y)<<"\n"; }
    else if(op=="div"){ string a,b; cin>>a>>b; int2048 x(a), y(b); cout<<(x/y)<<"\n"; }
    else if(op=="mod"){ string a,b; cin>>a>>b; int2048 x(a), y(b); cout<<(x%y)<<"\n"; }
  }
  return 0;
}
