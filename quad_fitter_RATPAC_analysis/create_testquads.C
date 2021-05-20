#include <iostream>
#include <vector>
#include <algorithm>
#include <random>

typedef std::vector<int> QUAD;

void print(std::vector<int> const &v)
{
    for (int i: v) {
        std::cout << i << " ";
    }
    std::cout << std::endl;
}
/*
void printQUAD(std::vector<vector<int>> const &q)
{
    for (int i=0; i<q.size();i++) {
        std::cout << i << " "<<q[i];
    }
    std::cout << std::endl;
}
*/
 
int main()
{
  std::vector<int> v = { 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,14, 15, 16 };
 
  std::shuffle(std::begin(v), std::end(v), std::default_random_engine());
 
  print(v);

  int nquad = v.size()/4;
  std::vector<QUAD> vec; vec.reserve(256);
  for (int i = 0; i < nquad; ++i ){
    auto start = v.begin() + 4*i;
    auto stop = start + 4;
    QUAD thequad = QUAD(start, stop);
    cout<<"thequad.size() "<<thequad.size()<<endl;
    cout<<"thequad[i] "<<thequad[i]<<endl;
    vec.push_back(thequad);
  }
  for(int k = 0; k < vec.size(); ++k ){
     for (int i = 0; i < nquad; ++i ){
     cout<<"k: "<<k<<" i "<<i<<" vec[k][i] "<<vec[k][i]<<endl;
  }}
   
  //printQUAD(vec); 

//  std::vector<vector<int>> vec; vec.reserve(256);

  return 0;

 
 
 
}
