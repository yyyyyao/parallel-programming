#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>

using namespace std;
int main(int argc, char** argv)
{
  double cu = 1.0;
  double cv = 1.0;
  double long sum = 0.0;

  int before[8][8] ={
  {57,49,44,39,33,28,23,22},
  {55,45,39,33,28,21,21,19},
  {55,44,37,31,20,17,17,16},
  {58,44,37,29,17,15,16,14},
  {66,46,31,22,16,14,12,10},
  {71,49,32,24,14,12,21,19},
  {69,50,30,22,12,4,-1,-2},
  {62,42,25,17,7,1,-16,-16}

  };
  int after[8][8];
  memset(after, 0, sizeof(int)*64);

  for(int v =0; v < 8; v++)
  {
    if(v) cv =1.0;
    if(!v) cv = 1 / sqrt(2.0);
    for(int u = 0; u < 8; u++)
    {
      if(u) cu =1.0;
      if(!u) cu = 1 / sqrt(2.0);
      sum = 0.0;
      for(int y = 0; y < 8; y++)
      {
        for(int x = 0; x < 8; x++)
        {
          double long a = cos((2 * x + 1)*u*M_PI/16) * cos(( 2 * y + 1)*v*M_PI/16);
          double long d = before[y][x];
          sum += d * a;
        }
      }
      after[v][u] = int(floor(cu*cv*sum/4 + 0.5));
    }
  }

  for(int j = 0; j < 8; j++)
  {
    for(int i = 0; i < 8; i++)
    {
      cout << after[j][i] << ",";
    }
  cout << endl;
  }
    return 0;
}

