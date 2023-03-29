#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

//#include "global-p.h"
//#include "struct-p.h"

static unsigned short gf[ORD], fg[ORD];

// nomal bases
// unsigned short gf[M]={0,1,2,4,8,9,11,15,7,14,5,10,13,3,6,12};
// unsigned short fg[M]={0,1,2,13,3,10,14,8,4,5,11,6,15,12,9,7};

// sage比較用
// unsigned short gf[16]={0,1,2,4,8,3,6,12,11,5,10,7,14,15,13,9};
// unsigned short fg[16]={0,1,2,5,3,9,6,11,4,15,10,8,7,14,12,13};

// unsigned short gf[M]={0,1,2,4,8,16,32,64,128,29,58,116,232,205,135,19,38,76,152,45,90,180,117,234,201,143,3,6,12,24,48,96,192,157,39,78,156,37,74,148,53,106,212,181,119,238,193,159,35,70,140,5,10,20,40,80,160,93,186,105,210,185,111,222,161,95,190,97,194,153,47,94,188,101,202,137,15,30,60,120,240,253,231,211,187,107,214,177,127,254,225,223,163,91,182,113,226,217,175,67,134,17,34,68,136,13,26,52,104,208,189,103,206,129,31,62,124,248,237,199,147,59,118,236,197,151,51,102,204,133,23,46,92,184,109,218,169,79,158,33,66,132,21,42,84,168,77,154,41,82,164,85,170,73,146,57,114,228,213,183,115,230,209,191,99,198,145,63,126,252,229,215,179,123,246,241,255,227,219,171,75,150,49,98,196,149,55,110,220,165,87,174,65,130,25,50,100,200,141,7,14,28,56,112,224,221,167,83,166,81,162,89,178,121,242,249,239,195,155,43,86,172,69,138,9,18,36,72,144,61,122,244,245,247,243,251,235,203,139,11,22,44,88,176,125,250,233,207,131,27,54,108,216,173,71,142};
// unsigned short fg[M]={0,1,2,26,3,51,27,199,4,224,52,239,28,105,200,76,5,101,225,15,53,142,240,130,29,194,106,249,201,9,77,114,6,139,102,48,226,37,16,34,54,148,143,219,241,19,131,70,30,182,195,126,107,40,250,186,202,155,10,121,78,229,115,167,7,192,140,99,103,222,49,254,227,153,38,180,17,146,35,137,55,209,149,207,144,151,220,190,242,211,20,93,132,57,71,65,31,67,183,164,196,73,127,111,108,59,41,85,251,134,187,62,203,95,156,160,11,22,122,44,79,213,230,173,116,244,168,88,8,113,193,248,141,129,100,14,104,75,223,238,50,198,255,25,228,166,154,120,39,185,181,125,18,69,147,218,36,33,138,47,56,64,210,92,150,189,208,206,145,136,152,179,221,253,191,98,243,87,212,172,21,43,94,159,133,61,58,84,72,110,66,163,32,46,68,217,184,124,165,119,197,24,74,237,128,13,112,247,109,162,60,83,42,158,86,171,252,97,135,178,188,205,63,91,204,90,96,177,157,170,161,82,12,246,23,236,123,118,45,216,80,175,214,234,231,232,174,233,117,215,245,235,169,81,89,176};

// unsigned short gf[27]={0,1,3,9,5,15,23,13,17,20,4,12,14,11,2,6,18,7,21,16,26,22,10,8,24,25,19,};
// unsigned short fg[27]={0,1,14,2,10,4,15,17,23,3,22,13,11,7,12,5,19,8,16,26,9,18,21,6,24,25,20,};
unsigned char tmp[E * K];
unsigned char pub[E * K];
unsigned char BH[E * K];
short mat[K * E][8192];
short ma[K][8192];
unsigned short P[M]={0};
unsigned short inv_P[M]={0};
//short uu;
short bm[K * E][8192];
short bm2[K * E][8192];
unsigned char aa[64]; //={ 148, 246, 52, 251, 16, 194, 72, 150, 249, 23, 90, 107, 151, 42, 154, 124, 48, 58, 30, 24, 42, 33, 38, 10, 115, 41, 164, 16, 33, 32, 252, 143, 86, 175, 8, 132, 103, 231, 95, 190, 61, 29, 215, 75, 251, 248, 72, 48, 224, 200, 147, 93, 112, 25, 227, 223, 206, 137, 51, 88, 109, 214, 17, 172};

#define I8T char
#define U8C(v) (v##U)

#define U8V(v) ((unsigned char)(v)&U8C(0xFF))
#define ROTL8(v, n) \
  (U8V((v) << (n)) | ((v) >> (8 - (n))))

#define R(x, n) (((x) << (n)) | ((x) >> (32 - (n))))

unsigned int rotate_left(unsigned int x, int n)
{
  assert(0 < n && n < 32);
  return (x << n) | (x >> (32 - n));
}

#define str_length 128
#define password_length 256

char password[password_length + 1];
#define SIZE_OF_ARRAY(array) (sizeof(array) / sizeof(array[0]))
#define SWAP(type, a, b) \
  {                      \
    type work = a;       \
    a = b;               \
    b = work;            \
  }

/*
    Fisher-Yates shuffle による方法
    配列の要素をランダムシャッフルする
*/
void random_shuffle(unsigned short *array, size_t size)
{
  for (size_t i = size; i > 1; --i)
  {
    size_t a = i - 1;
    size_t b = rand() % i;
    SWAP(int, array[a], array[b]);
  }
}

unsigned long xor128(void)
{
  unsigned int a = 0;

  static unsigned long x = 123456789, y = 362436069, z = 521288629, w = 88675123;
  unsigned long t;

  a = rand();
  t = x ^ (a << 11);
  a = y;
  y = z;
  z = w;
  return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

void seed(void)
{
  /*
   * 変数宣言
   */
  char str[str_length + 1];
  time_t t;
  int i, j, k, rnd;

  /*
   * 乱数の初期化
   */
  srand(clock() + time(&t));

  /*
   * 乱数生成とパスワードの生成
   */
  for (i = 0; i < str_length; i++)
  {
    for (j = 0; j < 2; j++)
    {
      k = i * 2 + j;
      do
      {
        rnd = rand();
        password[k] = str[i] + rnd;
      } while (!isalnum(password[k]));
    }
  }

  /*
   * NULL文字の挿入
   */
  password[password_length] = '\0';

  /*
   * パスワードの出力
   */
  //    printf("生成パスワード：%s",password);

  return;
}

/*
int mlt(int x, int y){

    if(x==0||y==0)
        return 0;

  return (x*y)%P;
}
*/

int mlt(int x, int y)
{

  if (x == 0 || y == 0)
    return 0;

  return ((x + y - 2) % (ORD - 1)) + 1;
}

int mltn(int n, int x)
{
  int i, j;

  if (n == 0)
    return 1;
  i = x;
  for (j = 0; j < n - 1; j++)
    i = mlt(i, x);

  return i;
}

int mltnc(int n, int x)
{
  int ret = 1;
  while (n > 0)
  {
    if (n & 1)
      ret = mlt(ret, x); // n の最下位bitが 1 ならば x^(2^i) をかける
    x = mlt(x, x);
    n >>= 1; // n を1bit 左にずらす
  }
  return ret;
}

int mltn2(int n, int x)
{
  int i, j;

  if (n == 0)
    return 1;
  i = x;
  for (j = 0; j < n - 1; j++)
    i = i * x % Pr;

  return i;
}

unsigned int mltn0(unsigned int n, unsigned int u)
{
  if (n % N == 0)
    return 1;
  if (u == 0)
    return 0;
  return (u * n - n) % (N - 1) + 1;
}

