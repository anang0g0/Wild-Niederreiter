#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

// #define Pr 13
#define F K *E

short oinb(short o)
{
    if (ORD == 0)
        return 0;
    if (ORD == 1)
        return 1;
    for (int i = 0; i < Pr; i++)
    {
        for (int j = 0; j < Pr; j++)
        {
            if ((o * i) % Pr == 1)
                return i;
        }
    }
}

// short fu(short d)
//{
//     return (Pr - d) % Pr;
// }

int genS(CTX a, CTX *inv_a)
{
    // short a.x[F][F] = {{0,2},{2,2}}; //{{2,1,1,1},{1,1,2,2},{2,0,1,1},{1,2,1,1}}; //入力用の配列
    static CTX b = {0};
    static unsigned c[F][F] = {0};
    // static CTX inv_a; // ここに逆行列が入る
    short buf, A = 0; // 一時的なデータを蓄える

    srand(clock());
label:
    b = a;
    // 単位行列を作る
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            inv_a->x[i][j] = (i == j) ? 1 : 0;
        }
    }
    // 掃き出し法
    for (int i = 0; i < F; i++)
    {
        printf("n=%d\n", i);
        if (a.x[i][i] == 0)
        {
            // a.x[i][i]=1;
            // printf("aa.x[%d][%d]=%d\n", i, i, a.x[i][i]);
            int s = i;
            // short tmp[F] = {0}, tt[F] = {0};
            while (a.x[s][i] == 0)
            {
                s++;
                if (s == F)
                {
                    // a.b = 0;
                    printf("線形従属\n");
                    exit(1);
                    return 0;
                }
            }
            // then 0
            for (int u = 0; u < F; u++)
            {
                a.x[i][u] = (a.x[i][u] + a.x[s][u]) % Pr;
                inv_a->x[i][u] = (inv_a->x[i][u] + inv_a->x[s][u]) % Pr;
            }
        }
        A = inv_a->x[0][0];
        // i++;
        // int s = 0;

        buf = oinb(a.x[i][i]);
        if (buf == 0)
        {
            printf("buf baka\n");
            exit(1);
        }
        // printf("buf=%d\n", buf);
        //  for(k=0;k<F;k++)
        if (buf > 0)
        {

            for (int j = 0; j < F; j++)
            {
                a.x[i][j] = a.x[i][j] * buf % Pr;
                inv_a->x[i][j] = inv_a->x[i][j] * buf % Pr;
                // printf("%d %d %d\n",i,j,a.x[i][j]);
            }
            for (int j = 0; j < F; j++)
            {
                // printf("%d:%d ",i,a.x[0][j]);
                if (i != j)
                {
                    buf = a.x[j][i];
                    for (int k = 0; k < F; k++)
                    {
                        a.x[j][k] = (a.x[j][k] - (a.x[i][k] * buf) % Pr);
                        if (a.x[j][k] < 0)
                            a.x[j][k] += Pr;
                        // if(inv_a.x[i][k]!=0)
                        inv_a->x[j][k] = (inv_a->x[j][k] - (inv_a->x[i][k] * buf) % Pr);
                        if (inv_a->x[j][k] < 0)
                            inv_a->x[j][k] += Pr;
                    }

                    if (a.x[i][i] != 1)
                    {
                        printf("%d:%d ", i, a.x[i][i]);
                        exit(1);
                    }
                    // printf("\n");
                }
            }
        }
    }

    printf("検算\n");
    unsigned tmp = 0;
    unsigned z[F][F] = {0};
    int count = 0;
    int ca = 0;

    // 検算

    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            tmp = 0;
            for (int k = 0; k < F; k++)
                tmp += ((b.x[i][k]) * (inv_a->x[k][j])) % Pr;
            z[i][j] = tmp % Pr;
            // printf("%d ", z[i][j]);
            if (z[i][j] == 0 && i != j)
            {
                ca++;
            }
        }
        printf("\n");
        if (z[i][i] == 1)
        {
            count++;
        }
        else if (z[i][i] != 1)
        {
            // inv_a.b = 0;
            printf("baka\n");
            exit(1);
            // return inv_a;
        }
    }
    if (ca == F * F - F && count == F)
    {
        // inv_a.b = 1;
        return 1; // inv_a;
    }
    else
    {
        for (int i = 0; i < F; i++)
        {
            for (int j = 0; j < F; j++)
                printf("%s,", b.x[i][j]);
            printf("\n");
        }
        printf("\n");
        exit(1);
        // inv_a.b = 0;
        return 0; // inv_a;
    }

    printf("もとの行列\n");
    exit(1);
    count = 0;
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            printf("%d ", b.x[i][j] % Pr);
        }
        printf("\n");
    }
    printf("\n");
    // int count = 0;
    for (int i = 0; i < F; i++)
    {
        for (int j = 0; j < F; j++)
        {
            printf("%d", z[i][j] % Pr);
        }
        printf("\n");
        if (z[i][i] % Pr == 1)
            count++;
    }

    if (count != F)
    {
        // inv_a.b = 0;
        // exit(1);
        return 0; // inv_a;
    }

    // exit(1);
    // inv_a.b = 1;
    return 0; // inv_a;
}

void kenzan(CTX a, CTX *inv_a, int p)
{
    unsigned tmp = 0;
    static CTX z = {0};
    // 検算
    if (p == 1)
    {
        for (int i = 0; i < F; i++)
        {
            for (int j = 0; j < F; j++)
            {
                tmp = 0;
                for (int k = 0; k < F; k++)
                    tmp += ((a.x[i][k] % Pr) * (inv_a->x[k][j] % Pr)) % Pr;
                fls.x[i][j] = tmp % Pr;
            }
        }
        return; // z;
    }

    if (p == 2)
    {
        for (int i = 0; i < F; i++)
        {
            for (int j = 0; j < M; j++)
            {
                tmp = 0;
                for (int k = 0; k < F; k++)
                    tmp += ((a.x[i][k] % Pr) * (inv_a->x[k][j] % Pr)) % Pr;
                fls.x[i][j] = tmp % Pr;
                // printf("%d ", z.x[i][j]);
            }
            // printf("\n");
        }
        return; // z;
    }
}
