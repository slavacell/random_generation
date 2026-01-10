#include <iostream>
#include <vector>
#include <cmath>
//#include <math.h>
//#include <string>

#define eps (0.1e-9)
#define pi (3.141592653589793)

using namespace std;

double rand_1(double x)
{
    x = (3.9 + (0.5 - x) / 10) * x * (10 - x);
    x *= 10;
    return fabs(x - int(x));
}

double rand_2(double x, double y, double seed, double k)
{
    k = sin(109 * x + 113 * y + 103 * k);
    k *= cos(seed) * pi;
    return rand_1(fabs(k - int(k)));
}

double erfinv(double a)
{
    double p, r, t;
    t = fmaf(a, 0.0 - a, 1.0);
    t = log(t);
    if (fabs(t) > 6.125) { // maximum ulp error = 2.35793
        p =              3.03697567e-10; //  0x1.4deb44p-32
        p = fmaf (p, t,  2.93243101e-8); //  0x1.f7c9aep-26
        p = fmaf (p, t,  1.22150334e-6); //  0x1.47e512p-20
        p = fmaf (p, t,  2.84108955e-5); //  0x1.dca7dep-16
        p = fmaf (p, t,  3.93552968e-4); //  0x1.9cab92p-12
        p = fmaf (p, t,  3.02698812e-3); //  0x1.8cc0dep-9
        p = fmaf (p, t,  4.83185798e-3); //  0x1.3ca920p-8
        p = fmaf (p, t, -2.64646143e-1); // -0x1.0eff66p-2
        p = fmaf (p, t,  8.40016484e-1); //  0x1.ae16a4p-1
    } else { // maximum ulp error = 2.35002
        p =              5.43877832e-9;  //  0x1.75c000p-28
        p = fmaf (p, t,  1.43285448e-7); //  0x1.33b402p-23
        p = fmaf (p, t,  1.22774793e-6); //  0x1.499232p-20
        p = fmaf (p, t,  1.12963626e-7); //  0x1.e52cd2p-24
        p = fmaf (p, t, -5.61530760e-5); // -0x1.d70bd0p-15
        p = fmaf (p, t, -1.47697632e-4); // -0x1.35be90p-13
        p = fmaf (p, t,  2.31468678e-3); //  0x1.2f6400p-9
        p = fmaf (p, t,  1.15392581e-2); //  0x1.7a1e50p-7
        p = fmaf (p, t, -2.32015476e-1); // -0x1.db2aeep-3
        p = fmaf (p, t,  8.86226892e-1); //  0x1.c5bf88p-1
    }
    r = a * p;
    return r;
}

struct perlin
{
    int m, n;
    vector<vector<vector<double>>> noise;
    vector<int> grid;
    vector<double> power;

    perlin(int M, int N, vector<int> G) : m(M), n(N), noise(M + 1, vector<vector<double>>(N + 1, vector<double>(N + 1, 0))), grid(G), power(m, double(1.0 / m)) {}
    perlin(int M, int N, vector<int> G, vector<double> P) : m(M), n(N), noise(M + 1, vector<vector<double>>(N + 1, vector<double>(N + 1, 0))), grid(G), power(P) {}

    void generate(double x, double y, double seed, double k)
    {
        for (int i = 0; i < m; i++)
        {
            int dn = n / grid[i];
            double a = rand_2(x, y, seed, k * (1 + sin(i)));
            for (int j = 0; j < grid[i]; j++)
            {
                for (int jj = 0; jj < grid[i]; jj++)
                {
                    a = rand_1(a);
                    noise[i][j * dn][jj * dn] = a;
                }
            }

            a = rand_2(x + 1, y, seed, k * (1 + sin(i)));
            for (int j = 0; j < grid[i]; j++)
            {
                a = rand_1(a);
                noise[i][n][j * dn] = a;
            }

            a = rand_2(x + 1, y + 1, seed, k * (1 + sin(i)));
            a = rand_1(a);
            noise[i][n][n] = a;

            a = rand_2(x, y + 1, seed, k * (1 + sin(i)));
            for (int j = 0; j < grid[i]; j++)
            {
                a = rand_1(a);
                noise[i][j * dn][n] = a;
                for (int jj = 0; jj < grid[i] - 1; jj++) {a = rand_1(a);}
            }
        }
    }

    void fill()
    {
        int i_, j_, dn;
        double a1, a2, a3, a4;
        for (int k = 0; k < m; k++)
        {
            dn = n / grid[k];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    a1 = noise[k][int(i / dn) * dn][int(j / dn) * dn];
                    a2 = noise[k][int(i / dn) * dn][int(j / dn) * dn + dn];
                    a3 = noise[k][int(i / dn) * dn + dn][int(j / dn) * dn];
                    a4 = noise[k][int(i / dn) * dn + dn][int(j / dn) * dn + dn];
                    i_ = dn - (i % dn); j_ = dn - (j % dn);
                    noise[k][i][j] = (a1 * i_ * j_ + a2 * i_ * (j % dn) + a3 * j_ * (i % dn) + a4 * (i % dn) * (j % dn)) / (dn * dn);
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                for (int k = 0; k < m; k++)
                {
                    noise[m][i][j] += noise[k][i][j] * power[k];
                }
            }
        }

    }

    void clear()
    {
        noise = vector<vector<vector<double>>>(m + 1, vector<vector<double>>(n + 1, vector<double>(n + 1, 0)));
    }

    void remake(int M, int N, vector<int> G)
    {
        m = M; n = N; grid = G;
        noise = vector<vector<vector<double>>>(M + 1, vector<vector<double>>(N + 1, vector<double>(N + 1, 0)));
        power = vector<double>(m, double(1.0 / m));
    }
    void remake(int M, int N, vector<int> G, vector<double> P)
    {
        m = M; n = N; grid = G; power = P;
        noise = vector<vector<vector<double>>>(M + 1, vector<vector<double>>(N + 1, vector<double>(N + 1, 0)));
    }

    double p(double P)
    {
        if (m == 1 && grid[0] == n){return P;}

        double d = 0.1 / m;
        return sqrt(2 * d) * erfinv(2 * P - 1) + 0.5;
    }

    double get(double i, double j)
    {
        if (0 <= i && i < n && 0 <= j && j < n){return noise[m][i][j];}
        else {return 0;}
    }

    void print(int k, int full)
    {
        if (full)
        {
            for (int i = 0; i < n + 1; i++)
            {
                cout << endl;
                for (int j = 0; j < n + 1; j++)
                {
                    cout << double(noise[k][i][j]) << " ";
                }
            }
            cout << endl;
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                cout << endl;
                for (int j = 0; j < n; j++)
                {
                    cout << double(noise[k][i][j]) << " ";
                }
            }
            cout << endl;
        }
    }
};


void cover(vector<vector<int>>& a, int xa, int ya, vector<vector<int>> b, int xb, int yb, int X, int Y)
{
    for (int i = 0; i < xb; i++)
    {
        for (int j = 0; j < yb; j++)
        {
            if (b[i][j] != 0)
            {
                if (0 <= X + i && X + i < xa && 0 <= Y + j && Y + j < ya)
                {
                    a[X + i][Y + j] = b[i][j];
                }
            }
        }
    }
}

void replace(vector<vector<int>>& a, int xa, int ya, int A, int B, int X, int Y)
{
    if (0 <= X && X < xa && 0 <= Y && Y < ya)
    {
        if (a[X][Y] == A)
        {
            a[X][Y] = B;
        }
    }
}

void generate(vector<vector<int>>& M, int dx, int dy, int x, int y, int seed)
{
    for (int i = 0; i < dx; i++)
    {
        for (int j = 0; j < dy; j++)
        {
            perlin p(3, 24, {3, 4, 6}, {0.4, 0.3, 0.3});
            p.generate(x / 24 + i, y / 24 + j, seed, 0.1);
            p.fill();
            double p0 = p.p(0.3);
            for (int ii = 0; ii < 24; ii++)
            {
                for (int jj = 0; jj < 24; jj++)
                {
                    if (p.get(ii, jj) < p0){M[24 * i + ii][24 * j + jj] = 249;}
                    else{M[24 * i + ii][24 * j + jj] = 250;}
                }
            }
        }
    }
    for (int i = 0; i < dx; i++)
    {
        for (int j = 0; j < dy; j++)
        {
            perlin p(3, 24, {3, 4, 6}, {0.4, 0.3, 0.3});
            p.generate(x / 24 + i, y / 24 + j, seed, 0.1);
            p.fill();
            double p0 = p.p(0.3);

            perlin a(1, 24, {24});
            a.generate(x / 24 + i, y / 24 + j, seed, 0.2);
            a.fill();
            double a0 = a.p(0.15);

            perlin b(3, 24, {24, 3, 2}, {0.5, 0.25, 0.25});
            b.generate(x / 24 + i, y / 24 + j, seed, 0.2);
            b.fill();
            double b0 = b.p(0.05);

            for (int ii = 0; ii < 24; ii++)
            {
                for (int jj = 0; jj < 24; jj++)
                {
                    if (p.get(ii, jj) < p0)
                    {
                        if (a.get(ii, jj) < a0){M[24 * i + ii][24 * j + jj] = 245;}
                    }
                    else
                    {
                        if (b.get(ii, jj) < b0){M[24 * i + ii][24 * j + jj] = 245;}
                    }
                }
            }
        }
    }

    vector<vector<int>> norm{{0, 1, -1, 0}, {1, 0, 0, -1}};

    for (int i = -1; i < dx + 1; i++)
    {
        for (int j = -1; j < dy + 1; j++)
        {
            perlin p(3, 24, {3, 4, 6}, {0.4, 0.3, 0.3});
            p.generate(x / 24 + i, y / 24 + j, seed, 0.1);
            p.fill();
            double p0 = p.p(0.3);

            double a = rand_2(x / 24 + i, y / 24 + j, seed, 0.3);
            if (a < 0.47)
            {
                double X = 24 * rand_1(a), Y = 24 * rand_1(rand_1(a));
                X = int(X); Y = int(Y);
                cover(M, 24 * dx, 24 * dy, {{219, 219, 219, 219, 219},
                                            {219, 0, 0, 0, 219},
                                            {219, 0, 0, 0, 219},
                                            {219, 0, 0, 0, 219},
                                            {219, 219, 219, 219, 219}}, 5, 5, 24 * i + X - 2, 24 * j + Y - 2);
                for (int ii = 0; ii < 3; ii++)
                {
                    for (int jj = 0; jj < 3; jj++)
                    {
                        replace(M, 24 * dx, 24 * dy, 245, p.get(X - 1 + ii, Y - 1 + jj) < p0 ? 249 : 250, 24 * i + X - 1 + ii, 24 * j + Y - 1 + jj);
                    }
                }
                if (p.get(X, Y) >= p0)
                {
                    cover(M, 24 * dx, 24 * dy, {{245, 245, 245},
                                            {245, 0, 245},
                                            {245, 245, 245}}, 3, 3, 24 * i + X - 1, 24 * j + Y - 1);
                }
                double b = 4 * rand_2(X, Y, seed, 0.4);
                replace(M, 24 * dx, 24 * dy, 219, p.get(X + norm[0][b] * 2, Y + norm[1][b] * 2) < p0 ? 249 : 250, 24 * i + X + norm[0][b] * 2, 24 * j + Y + norm[1][b] * 2);



            }
        }
    }
}


int main()
{
/*
    for (int i = 0; i < 10000000; i++)
    {
        double a = rand_2(1000, 1000, i, 0.3),
        b = rand_2(1000, 1001, i, 0.3),
        c = rand_2(1001, 1000, i, 0.3),
        d = rand_2(1001, 1001, i, 0.3);
        if (a < 0.47 && b < 0.47 && c < 0.47 && d < 0.47)
        {
            double Xa = 24 * rand_1(a), Ya = 24 * rand_1(rand_1(a)); Xa = int(Xa); Ya = int(Ya);
            double Xb = 24 * rand_1(b), Yb = 24 * rand_1(rand_1(b)); Xb = int(Xb); Yb = int(Yb);
            double Xc = 24 * rand_1(c), Yc = 24 * rand_1(rand_1(c)); Xc = int(Xc); Yc = int(Yc);
            double Xd = 24 * rand_1(d), Yd = 24 * rand_1(rand_1(d)); Xd = int(Xd); Yd = int(Yd);
            if (Xa > 17 && Ya > 17 && Xb > 17 && Yb < 6 && Xc < 6 && Yc > 17 && Xd < 6 && Yd < 6){cout << i << endl;}
        }
    }
*/


    while (true)
    {
        int x, y, seed, dx, dy, s;
        cout << "dx dy x y seed /\n6 3 0 0 12346789 0\n";
        cin >> dx >> dy >> x >> y >> seed >> s;

        vector<vector<int>> M(dx * 24, vector<int>(dy * 24, 32));


        generate(M, dx, dy, x, y, seed);


        for (int i = 0; i < dx * 24; i++)
        {
            for (int j = 0; j < dy * 24; j++)
            {
                cout << char(M[i][j]);
                if (s && j % 24 == 23){cout << " ";}
            }
            cout << endl;
            if (s && i % 24 == 23){cout << endl;}
        }
    }







    //perlin P(3, 16, {8, 2, 1}, {0.16, 0.42, 0.42});

    //p.generate(1488, 1488, 1488, 1);
    //p.generate(1, 1, 1, 1);
    //p.fill();

    //double p1 = p.p(0.333);
    /*vector<vector<perlin>> p(20, vector<perlin>(5, P));
    for (int i = 0; i < 20; i++)
    {
        for (int j = 0; j < 5; j++)
        {
            p[i][j].generate(i, j, 1488, 1);
            //p[i][j].print(0, 1);
            p[i][j].fill();
        }
    }
    double p1 = p[0][0].p(0.25), p2 = p[0][0].p(0.5), p3 = p[0][0].p(0.75), pp;

    for (int i = 0; i < 20 * 16; i++)
    {
        for (int j = 0; j < 5 * 16; j++)
        {
            pp = p[i / 16][j / 16].noise[3][i % 16][j % 16];
            if (pp < p1){cout << "# ";}
            else if (pp < p2){cout << "* ";}
            else if (pp < p3){cout << ". ";}
            else {cout << "  ";}
        }
        cout << endl;
    }*/


    /*for (int i = 0; i < 16; i++)
    {
        for (int j = 0; j < 16; j++)
        {
            if (p.noise[3][i][j] < 0.333) {cout << ". ";}
            else if (p.noise[3][i][j] < 00.667) {cout << "+ ";}
            else {cout << "- ";}

            //cout << " ";
        }
        cout << endl;
    }*/






    /*vector<int> b(1000, 0);
    double a = 0.6, s = 0;

    for (int i = 0; i < 100000000; i++)
    {
        a = rand_1(a);
        //cout << a << endl;
        b[a * 1000]++;
    }
    for (int i = 0; i < 1000; i++)
    {
        cout << b[i] << "  ";
        s += abs(b[i] - 100000);
    }
    cout << endl << s / 1000 << endl;*/
    return 0;
}
