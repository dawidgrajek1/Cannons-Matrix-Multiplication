#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define N 1000
#define P 4
#define PP 2

float a_global[N][N], b_global[N][N], c_global[N][N], c_sek[N][N];
float a[N / PP][N / PP], b[N / PP][N / PP], c[N / PP][N / PP];
float aa[N / PP][N / PP], bb[N / PP][N / PP];
float (*psa)[N / PP], (*psb)[N / PP], (*pra)[N / PP], (*prb)[N / PP];
float tmpa[N / PP][N / PP], tmpb[N / PP][N / PP], tmpc[N / PP][N / PP];

double startwtime1, startwtime2, endwtime;

int modulo(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

int shift_left(int rank, int p, int pp, int times)
{
    if (times == 0)
        return rank;
    int out = rank;
    for (int i = 0; i < times; i++)
        out = modulo(out - 1, pp) + ((out / pp) * pp);
    return out;
}
int shift_right(int rank, int p, int pp, int times)
{
    if (times == 0)
        return rank;
    int out = rank;
    for (int i = 0; i < times; i++)
        out = ((out + 1) % pp) + ((out / pp) * pp);
    return out;
}
int shift_up(int rank, int p, int pp, int times)
{
    if (times == 0)
        return rank;
    int out = rank;
    for (int i = 0; i < times; i++)
        out = (out + p - pp) % p;
    return out;
}
int shift_down(int rank, int p, int pp, int times)
{
    if (times == 0)
        return rank;
    int out = rank;
    for (int i = 0; i < times; i++)
        out = (out + PP) % P;
    return out;
}

int main(int argc, char **argv)
{
    FILE *plik;
    FILE *plik_out;

    int my_rank, ncpus;
    int row, col, mod = 0;
    int data_received = -1;
    int tag = 101;
    int koniec;

    MPI_Status statSend[2], statRecv[2], statRecvCollect[P];
    MPI_Request reqSend[2], reqRecv[2], reqSendCollect[P], reqRecvCollect[P];

    MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpus);

    row = my_rank / PP;
    col = my_rank % PP;

    int upper_neighbor = (my_rank + P - PP) % P;
    int lower_neighbor = (my_rank + PP) % P;
    int right_neighbor = ((my_rank + 1) % PP) + ((my_rank / PP) * PP);
    int left_neighbor = modulo(my_rank - 1, PP) + ((my_rank / PP) * PP);

    if (my_rank == 0)
        printf("obliczenia metoda Cannona dla tablicy %d x %d elementow \n", N, N);

    if (my_rank == 0)
        startwtime1 = MPI_Wtime(); // czas w sekundach

    // wczytanie danych przez proces rank=0
    if (my_rank == 0)
    {
        plik = fopen("liczby.txt", "r");
        if (plik == NULL)
        {
            printf("Blad otwarcia pliku \"liczby.txt\"\n");
            koniec = 1;
            MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Finalize();
            exit(0);
        }
        else
        {
            koniec = 0;
            MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Bcast(&koniec, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (koniec)
        {
            MPI_Finalize();
            exit(0);
        }
    }

    if (ncpus != P)
    {
        if (my_rank == 0)
            printf("wywolano obliczenia iloczynu macierzy metoda cannona na %d procesach - uruchom mpiexec -n %d matrixmult\n", ncpus, P);
        MPI_Finalize();
        exit(0);
    }

    if (my_rank == 0)
    {
        // odczyt danych wejściowych tablicy A
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                fscanf(plik, "%f", &a_global[i][j]);
            }
        // odczyt danych wejściowych tablicy B
        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
            {
                fscanf(plik, "%f", &b_global[i][j]);
            }

        // proces 0 przypisuje początkowe wartości tablic A i B do tablic a i b
        for (int i = 0; i < N / PP; i++)
            for (int j = 0; j < N / PP; j++)
            {
                a[i][j] = a_global[i][j];
                b[i][j] = b_global[i][j];
            }

        // rozesłanie tablicy a zgodnie z dystrybucją początkową tablicy A
        // rozesłanie tablicy a zgodnie z dystrybucją początkową tablicy B
        for (int i = 1; i < P; i++)
        {
            for (int j = 0; j < N / PP; j++)
                for (int k = 0; k < N / PP; k++)
                {
                    tmpa[j][k] = a_global[(shift_right(i, P, PP, i / PP) / PP) * (N / PP) + j][(shift_right(i, P, PP, i / PP) % PP) * (N / PP) + k];
                    tmpb[j][k] = b_global[(shift_down(i, P, PP, i % PP) / PP) * (N / PP) + j][(shift_down(i, P, PP, i % PP) % PP) * (N / PP) + k];
                }
            MPI_Isend(tmpa, N * N / P, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqSend[0]);
            MPI_Isend(tmpb, N * N / P, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqSend[1]);
        }
    }
    else
    {
        MPI_Irecv(a, N * N / P, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[0]);
        MPI_Irecv(b, N * N / P, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqRecv[1]);
        MPI_Waitall(2, reqRecv, statRecv);
    }
    pra = aa;
    prb = bb;
    psa = a;
    psb = b;

    // przygotowanie lokalnej tablicy wynikowej
    for (int i = 0; i < N / PP; i++)
        for (int j = 0; j < N / PP; j++)
        {
            c[i][j] = 0;
        }

    if (my_rank == 0)
        startwtime2 = MPI_Wtime(); // czas w sekundach

    // obliczenia iloczynu macierzy zgodnie z algorytmem Cannona
    for (int kk = 0; kk < PP; kk++) // KOLENA ITERACJA PRZETWARZANIA (OBLICZENIA I KOMUNIKACJA)
    {
        // OBLICZENIA
        for (int i = 0; i < N / PP; i++)
            for (int k = 0; k < N / PP; k++)
                for (int j = 0; j < N / PP; j++)
                    c[i][j] += psa[i][k] * psb[k][j];

        // KOMUNIKAJCA
        //  MPI_Irecv(adres, ile_słów, typ_danych, odbiorca/nadawca, znacznik, zakres_procesów, Id_komunikacji);
        MPI_Irecv(pra, N * N / PP / PP, MPI_FLOAT, right_neighbor, tag, MPI_COMM_WORLD, &reqRecv[0]);
        MPI_Irecv(prb, N * N / PP / PP, MPI_FLOAT, lower_neighbor, tag, MPI_COMM_WORLD, &reqRecv[1]);
        MPI_Isend(psa, N * N / PP / PP, MPI_FLOAT, left_neighbor, tag, MPI_COMM_WORLD, &reqSend[0]);
        MPI_Isend(psb, N * N / PP / PP, MPI_FLOAT, upper_neighbor, tag, MPI_COMM_WORLD, &reqSend[1]);

        // OCZEKIWANIE NA KOMUNIKACJĘ ASYNCHRONICZNĄ
        MPI_Waitall(2, reqSend, statSend);
        MPI_Waitall(2, reqRecv, statRecv);

        // ZMIANA OBSZARÓW DANYCH LICZONYCH I PRZESYAŁANYCH
        if (mod = ((mod + 1) % 2))
        {
            pra = a;
            prb = b;
            psa = aa;
            psb = bb;
        }
        else
        {
            pra = aa;
            prb = bb;
            psa = a;
            psb = b;
        }
    }

    if (my_rank == 0)
    {
        endwtime = MPI_Wtime();
        printf("Calkowity czas przetwarzania wynosi %f sekund\n", endwtime - startwtime1);
        printf("Calkowity czas obliczen wynosi %f sekund\n", endwtime - startwtime2);
    }

    MPI_Isend(c, N * N / PP / PP, MPI_FLOAT, 0, tag, MPI_COMM_WORLD, &reqSendCollect[my_rank]);

    // test poprawnosci wyniku
    if (my_rank == 0)
    {
        // odbiór wyników obliczeń równoległych do globalnej tablicy wynikowej Cglob
        for (int i = 0; i < P; i++)
        {
            MPI_Irecv(tmpc, N * N / PP / PP, MPI_FLOAT, i, tag, MPI_COMM_WORLD, &reqRecvCollect[i]);
            MPI_Wait(&reqRecvCollect[i], &statRecvCollect[i]);

            for (int j = 0; j < N / PP; j++)
                for (int k = 0; k < N / PP; k++)
                {
                    c_global[(i / PP) * (N / PP) + j][(i % PP) * (N / PP) + k] = tmpc[j][k];
                }
            // printf("odebrano wynik od procesu %d\n", i);
        }

        // obliczenia sekwencyjne mnożenia tablic CSek=A*B
        startwtime2 = MPI_Wtime(); // czas w sekundach
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                c_sek[i][j] = 0;
                for (int k = 0; k < N; k++)
                {
                    c_sek[i][j] += a_global[i][k] * b_global[k][j];
                }
            }
        }
        endwtime = MPI_Wtime();
        printf("Calkowity czas obliczen sekwencyjnych wynosi %f sekund\n", endwtime - startwtime2);

        // porównanie poprawności obliczeń (Csek, Cglob) przy uwzględniniu progu poprawności
        int errors = 0;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                if (c_sek[i][j] != c_global[i][j] && (c_sek[i][j] / c_global[i][j] < 0.9 || c_sek[i][j] / c_global[i][j] > 1.1))
                {
                    errors++;
                }
            }
        }
        printf("Errors: %d (%.2f%%)\n", errors, (float)errors / (N * N) * 100);

        // save to output
        plik_out = fopen("result.txt", "w");
        if (plik_out == NULL)
        {
            printf("Blad otwarcia pliku \"result.txt\"\n");
            exit(0);
        }
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                fprintf(plik_out, "%6.1f ", c_global[i][j]);
            }
            fprintf(plik_out, "\n");
        }
        fclose(plik_out);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}
