#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv)
{
    FILE *plik;
    int i, j;
    srand(time(NULL));

    int rozmiar = argv[1] ? atoi(argv[1]) : 2000;

    plik = fopen("liczby.txt", "w");
    if (plik == NULL)
    {
        printf("Blad otwarcia pliku \"liczby.txt\"\n");
        exit(0);
    }

    int counter = 0;
    for (i = 0; i < rozmiar * 2; i++)
    {
        for (j = 0; j < rozmiar; j++)
        {
            fprintf(plik, "%6.2f", (float)(rand() % 10) + 1.0);
        }
        fprintf(plik, "\n");
    }
    fclose(plik);
    return 0;
}