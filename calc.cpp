#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "matrix.h"

#define _DEBUG_

// input data in .txt file passed as arg of main

void ReadData( FILE* fin, double** pMatrix, double* b, int nDim );
int CreateVector( double** pVec, int nDim );
void PrintVector( double* pVec, int nDim );

int main( int argc, char* argv[] )
{
    if( argc != 2 )
    {
        printf( "Usage: %s <input_file>\n", argv[0] );
        return 1;
    }

    FILE* fin = NULL;

    if( !( fin = fopen( argv[1], "r" ) ) )
    {
        perror( "Error open input file\n" );
        return 2;
    }

    int nDim = 0;
    fscanf( fin, "%d\n", &nDim );
#ifdef _DEBUG_
    printf( "nDim = %d\n", nDim );
#endif

    double** pMatrix = NULL;
    if( !CreateMatrix( &pMatrix, nDim ) )
    {
        perror( "ERROR: Nie przydzielono pamieci dla macierzy!\n" );
        return 3;
    }

    double* pB = NULL;
    if( !CreateVector( &pB, nDim ) )
    {
        perror( "ERROR: Nie przydzielono pamieci dla wektora wyrazow wolnych!\n" );
        return 4;
    }

    ReadData( fin, pMatrix, pB, nDim );
    fclose( fin );

#ifdef _DEBUG_
    printf( "Macierz po  wczytaniu\n" );
    PrintMatrix( pMatrix, nDim ); 
    printf( "Wektor po  wczytaniu\n" );
    PrintVector( pB, nDim );
#endif

    if( nDim == 1 )
    {
        printf( "Wektor wynikowy:\n%lf", *pB / **pMatrix );
        DeleteMatrix( &pMatrix, nDim );
        free( pB );
        return 0;
    }

    double det = Det( pMatrix, nDim );

#ifdef _DEBUG_
    printf( "det = %lf\n", det );
#endif
    if( fabs( det ) < 1e-300 )
    {
        perror( "ERROR: Wyznacznik rowny zero!\n" );
        return 5;
    }

    double* pRes = NULL;
    if( !CreateVector( &pRes, nDim ) )
    {
        perror( "ERROR: Nie przydzielono pamieci dla wektora wynikowego!\n" );
        return 6;
    }

    double** pInv = NULL;
    if( !CreateMatrix( &pInv, nDim ) )
    {
        perror( "ERROR: Nie przydzielono pamieci dla macierzy odwrotnej!\n" );
        return 7;
    }

    InverseMatrix( pInv, pMatrix, nDim, det );
#ifdef _DEBUG_
    printf( "Macierz odwrotna\n" );
    PrintMatrix( pInv, nDim );
#endif

    LayoutEqu( pInv, pB, pRes, nDim );

    printf( "Wektor wynikowy:\n" );
    PrintVector( pRes, nDim );

    DeleteMatrix( &pInv, nDim );
    DeleteMatrix( &pMatrix, nDim );
    free( pB );
    free( pRes );

    return 0;
}

void ReadData( FILE* fin, double** pMatrix, double* b, int nDim )
{
    for( int i = 0; i < nDim; i++ )
    {
        double* p = *pMatrix++;
        for( int j = 0; j < nDim; j++ )
        {
            fscanf( fin, "%lf", p++ );
        }
        fscanf( fin, "%lf", b++ );
    }
}

int CreateVector( double** pVec, int nDim )
{
    *pVec = (double*)malloc( nDim * sizeof( double ) );
    if( !*pVec ) return 0;
    memset( *pVec, 0, nDim * sizeof( double ) );
    return 1;
}

void PrintVector( double* pVec, int nDim )
{
    for( int i = 0; i < nDim; i++ )
    {
        printf( "%lf\n", *pVec++ );
    }
}