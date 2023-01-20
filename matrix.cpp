#include "matrix.h"
#include <stdio.h>
#include <memory.h>
#include <stdlib.h>

void Complement( double** pTabO, double** pTabI, int nRow, int nCol, int nDim );
void ComplMatrix( double** pTabD, double** pTab, int nDim );

int CreateMatrix( double*** pTab, int nDim )
{
    double** p = *pTab = (double**)malloc( nDim * sizeof( double* ) );
    if( !*pTab ) return 0;
    memset( *pTab, 0, nDim * sizeof( double* ) );

    for( int i = 0; i < nDim; i++, p++ )
    {
        *p = (double*)malloc( nDim * sizeof( double ) );
        if( !*p ) return 0;
        memset( *p, 0, nDim * sizeof( double ) );
    }
    return 1;
}

void DeleteMatrix( double*** pTab, int nDim )
{
    double** p = *pTab;
    for( int i = 0; i < nDim; i++, p++ )
    {
        free( *p );
    }
    free( *pTab );
    *pTab = NULL;
}

void TransMatrix( double** pTab, int nDim )
{
    double** pW = pTab;

    for( int i = 0; i < nDim - 1; i++, pW++ )
    {
        double* pK = *pW + i + 1;

        for( int j = i + 1; j < nDim; j++, pK++ )
        {
            double temp = *pK;
            *pK = pTab[j][i];
            pTab[j][i] = temp;
        }
    }
}

void InverseMatrix( double** pInv, double** pTab, int nDim, double det )
{
    ComplMatrix( pInv, pTab, nDim );
    TransMatrix( pInv, nDim );
    for( int i = 0; i < nDim; i++ )
    {
        double* p = *pInv++;
        for( int j = 0; j < nDim; j++ )
        {
            *p++ /= det;
        }
    }
}

double Det( double** pTab, int nDim )
{
    if( nDim == 1 ) return **pTab;
    if( nDim == 2 ) return **pTab * pTab[1][1] - pTab[0][1] * pTab[1][0];

    double** pTabO = NULL;
    if( !CreateMatrix( &pTabO, nDim - 1 ) )
    {
        perror( "ERROR: Nie przydzielono pamieci dla macierzy!\n" );
        return 0;
    }

    double res = 0;
    double* p = *pTab;
    int sgn = 1;

    for( int i = 0; i < nDim; i++ )
    {
        Complement( pTabO, pTab, 0, i, nDim );
        res += *p++ * Det( pTabO, nDim - 1 ) * sgn;
        sgn = -sgn;
    }

    DeleteMatrix( &pTabO, nDim - 1 );

    return res;
}

void LayoutEqu( double** pInv, double* pB, double* pRes, int nDim )
{
    for( int i = 0; i < nDim; i++, pRes++ )
    {
        double* p1 = *pInv++;
        double* p2 = pB;

        for( int j = 0; j < nDim; j++ )
        {
            *pRes += *p1++ * *p2++;
        }
    }
}

void PrintMatrix( double** pTab, int nDim )
{
    for( int i = 0; i < nDim; i++ )
    {
        double* p = *pTab++;
        for( int j = 0; j < nDim; j++ )
        {
            printf( "%lf ", *p++ );
        }
        printf( "\n" );
    }
}

void Complement( double** pTabO, double** pTabI, int nRow, int nCol, int nDim )
{
    for( int i = 0; i < nDim; i++, pTabI++ )
    {
        if( i == nRow ) continue;

        double* pI = *pTabI;
        double* pO = *pTabO++;

        for( int j = 0; j < nDim; j++, pI++ )
        {
            if( j == nCol ) continue;
            *pO++ = *pI;
        }
    }
}

void ComplMatrix( double** pTabD, double** pTab, int nDim )
{
    double** pTabC = NULL;
    if( !CreateMatrix( &pTabC, nDim - 1 ) )
    {
        perror( "ERROR: Nie przydzielono pamieci dla macierzy!\n" );
        return;
    }

    for( int i = 0; i < nDim; i++ )
    {
        int sgn = ( i % 2 ) ? -1 : 1;
        double* p = *pTabD++;

        for( int j = 0; j < nDim; j++, p++ )
        {
            Complement( pTabC, pTab, i, j, nDim );
            *p = Det( pTabC, nDim - 1 ) * sgn;
            sgn = -sgn;
        }
    }

    DeleteMatrix( &pTabC, nDim - 1 );
}