#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>

#ifdef __WIN32__
#include<direct.h>
#else
#include <sys/stat.h>
#endif

template <typename T>
T* init1D(int m, T init)
{
    T* A = new T[m];
    for (int i = 0; i < m; i++)
    {
        A[i] = init;
    }
    return A;
}

template <typename T>
T** init2D(int m, int n, T init)
{
    T** A = new T*[m];
    for (int i = 0; i < m; i++)
    {
        A[i] = new T[n];
        for (int j = 0; j < n; j++)
        {
            A[i][j] = init;
        }
    }
    return A;
}

template <typename T>
void free1D(T* A)
{
    delete[] A;
}

template <typename T>
void free2D(T** A, int m)
{
    for (int i = 0; i < m; i++)
    {
        delete[] A[i];
    }
    delete[] A;
}

int MKDIR(std::string fname,int acc)
{

    #ifdef __WIN32__
                return mkdir(fname.c_str());
    #else
                return mkdir(fname.c_str(),acc);
    #endif
}
