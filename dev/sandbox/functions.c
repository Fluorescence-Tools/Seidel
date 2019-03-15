
#include <stdio.h>

typedef struct Options{
    char mode[50];
}Option;

void hello(Option option);

double dprod(double *x, int n);

void dcumsum(double *a, double *b, int n);

void
hello(Option option)
{
    printf("C says hello\n");
    printf("mode is %s\n", option.mode);
}


double 
dprod(double *x, int n)
{
    int i;
    double y = 1.0;
    
    for (i = 0; i < n; i++)
    {
        y *= x[i];
    }

    return y;
}

void
dcumsum(double *a, double *b, int n)
{
    int i;
    
    b[0] = a[0];
    for (i = 1; i < n; i++)
    {
        b[i] = a[i] + b[i-1];
    }
}

#define SMBUS_API __declspec(dllexport)
#define SMB_MAX_DATA_SIZE 5

typedef void* SMBUS_HANDLE;

typedef struct _SMB_REQUEST
{
    unsigned char Address;
    unsigned char Command;
    unsigned char BlockLength;
    unsigned char Data[SMB_MAX_DATA_SIZE];
} SMB_REQUEST;

SMBUS_API int SmBusReadByte(SMBUS_HANDLE handle,SMB_REQUEST *request)
{
    unsigned char i;
    for(i = 0; i < request->BlockLength; i++)
        request->Data[i] = i;
    return request->BlockLength;
}

SMBUS_API SMBUS_HANDLE OpenSmbus(void)
{
    return (void*)0x12345678;
}