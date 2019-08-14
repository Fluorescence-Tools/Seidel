

#define SMBUS_API __declspec(dllexport)
#define SMB_MAX_DATA_SIZE 5
#include <stdlib.h>

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

#include <stdlib.h>

__declspec(dllexport) void __stdcall Foo(unsigned char** ppMem, int* pSize)
{
    char i;
    *pSize = 4;
    *ppMem = (unsigned char*)malloc(*pSize);
    for(i = 0; i < *pSize; i++)
        (*ppMem)[i] = i;
}