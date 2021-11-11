/*! \file magma.c Прямая реализация алгоритма по ГОСТ Р 34.12?2015 п. 5 Алгоритм блочного шифрования с длиной блока n = 64 бит */

#include "magma.h"

#include <stdio.h>
#include <inttypes.h>

/// 5.1.1 Нелинейное биективное преобразование
const static uint8_t p[8][16] =
{
    {12, 4, 6, 2, 10, 5, 11, 9, 14, 8, 13, 7, 0, 3, 15, 1},
    {6, 8, 2, 3, 9, 10, 5, 12, 1, 14, 4, 7, 11, 13, 0, 15},
    {11, 3, 5, 8, 2, 15, 10, 13, 14, 1, 7, 4, 12, 9, 6, 0},
    {12, 8, 2, 1, 13, 4, 15, 6, 7, 0, 10, 5, 3, 14, 9, 11},
    {7, 15, 5, 10, 8, 1, 6, 13, 0, 9, 3, 14, 11, 4, 2, 12},
    {5, 13, 15, 6, 9, 2, 12, 10, 11, 7, 8, 1, 4, 3, 14, 0},
    {8, 14, 2, 5, 6, 9, 1, 12, 15, 4, 11, 0, 13, 10, 3, 7},
    {1, 7, 14, 13, 0, 5, 8, 3, 4, 15, 10, 6, 9, 12, 11, 2}
};

/// 5.2 Преобразования
/// Преобразование t: V32 ? V32
static inline uint32_t t (uint32_t a)
{
    return (uint32_t)p[0][a&0xF] |
           ((uint32_t)p[1][(a>>4)&0xF] << 4) |
           ((uint32_t)p[2][(a>>8)&0xF] << 8) |
           ((uint32_t)p[3][(a>>12)&0xF] << 12) |
           ((uint32_t)p[4][(a>>16)&0xF] << 16) |
           ((uint32_t)p[5][(a>>20)&0xF] << 20) |
           ((uint32_t)p[6][(a>>24)&0xF] << 24) |
           ((uint32_t)p[7][(a>>28)&0xF] << 28);
}
/// g[k]: V32 ? V32
static inline uint32_t g (uint32_t a, uint32_t k)
{
    uint32_t v = t(a + k);
    return (v << 11) | (v >> (32 - 11));
}
/// G[k]: V32 ? V32 ? V32 ? V32
static inline void G(uint32_t* a1, uint32_t* a0, uint32_t k)
{
    uint32_t a0n = (*a1) ^ g((*a0), k);
    *a1 = *a0;
    *a0 = a0n;
}
///G?[k]: V32 ? V32 ? V6
static inline uint64_t Gz(uint32_t* a1, uint32_t* a0, uint32_t k)
{
   G(a1,a0,k);
   return ((uint64_t)(*a0) << 32) | (*a1);
}

///5.3 Алгоритм развертывания ключа
static const uint32_t K[8] =
{
   0xffeeddcc,
   0xbbaa9988,
   0x77665544,
   0x33221100,
   0xf0f1f2f3,
   0xf4f5f6f7,
   0xf8f9fafb,
   0xfcfdfeff,
};

static inline uint32_t k(uint8_t i)
{
    return (i < 24) ? K[i&7] : K[(~i)&7];
}
///5.4.1 Алгоритм зашифрования
uint64_t E(uint64_t a)
{
    uint32_t a1 = a >> 32;
    uint32_t a0 = a;
    uint32_t i;
    for (i = 0; i < 31; i++)
        G(&a1, &a0, k(i));
    return Gz(&a1, &a0, k(31));
}
///5.4.2 Алгоритм расшифрования
uint64_t D(uint64_t a)
{
    uint32_t a1 = a >> 32;
    uint32_t a0 = a;
    uint32_t i;
    for (i = 0; i < 31; i++)
        G(&a1, &a0, k(31 - i));
    return Gz(&a1, &a0, k(0));
}

/// Контрольные примеры
/// А.2 Алгоритм блочного шифрования с длиной блока n = 64 бит
void magmaCheck()
{
    ///А.2.1 Преобразование t
    {
        printf("А.2.1 Преобразование t\n");
        #define t_check(a,b) {uint32_t _t = t(a); printf("t(%08x)=%08x \t - [%s]\n", (a), _t, (_t == (b))?"Ok":"Error");}
        t_check(0xfdb97531, 0x2a196f34);
        t_check(0x2a196f34, 0xebd9f03a);
        t_check(0xebd9f03a, 0xb039bb3d);
        t_check(0xb039bb3d, 0x68695433);
        #undef t_check
    }
    ///A.2.2 Преобразование g
    {
        printf("A.2.2 Преобразование g\n");
        #define g_check(k,a,b) {uint32_t _g = g(a,k); printf("g[%08x](%08x)=%08x \t - [%s]\n",(k),(a), _g, (_g == (b))?"Ok":"Error");}
        g_check(0x87654321, 0xfedcba98, 0xfdcbc20c);
        g_check(0xfdcbc20c, 0x87654321, 0x7e791a4b);
        g_check(0x7e791a4b, 0xfdcbc20c, 0xc76549ec);
        g_check(0xc76549ec, 0x7e791a4b, 0x9791c849);
        #undef g_check
    }
    /// A.2.3 Алгоритм развертывания ключа
    {
        printf("A.2.3 Алгоритм развертывания ключа\n");
        uint8_t i,j;
        for(j = 0; j < 8; j++)
        {
            for(i = 0; i < 4; i++)
            {
                uint8_t id = i*8 + j;
                printf("K%u = %08x\t", id + 1, k(id));
            }
            printf("\n");
        }
    }
    /// A.2.4 Алгоритм зашифрования
    {
        printf("A.2.4 Алгоритм зашифрования\n");
        uint64_t a = 0xfedcba9876543210;
        uint64_t b = E(a);
        printf("a=%016"PRIx64"\tb=%016"PRIx64"\t[%s]\n",a,b,(b == 0x4ee901e5c2d8ca3d)?"Ok":"Error");
    }
    /// A.2.5 Алгоритм расшифрования
    {
        printf("A.2.5 Алгоритм расшифрования\n");
        uint64_t b =  0x4ee901e5c2d8ca3d;
        uint64_t a = D(b);
        printf("b=%016"PRIx64"\ta=%016"PRIx64"\t[%s]\n",b,a,(a == 0xfedcba9876543210)?"Ok":"Error");
    }
}
