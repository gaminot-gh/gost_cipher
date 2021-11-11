#include "kuznechik.h"

#include <stdio.h>
#include <inttypes.h>

///4.1.1 Нелинейное биективное преобразование
static const uint8_t p[256] =
{
    252, 238, 221, 17, 207, 110, 49, 22, 251, 196, 250, 218, 35, 197, 4, 77, 233,
    119, 240, 219, 147, 46, 153, 186, 23, 54, 241, 187, 20, 205, 95, 193, 249, 24, 101,
    90, 226, 92, 239, 33, 129, 28, 60, 66, 139, 1, 142, 79, 5, 132, 2, 174, 227, 106, 143,
    160, 6, 11, 237, 152, 127, 212, 211, 31, 235, 52, 44, 81, 234, 200, 72, 171, 242, 42,
    104, 162, 253, 58, 206, 204, 181, 112, 14, 86, 8, 12, 118, 18, 191, 114, 19, 71, 156,
    183, 93, 135, 21, 161, 150, 41, 16, 123, 154, 199, 243, 145, 120, 111, 157, 158, 178,
    177, 50, 117, 25, 61, 255, 53, 138, 126, 109, 84, 198, 128, 195, 189, 13, 87, 223,
    245, 36, 169, 62, 168, 67, 201, 215, 121, 214, 246, 124, 34, 185, 3, 224, 15, 236,
    222, 122, 148, 176, 188, 220, 232, 40, 80, 78, 51, 10, 74, 167, 151, 96, 115, 30, 0,
    98, 68, 26, 184, 56, 130, 100, 159, 38, 65, 173, 69, 70, 146, 39, 94, 85, 47, 140, 163,
    165, 125, 105, 213, 149, 59, 7, 88, 179, 64, 134, 172, 29, 247, 48, 55, 107, 228, 136,
    217, 231, 137, 225, 27, 131, 73, 76, 63, 248, 254, 141, 83, 170, 144, 202, 216, 133,
    97, 32, 113, 103, 164, 45, 43, 9, 91, 203, 155, 37, 208, 190, 229, 108, 82, 89, 166,
    116, 210, 230, 244, 180, 192, 209, 102, 175, 194, 57, 75, 99, 182
};

#define GF8_POLY 0b11000011 ///< p(x) = x8 + x7 + x6 + x + 1
static inline uint8_t gf8m(uint8_t a, uint8_t b)
{
    if(!a || !b) return 0;
    uint16_t r = 0;
    uint8_t i;
    for(i = 0; i < 8; i++)
    {
        r ^= (((uint8_t)(~((b & 1) - 1))) & a);
        b = b >> 1;
        a = (a << 1) ^ ((uint8_t)(~((a >> 7) - 1)) & GF8_POLY);
    }
    return r;
}
/// 4.1.2 Линейное преобразование
/// l: V8^16 ? V8
static inline uint8_t l(uint8_t a[16])
{
    return
        gf8m(148, a[15]) ^
        gf8m(32, a[14]) ^
        gf8m(133, a[13]) ^
        gf8m(16, a[12]) ^
        gf8m(194, a[11]) ^
        gf8m(192, a[10]) ^
        gf8m(1, a[9]) ^
        gf8m(251, a[8]) ^
        gf8m(1, a[7]) ^
        gf8m(192, a[6]) ^
        gf8m(194, a[5]) ^
        gf8m(16, a[4]) ^
        gf8m(133, a[3]) ^
        gf8m(32, a[2]) ^
        gf8m(148, a[1]) ^
        gf8m(1, a[0]);
}
/// 4.2 Преобразования
///X[k]: V128 ? V128
static inline uint64_t* X(uint64_t* a, const uint64_t k[2])
{
    a[0] ^= k[0];
    a[1] ^= k[1];
    return a;
}
///S: V128 ? V128
static inline uint64_t* S(uint64_t* a)
{
    a[0] =  (uint64_t)p[(uint8_t)a[0]] |
            (((uint64_t)p[(uint8_t)(a[0] >> 8)]) << 8) |
            (((uint64_t)p[(uint8_t)(a[0] >> 16)]) << 16) |
            (((uint64_t)p[(uint8_t)(a[0] >> 24)]) << 24) |
            (((uint64_t)p[(uint8_t)(a[0] >> 32)]) << 32) |
            (((uint64_t)p[(uint8_t)(a[0] >> 40)]) << 40) |
            (((uint64_t)p[(uint8_t)(a[0] >> 48)]) << 48) |
            (((uint64_t)p[(uint8_t)(a[0] >> 56)]) << 56);
    a[1] =  (uint64_t)p[(uint8_t)a[1]] |
            (((uint64_t)p[(uint8_t)(a[1] >> 8)]) << 8) |
            (((uint64_t)p[(uint8_t)(a[1] >> 16)]) << 16) |
            (((uint64_t)p[(uint8_t)(a[1] >> 24)]) << 24) |
            (((uint64_t)p[(uint8_t)(a[1] >> 32)]) << 32) |
            (((uint64_t)p[(uint8_t)(a[1] >> 40)]) << 40) |
            (((uint64_t)p[(uint8_t)(a[1] >> 48)]) << 48) |
            (((uint64_t)p[(uint8_t)(a[1] >> 56)]) << 56);
    return a;
}
///S^-1: V128 ? V128
static uint8_t p_inv[256];
static void pInvCreate()
{
   uint16_t i;
   for(i = 0; i < 256; i++)
        p_inv[p[i]] = i;
}
static inline void S_inv(uint64_t* a)
{
    a[0] =  (uint64_t)p_inv[(uint8_t)a[0]] |
            (((uint64_t)p_inv[(uint8_t)(a[0] >> 8)]) << 8) |
            (((uint64_t)p_inv[(uint8_t)(a[0] >> 16)]) << 16) |
            (((uint64_t)p_inv[(uint8_t)(a[0] >> 24)]) << 24) |
            (((uint64_t)p_inv[(uint8_t)(a[0] >> 32)]) << 32) |
            (((uint64_t)p_inv[(uint8_t)(a[0] >> 40)]) << 40) |
            (((uint64_t)p_inv[(uint8_t)(a[0] >> 48)]) << 48) |
            (((uint64_t)p_inv[(uint8_t)(a[0] >> 56)]) << 56);
    a[1] =  (uint64_t)p_inv[(uint8_t)a[1]] |
            (((uint64_t)p_inv[(uint8_t)(a[1] >> 8)]) << 8) |
            (((uint64_t)p_inv[(uint8_t)(a[1] >> 16)]) << 16) |
            (((uint64_t)p_inv[(uint8_t)(a[1] >> 24)]) << 24) |
            (((uint64_t)p_inv[(uint8_t)(a[1] >> 32)]) << 32) |
            (((uint64_t)p_inv[(uint8_t)(a[1] >> 40)]) << 40) |
            (((uint64_t)p_inv[(uint8_t)(a[1] >> 48)]) << 48) |
            (((uint64_t)p_inv[(uint8_t)(a[1] >> 56)]) << 56);
}
///R: V128 ? V128
static inline void R(uint64_t* a)
{
    uint8_t _a[16] =
                    {
                        (uint8_t)a[0],
                        (uint8_t)(a[0] >> 8),
                        (uint8_t)(a[0] >> 16),
                        (uint8_t)(a[0] >> 24),
                        (uint8_t)(a[0] >> 32),
                        (uint8_t)(a[0] >> 40),
                        (uint8_t)(a[0] >> 48),
                        (uint8_t)(a[0] >> 56),
                        (uint8_t)a[1],
                        (uint8_t)(a[1] >> 8),
                        (uint8_t)(a[1] >> 16),
                        (uint8_t)(a[1] >> 24),
                        (uint8_t)(a[1] >> 32),
                        (uint8_t)(a[1] >> 40),
                        (uint8_t)(a[1] >> 48),
                        (uint8_t)(a[1] >> 56),
                    };
    a[0] = ((uint64_t)_a[1]) |
           (((uint64_t)_a[2]) << 8) |
           (((uint64_t)_a[3]) << 16) |
           (((uint64_t)_a[4]) << 24) |
           (((uint64_t)_a[5]) << 32) |
           (((uint64_t)_a[6]) << 40) |
           (((uint64_t)_a[7]) << 48) |
           (((uint64_t)_a[8]) << 56);
    a[1] = ((uint64_t)_a[9]) |
           (((uint64_t)_a[10]) << 8) |
           (((uint64_t)_a[11]) << 16) |
           (((uint64_t)_a[12]) << 24) |
           (((uint64_t)_a[13]) << 32) |
           (((uint64_t)_a[14]) << 40) |
           (((uint64_t)_a[15]) << 48) |
           (((uint64_t)l(_a)) << 56);
}
///L: V128 ? V128
static inline uint64_t* L(uint64_t* a)
{
    uint8_t i;
    for(i = 0; i < 16; i++)
        R(a);
    return a;
}
///R-1: V128 ? V128
static inline void R_inv(uint64_t* a)
{
    uint8_t _a[16] =
                    {
                        (uint8_t)(a[1] >> 56),
                        (uint8_t)a[0],
                        (uint8_t)(a[0] >> 8),
                        (uint8_t)(a[0] >> 16),
                        (uint8_t)(a[0] >> 24),
                        (uint8_t)(a[0] >> 32),
                        (uint8_t)(a[0] >> 40),
                        (uint8_t)(a[0] >> 48),
                        (uint8_t)(a[0] >> 56),
                        (uint8_t)a[1],
                        (uint8_t)(a[1] >> 8),
                        (uint8_t)(a[1] >> 16),
                        (uint8_t)(a[1] >> 24),
                        (uint8_t)(a[1] >> 32),
                        (uint8_t)(a[1] >> 40),
                        (uint8_t)(a[1] >> 48),
                    };
    a[0] = ((uint64_t)l(_a)) |
           (((uint64_t)_a[1]) << 8) |
           (((uint64_t)_a[2]) << 16) |
           (((uint64_t)_a[3]) << 24) |
           (((uint64_t)_a[4]) << 32) |
           (((uint64_t)_a[5]) << 40) |
           (((uint64_t)_a[6]) << 48) |
           (((uint64_t)_a[7]) << 56);
    a[1] = ((uint64_t)_a[8]) |
           (((uint64_t)_a[9]) << 8) |
           (((uint64_t)_a[10]) << 16) |
           (((uint64_t)_a[11]) << 24) |
           (((uint64_t)_a[12]) << 32) |
           (((uint64_t)_a[13]) << 40) |
           (((uint64_t)_a[14]) << 48) |
           (((uint64_t)_a[15]) << 56);
}
///L^-1: V128 ? V128
static inline uint64_t* L_inv(uint64_t* a)
{
    uint8_t i;
    for(i = 0; i < 16; i++)
        R_inv(a);
    return a;
}
///F [k]: V128 ? V128 ? V128 ? V128
static inline void F(uint64_t* a1, uint64_t* a0, const uint64_t k[2])
{
    uint64_t _a0[2] = {a0[0], a0[1]};
    a0[0] = a1[0]; a0[1] = a1[1];
    X(L(S(X(a1,k))), _a0);

}
/// 4.3 Алгоритм развертывания ключа
static uint64_t C[32][2];
static void cValCreate()
{
    uint64_t i;
    for(i = 0; i < 32; i++)
    {
        C[i][0] = i+1;
        C[i][1] = 0;
        L(C[i]);
    }
}

static const uint64_t K[2][2] =
{
    {0x0011223344556677, 0x8899aabbccddeeff},
    {0x0123456789abcdef, 0xfedcba9876543210}
};
static uint64_t Kexp[10][2];
static void keyExpand()
{
    Kexp[0][0] = K[0][0];
    Kexp[0][1] = K[0][1];
    Kexp[1][0] = K[1][0];
    Kexp[1][1] = K[1][1];
    uint8_t i;
    uint8_t j;
    for(i = 2; i < 10; i+=2)
    {
        Kexp[i][0] = Kexp[i-2][0];
        Kexp[i][1] = Kexp[i-2][1];
        Kexp[i|1][0] = Kexp[i-1][0];
        Kexp[i|1][1] = Kexp[i-1][1];
        for(j = 0; j < 8; j++)
        {
            F(Kexp[i], Kexp[i|1], C[(((i-2)<<2)) | j]);
       }
    }
}
///4.4.1 Алгоритм зашифрования
static inline uint64_t* E(uint64_t* a)
{
    uint8_t i;
    for(i = 0; i < 9; i++)
        L(S(X(a,Kexp[i])));
    return X(a,Kexp[9]);
}
///4.4.2 Алгоритм расшифрования
static inline uint64_t* D(uint64_t* a)
{
    uint8_t i;
    for(i = 0; i < 9; i++)
        S_inv(L_inv(X(a,Kexp[9-i])));
    return X(a,Kexp[0]);
}

/// Контрольные примеры
/// А.1 Алгоритм блочного шифрования с длиной блока n = 128 бит
void kuznechikCheck()
{
    /// Синтез обратной перестновки для S^-1 (S_inv)
    pInvCreate();
    /// Синтез констант для развертывания ключа
    cValCreate();
    /// Развертывание ключа
    keyExpand();
    /// А.1.1 Преобразование S
    {
        printf("А.1.1 Преобразование S\n");
        #define S_check(a1,a0,b1,b0) {uint64_t _a[2] = {a0, a1}; uint64_t _b[2] = {a0, a1}; uint64_t __b[2] = {b0, b1};\
                              S(_b);\
                              printf("S(%016"PRIx64"%016"PRIx64")\t=\t%016"PRIx64"%016"PRIx64"\t - [%s]\n",\
                              _a[1],_a[0], _b[1], _b[0], ((_b[0] == __b[0])&&(_b[1] == __b[1]))?"Ok":"Error");}
        S_check(0xffeeddccbbaa9988, 0x1122334455667700, 0xb66cd8887d38e8d7, 0x7765aeea0c9a7efc)
        S_check(0xb66cd8887d38e8d7, 0x7765aeea0c9a7efc, 0x559d8dd7bd06cbfe, 0x7e7b262523280d39)
        S_check(0x559d8dd7bd06cbfe, 0x7e7b262523280d39, 0x0c3322fed531e463, 0x0d80ef5c5a81c50b)
        S_check(0x0c3322fed531e463, 0x0d80ef5c5a81c50b, 0x23ae65633f842d29, 0xc5df529c13f5acda)
        #undef S_check
        printf("В дополнение к А.1.1 преобразование S^-1\n");
        #define S_inv_check(a1,a0,b1,b0) {uint64_t _a[2] = {a0, a1}; uint64_t _b[2] = {a0, a1}; uint64_t __b[2] = {b0, b1};\
                              S_inv(_b);\
                              printf("S^-1(%016"PRIx64"%016"PRIx64")\t=\t%016"PRIx64"%016"PRIx64"\t - [%s]\n",\
                              _a[1],_a[0], _b[1], _b[0], ((_b[0] == __b[0])&&(_b[1] == __b[1]))?"Ok":"Error");}
        S_inv_check(0x23ae65633f842d29, 0xc5df529c13f5acda, 0x0c3322fed531e463, 0x0d80ef5c5a81c50b)
        S_inv_check(0x0c3322fed531e463, 0x0d80ef5c5a81c50b, 0x559d8dd7bd06cbfe, 0x7e7b262523280d39)
        S_inv_check(0x559d8dd7bd06cbfe, 0x7e7b262523280d39, 0xb66cd8887d38e8d7, 0x7765aeea0c9a7efc)
        S_inv_check(0xb66cd8887d38e8d7, 0x7765aeea0c9a7efc, 0xffeeddccbbaa9988, 0x1122334455667700)
        #undef S_inv_check

    }
    ///А.1.2 Преобразование R
    {
        printf("А.1.2 Преобразование R\n");
        #define R_check(a1,a0,b1,b0) {uint64_t _a[2] = {a0, a1}; uint64_t _b[2] = {a0, a1}; uint64_t __b[2] = {b0, b1};\
                      R(_b);\
                      printf("R(%016"PRIx64"%016"PRIx64")\t=\t%016"PRIx64"%016"PRIx64"\t - [%s]\n",\
                      _a[1],_a[0], _b[1], _b[0], ((_b[0] == __b[0])&&(_b[1] == __b[1]))?"Ok":"Error");}

        R_check(0x0000000000000000, 0x0000000000000100, 0x9400000000000000, 0x0000000000000001)
        R_check(0x9400000000000000, 0x0000000000000001, 0xa594000000000000, 0x0000000000000000)
        R_check(0xa594000000000000, 0x0000000000000000, 0x64a5940000000000, 0x0000000000000000)
        R_check(0x64a5940000000000, 0x0000000000000000, 0x0d64a59400000000, 0x0000000000000000)
        #undef R_check
        printf("В дополнение к А.1.2 преобразование R^-1\n");
        #define R_inv_check(a1,a0,b1,b0) {uint64_t _a[2] = {a0, a1}; uint64_t _b[2] = {a0, a1}; uint64_t __b[2] = {b0, b1};\
              R_inv(_b);\
              printf("R^-1(%016"PRIx64"%016"PRIx64")\t=\t%016"PRIx64"%016"PRIx64"\t - [%s]\n",\
              _a[1],_a[0], _b[1], _b[0], ((_b[0] == __b[0])&&(_b[1] == __b[1]))?"Ok":"Error");}

        R_inv_check(0x0d64a59400000000, 0x0000000000000000, 0x64a5940000000000, 0x0000000000000000)
        R_inv_check(0x64a5940000000000, 0x0000000000000000, 0xa594000000000000, 0x0000000000000000)
        R_inv_check(0xa594000000000000, 0x0000000000000000, 0x9400000000000000, 0x0000000000000001)
        R_inv_check(0x9400000000000000, 0x0000000000000001, 0x0000000000000000, 0x0000000000000100)
        #undef R_inv_check

    }
    ///А.1.3 Преобразование L
    {
        printf("А.1.3 Преобразование L\n");
        #define L_check(a1,a0,b1,b0) {uint64_t _a[2] = {a0, a1}; uint64_t _b[2] = {a0, a1}; uint64_t __b[2] = {b0, b1};\
              L(_b);\
              printf("L(%016"PRIx64"%016"PRIx64")\t=\t%016"PRIx64"%016"PRIx64"\t - [%s]\n",\
              _a[1],_a[0], _b[1], _b[0], ((_b[0] == __b[0])&&(_b[1] == __b[1]))?"Ok":"Error");}

        L_check(0x64a5940000000000,0x0000000000000000, 0xd456584dd0e3e84c, 0xc3166e4b7fa2890d)
        L_check(0xd456584dd0e3e84c,0xc3166e4b7fa2890d, 0x79d26221b87b584c, 0xd42fbc4ffea5de9a)
        L_check(0x79d26221b87b584c,0xd42fbc4ffea5de9a, 0x0e93691a0cfc6040, 0x8b7b68f66b513c13)
        L_check(0x0e93691a0cfc6040,0x8b7b68f66b513c13, 0xe6a8094fee0aa204, 0xfd97bcb0b44b8580)
        #undef L_check
        printf("В дополнение к А.1.3 преобразование L^-1\n");
        #define L_inv_check(a1,a0,b1,b0) {uint64_t _a[2] = {a0, a1}; uint64_t _b[2] = {a0, a1}; uint64_t __b[2] = {b0, b1};\
              L_inv(_b);\
              printf("L^-1(%016"PRIx64"%016"PRIx64")\t=\t%016"PRIx64"%016"PRIx64"\t - [%s]\n",\
              _a[1],_a[0], _b[1], _b[0], ((_b[0] == __b[0])&&(_b[1] == __b[1]))?"Ok":"Error");}

        L_inv_check(0xe6a8094fee0aa204, 0xfd97bcb0b44b8580, 0x0e93691a0cfc6040,0x8b7b68f66b513c13)
        L_inv_check(0x0e93691a0cfc6040, 0x8b7b68f66b513c13, 0x79d26221b87b584c,0xd42fbc4ffea5de9a)
        L_inv_check(0x79d26221b87b584c, 0xd42fbc4ffea5de9a, 0xd456584dd0e3e84c,0xc3166e4b7fa2890d)
        L_inv_check(0xd456584dd0e3e84c, 0xc3166e4b7fa2890d, 0x64a5940000000000,0x0000000000000000)
        #undef L_inv_check
    }
    ///А.1.4 Алгоритм развертывания ключа
    {
        printf("А.1.4 Алгоритм развертывания ключа\n");
        printf("K=%016"PRIx64"%016"PRIx64"%016"PRIx64"%016"PRIx64"\n",
               K[0][1],K[0][0],K[1][1],K[1][0]);
        printf("C1=%016"PRIx64"%016"PRIx64" - [%s]\n",
               C[0][1],C[0][0], (C[0][1] == 0x6ea276726c487ab8)&&(C[0][0] == 0x5d27bd10dd849401)?"Ok":"Error");
        #define K_check(i,k1,k0) {uint64_t _k[2] = {k0, k1}; uint8_t _i = i - 1;\
          printf("K%u = %016"PRIx64"%016"PRIx64"\t - [%s]\n",\
          i,Kexp[_i][1],Kexp[_i][0], \
          ((_k[1] == Kexp[_i][1])&&(_k[0] == Kexp[_i][0]))?"Ok":"Error");}
        K_check(1 , 0x8899aabbccddeeff, 0x0011223344556677)
        K_check(2 , 0xfedcba9876543210, 0x0123456789abcdef)
        K_check(3 , 0xdb31485315694343, 0x228d6aef8cc78c44)
        K_check(4 , 0x3d4553d8e9cfec68, 0x15ebadc40a9ffd04)
        K_check(5 , 0x57646468c44a5e28, 0xd3e59246f429f1ac)
        K_check(6 , 0xbd079435165c6432, 0xb532e82834da581b)
        K_check(7 , 0x51e640757e8745de, 0x705727265a0098b1)
        K_check(8 , 0x5a7925017b9fdd3e, 0xd72a91a22286f984)
        K_check(9 , 0xbb44e25378c73123, 0xa5f32f73cdb6e517)
        K_check(10, 0x72e9dd7416bcf45b, 0x755dbaa88e4a4043)
        #undef K_check
    }
    ///А.1.5 Алгоритм зашифрования
    {
        printf("А.1.5 Алгоритм зашифрования\n");
        uint64_t a[2] = {0xffeeddccbbaa9988, 0x1122334455667700};
        printf("a=%016"PRIx64"%016"PRIx64"\n", a[1], a[0]);
        E(a);
        uint64_t b[2] = {0x5a468d42b9d4edcd, 0x7f679d90bebc2430};
        printf("b=%016"PRIx64"%016"PRIx64" - [%s]\n", a[1], a[0],
                ((a[1] == b[1])&&(a[0] == b[0]))?"Ok":"Error");
    }
    ///А.1.6 Алгоритм расшифрования
    {
        printf("А.1.6 Алгоритм расшифрования\n");
        uint64_t b[2] = {0x5a468d42b9d4edcd, 0x7f679d90bebc2430};
        printf("b=%016"PRIx64"%016"PRIx64"\n", b[1], b[0]);
        D(b);
        uint64_t a[2] = {0xffeeddccbbaa9988, 0x1122334455667700};
        printf("a=%016"PRIx64"%016"PRIx64" - [%s]\n", b[1], b[0],
                ((a[1] == b[1])&&(a[0] == b[0]))?"Ok":"Error");
    }
}
