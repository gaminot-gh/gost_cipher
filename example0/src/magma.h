/*! \file magma.h Прямая реализация алгоритма по ГОСТ Р 34.12-2015 п. 5 Алгоритм блочного шифрования с длиной блока n = 64 бит */
#ifndef MAGMA_H
#define MAGMA_H
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

/// Контрольные примеры
/// А.2 Алгоритм блочного шифрования с длиной блока n = 64 бит
void magmaCheck();


#ifdef __cplusplus
};
#endif

#endif
