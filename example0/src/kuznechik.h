/// \file kuznechik.h  Прямая реализация алгоритма по ГОСТ Р 34.12-2015 п. 4 Алгоритм блочного шифрования с длиной блока n = 128 бит*/
#ifndef KUZNECHIK_H
#define KUZNECHIK_H
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{
#endif

/// Контрольные примеры
/// А.2 Алгоритм блочного шифрования с длиной блока n = 64 бит
void kuznechikCheck();


#ifdef __cplusplus
}
#endif





#endif
