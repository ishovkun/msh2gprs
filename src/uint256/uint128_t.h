// PUBLIC IMPORT HEADER
#ifndef _UINT128_H_
#define _UINT128_H_
#include "uint256_t_config.include"
#define UINT256_T_EXTERN _UINT256_T_IMPORT
#include "uint128_t.include"

namespace std
{

template <>
struct hash<uint128_t>
{
  uint128_t operator()(uint128_t key) const {
    return key;
  }
};

}


#endif
