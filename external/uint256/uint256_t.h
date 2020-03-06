// PUBLIC IMPORT HEADER
#ifndef _UINT256_H_
#define _UINT256_H_
#include "uint256_t_config.include"
#define UINT256_T_EXTERN _UINT256_T_IMPORT
#include "uint128_t.include"
#include "uint256_t.include"


namespace std
{

template <>
struct hash<uint256_t>
{
  uint256_t operator()(uint256_t key) const {
    return key;
  }
};

}

#endif
