# uint256/CMakeLists.txt

ADD_LIBRARY(uint256
  # headers
  uint128_t.build
  uint128_t.include
  uint256_t_config.include
  uint256_t.include
  uint256_t.h
  uint128_t.h
  # implementation
  uint128_t.cpp
  uint256_t.cpp
)

SET_TARGET_PROPERTIES (
    uint256
    PROPERTIES LINKER_LANGUAGE CXX
)
