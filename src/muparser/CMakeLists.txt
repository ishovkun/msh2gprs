# muparser/CMakeLists.txt


ADD_LIBRARY(muparser
  muParserBase.cpp
  muParserBase.h
  muParserBytecode.cpp
  muParserBytecode.h
  muParserCallback.cpp
  muParserCallback.h
  muParser.cpp
  muParserDef.h
  muParserDLL.cpp
  muParserDLL.h
  muParserError.cpp
  muParserError.h
  muParserFixes.h
  muParser.h
  muParserInt.cpp
  muParserInt.h
  muParserStack.h
  muParserTemplateMagic.h
  muParserTest.cpp
  muParserTest.h
  muParserToken.h
  muParserTokenReader.cpp
  muParserTokenReader.h
  )

SET_TARGET_PROPERTIES (
    muparser
    PROPERTIES LINKER_LANGUAGE CXX
)
