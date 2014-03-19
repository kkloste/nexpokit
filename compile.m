%mex -g CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x -D_GLIBCXX_DEBUG" -I. -largeArrayDims gsqres_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims gsqres_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -largeArrayDims gexpm_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims gexpmq_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims gexpm_hash_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -largeArrayDims expm_svec_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims hpexpm_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims gexpmq_mex.cpp


