if ismac
mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims -I. gsqres_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -largeArrayDims -I. gexpm_hash_mex.cpp
else
mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims gsqres_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -I. -largeArrayDims gexpm_hash_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -Wall -std=c++0x" -largeArrayDims expm_svec_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims hpexpm_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims gexpmq_mex.cpp
end


