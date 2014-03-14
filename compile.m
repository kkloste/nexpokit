
if ismac
mex -O -largeArrayDims gsqres_mex.cpp
else
%mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims gsqres_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims gexpmq_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims gexpm_ord_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims gexpm_mex.cpp
mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims gexpm_hash_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims expm_svec_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -std=c++0x" -largeArrayDims expm_tay_mex.cpp
%mex -O CXXFLAGS="\$CXXFLAGS -Wall" -largeArrayDims hpexpm_mex.cpp
end