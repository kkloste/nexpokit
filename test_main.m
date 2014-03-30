function test_main
curwd = pwd;
mydir = fileparts(mfilename('fullpath'));
cd(fullfile(mydir,'test'));
try
    test_expm_svec
    test_gexpm_hash
    test_gsqres
    test_kmatexp
catch me
    cd(curwd)
    rethrow(me);
end