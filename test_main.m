function test_main
curwd = pwd;
mydir = fileparts(mfilename('fullpath'));
cd(fullfile(mydir,'test'));
try
    test_expmimv
    test_gexpm
    test_gexpmq
    test_kmatexp
catch me
    cd(curwd)
    rethrow(me);
end