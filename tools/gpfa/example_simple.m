load('mat_sample/sample_dat');
runIdx = 1;
method = 'gpfa';
xDim = 8;
kernSD = 30;
result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim,... 
                    'kernSDList', kernSD);