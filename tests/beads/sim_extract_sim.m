% simulation
xx = randn(100,30)*0.1 + 0.5;
figure;plot(xx');ylim([0 1]);
xxMean = mean(xx,1);
figure;plot(xxMean);ylim([0 1]);


