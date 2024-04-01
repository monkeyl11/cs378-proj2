
% taken from mathworks help page
fileID = fopen('output.txt','r');
formatSpec = '%d %d %f';
sizeA = [3 Inf];
A = fscanf(fileID,formatSpec,sizeA);
fclose(fileID);
A = A';
T = array2table(A);
T.Properties.VariableNames(1:3) = {'dim','ours','ref'};
figure
plot(T.dim, T.ours, T.dim, T.ref);
title('Performance Plot');
xlabel('Dimension (nxn special matrix)');
ylabel('GFLOPS');
legend({'our implementation','reference implementation'},'Location','southwest')

print( 'plot', '-dpng' );
