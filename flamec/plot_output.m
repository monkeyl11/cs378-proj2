
system("make driver");
system("driver.x");


% taken from mathworks help page
fileID = fopen('output.txt','r');
formatSpec = '%d %f %f';

fscanf(fileID,'')

while ~feof(fileID)
    name = fscanf(fileID,'%s', [1]);
    len = fscanf(fileID,'%d', [1]);

    sizeA = [3 len];
    A = fscanf(fileID,formatSpec,sizeA);
    A = A';

    T = array2table(A);
    T.Properties.VariableNames(1:3) = {'dim','ours','ref'};
    figure
    plot(T.dim, T.ours, T.dim, T.ref);
    title([name, 'performance plot']);
    xlabel('Dimension (nxn special matrix)');
    ylabel('Time taken to finish (seconds)');
    legend({'our implementation','reference implementation'},'Location','southwest')

    print( [name,'_plot'], '-dpng' );
end

fclose(fileID);