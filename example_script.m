iF1 = 2;
iF2 = 55;

file1 = 'Data/SCAPE/mesh' + string(num2str(iF1,'%03d')) + '.off';
file2 = 'Data/SCAPE/mesh' + string(num2str(iF2,'%03d')) + '.off';

[P,tau,C,X,Y] = smoothshells(file1,file2);