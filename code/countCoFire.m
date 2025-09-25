function N = countCoFire(datA, datB, times, coFwin)


datA = find(datA(times));
datB = find(datB(times));

if isempty(datA) || isempty(datB)
    N = 0;
else
    datAB = datA - datB';
    N = sum(datAB(:) >= -coFwin & datAB(:) <= coFwin);

end




return