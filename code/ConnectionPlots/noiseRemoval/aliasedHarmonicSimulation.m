fs=100000;
rsFs=500;

vec1 = (0:(1/fs):10)';
baseFreq = (60)*(2*pi);
% tmp1 = sin(vec1*baseFreq)+sin(vec1*baseFreq*2)/2+sin(vec1*baseFreq*3)/3+sin(vec1*baseFreq*4)/4+sin(vec1*baseFreq*5)/5+sin(vec1*baseFreq*6)/6+sin(vec1*baseFreq*7)/7;
tmp1 = sin(vec1*baseFreq*5)/5;

vec2=1:(fs/rsFs):numel(vec1);
vec3=(0:(1/fs):10)';
tmp2=tmp1(vec2);
figure();pwelch(tmp2,[],[],[],rsFs)
savepdf(gcf,'~/TaskAnalysis/lang/WN/noiseRemoval/aliasedHarmonics2.pdf')

figure();plot(tmp1)
hold on
scatter(vec2,tmp2)

