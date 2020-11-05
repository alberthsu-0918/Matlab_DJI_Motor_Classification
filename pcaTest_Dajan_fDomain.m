% clc;
% Dajan Test 2.1. Finding 10 peaks of wave's spectrum. 
% Apply PCA to data set of the peak and get new "projection".
% The result can be used as basics of classification.
% Refer to
% http://research.med.helsinki.fi/corefacilities/proteinchem/pca_introduction_basics.pdf
% and
% https://learnche.org/pid/latent-variable-modelling/principal-component-analysis/pca-example-analysis-of-spectral-data
% Albert Hsu  2017/01/26
% ¼Ú­C ©Ò¦VµL¼Ä

clear all;

L_r = 128000;
L = 4096; %load 4096 samples of wavefile
N = 40; %number of oberservation


x1 = audioread('1_pass.wav');
x2 = audioread('2_pass.wav');
x3 = audioread('3_pass.wav'); x3 = x3(1:128000);
x4 = audioread('4_pass.wav');
x5 = audioread('5_pass.wav');
x6 = audioread('6_pass.wav');
x7 = audioread('7_pass.wav');
x8 = audioread('8_pass.wav');
x9 = audioread('9_pass.wav');
x10 = audioread('10_pass.wav');
x11 = audioread('11_pass.wav');
x12 = audioread('12_pass.wav');
x13 = audioread('13_pass.wav');
x14 = audioread('14_pass.wav');
x15 = audioread('15_pass.wav');
x16 = audioread('16_pass.wav');
x17 = audioread('17_pass.wav');
x18 = audioread('18_pass.wav');
x19 = audioread('19_pass.wav');
x20 = audioread('20_pass.wav');
 
% x21 = audioread('1_fail.wav'); x21 = x21(1:L_r);
% x22 = audioread('2_fail.wav'); x22 = x22(1:L_r);
% x23 = audioread('3_fail.wav'); x23 = x23(1:L_r);
% x24 = audioread('4_fail.wav'); x24 = x24(1:L_r);
% x25 = audioread('5_fail.wav'); x25 = x25(1:L_r);
% x26 = audioread('6_fail.wav'); x26 = x26(1:L_r);
% x27 = audioread('7_fail.wav'); x27 = x27(1:L_r);
% x28 = audioread('8_fail.wav'); x28 = x28(1:L_r);
% x29 = audioread('9_fail.wav'); x29 = x29(1:L_r);
% x30 = audioread('21_fail.wav'); x30 = x30(1:L_r);
% x31 = audioread('11_fail.wav'); x31 = x31(1:L_r);
% x32 = audioread('12_fail.wav'); x32 = x32(1:L_r);
% x33 = audioread('13_fail.wav'); x33 = x33(1:L_r);
% x34 = audioread('14_fail.wav'); x34 = x34(1:L_r);
% x35 = audioread('15_fail.wav'); x35 = x35(1:L_r);
% x36 = audioread('16_fail.wav'); x36 = x36(1:L_r);
% x37 = audioread('17_fail.wav'); x37 = x37(1:L_r);
% x38 = audioread('18_fail.wav'); x38 = x38(1:L_r);
% x39 = audioread('19_fail.wav'); x39 = x39(1:L_r);
% x40 = audioread('20_fail.wav'); x40 = x40(1:L_r);

% x21 = audioread('79_fail.wav'); x21 = x21(1:L_r);
% x22 = audioread('80_fail.wav'); x22 = x22(1:L_r);
% x23 = audioread('81_fail.wav'); x23 = x23(1:L_r);
% x24 = audioread('82_fail.wav'); x24 = x24(1:L_r);
% x25 = audioread('83_fail.wav'); x25 = x25(1:L_r);
% x26 = audioread('84_fail.wav'); x26 = x26(1:L_r);
% x27 = audioread('85_fail.wav'); x27 = x27(1:L_r);
% x28 = audioread('86_fail.wav'); x28 = x28(1:L_r);
% x29 = audioread('87_fail.wav'); x29 = x29(1:L_r);
% x30 = audioread('88_fail.wav'); x30 = x30(1:L_r);
% x31 = audioread('89_fail.wav'); x31 = x31(1:L_r);
% x32 = audioread('90_fail.wav'); x32 = x32(1:L_r);
% x33 = audioread('91_fail.wav'); x33 = x33(1:L_r);
% x34 = audioread('92_fail.wav'); x34 = x34(1:L_r);
% x35 = audioread('93_fail.wav'); x35 = x35(1:L_r);
% x36 = audioread('95_fail.wav'); x36 = x36(1:L_r);
% x37 = audioread('96_fail.wav'); x37 = x37(1:L_r);
% x38 = audioread('97_fail.wav'); x38 = x38(1:L_r);
% x39 = audioread('98_fail.wav'); x39 = x39(1:L_r);
% x40 = audioread('99_fail.wav'); x40 = x40(1:L_r);

%file below are defitinite fail
x21 = audioread('1_fail.wav'); x21 = x21(1:L_r);
x22 = audioread('29_fail.wav'); x22 = x22(1:L_r);
x23 = audioread('3_fail.wav'); x23 = x23(1:L_r);
x24 = audioread('30_fail.wav'); x24 = x24(1:L_r);
x25 = audioread('5_fail.wav'); x25 = x25(1:L_r);
x26 = audioread('6_fail.wav'); x26 = x26(1:L_r);
x27 = audioread('7_fail.wav'); x27 = x27(1:L_r);
x28 = audioread('8_fail.wav'); x28 = x28(1:L_r);
x29 = audioread('26_fail.wav'); x29 = x29(1:L_r);
x30 = audioread('11_fail.wav'); x30 = x30(1:L_r);
x31 = audioread('27_fail.wav'); x31 = x31(1:L_r);
x32 = audioread('13_fail.wav'); x32 = x32(1:L_r);
x33 = audioread('68_fail.wav'); x33 = x33(1:L_r);
x34 = audioread('17_fail.wav'); x34 = x34(1:L_r);
x35 = audioread('19_fail.wav'); x35 = x35(1:L_r);
x36 = audioread('21_fail.wav'); x36 = x36(1:L_r);
x37 = audioread('72_fail.wav'); x37 = x37(1:L_r);
x38 = audioread('23_fail.wav'); x38 = x38(1:L_r);
x39 = audioread('24_fail.wav'); x39 = x39(1:L_r);
x40 = audioread('25_fail.wav'); x40 = x40(1:L_r);

%files below are hard to determine
% x21 = audioread('9_fail.wav'); x21 = x21(1:L_r);
% x22 = audioread('12_fail.wav'); x22 = x22(1:L_r);
% x23 = audioread('15_fail.wav'); x23 = x23(1:L_r);
% x24 = audioread('18_fail.wav'); x24 = x24(1:L_r);
% x25 = audioread('20_fail.wav'); x25 = x25(1:L_r);
% x26 = audioread('28_fail.wav'); x26 = x26(1:L_r);
% x27 = audioread('31_fail.wav'); x27 = x27(1:L_r);
% x28 = audioread('32_fail.wav'); x28 = x28(1:L_r);
% x29 = audioread('80_fail.wav'); x29 = x29(1:L_r);
% x30 = audioread('89_fail.wav'); x30 = x30(1:L_r);
% x31 = audioread('92_fail.wav'); x31 = x31(1:L_r);
% x32 = audioread('96_fail.wav'); x32 = x32(1:L_r);
% x33 = audioread('97_fail.wav'); x33 = x33(1:L_r);
% x34 = audioread('98_fail.wav'); x34 = x34(1:L_r);
% x35 = audioread('99_fail.wav'); x35 = x35(1:L_r);
% x36 = audioread('108_fail.wav'); x36 = x36(1:L_r);
% x37 = audioread('109_fail.wav'); x37 = x37(1:L_r);
% x38 = audioread('110_fail.wav'); x38 = x38(1:L_r);
% x39 = audioread('111_fail.wav'); x39 = x39(1:L_r);
% x40 = audioread('112_fail.wav'); x40 = x40(1:L_r);

x = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23 x24 x25 x26 x27 x28 x29 x30 x31 x32 x33 x34 x35 x36 x37 x38 x39 x40];

% normalize energy of each wav file
% for i = 1:N
%     x(:,i) = x(:,i)./sqrt(sum(x(:,i).*x(:,i)));
% end
FrameN = 3; %length(x1)/L;


h = hann(L);
%h = 1;
x_ps = zeros(L/2, N);
x_h = zeros(L, N);

% average FFT spectrum
for j = 1:N
    for i = 1:FrameN
        x_h(:,j) = x((i-1)*L+1 : L*i, j);
        x_h(:,j) = x_h(:,j).*h;
        x_ps_temp = abs(fft(x_h(:,j)));
        x_ps_temp = x_ps_temp(1:L/2);
        x_ps(:,j) = x_ps(:,j) + x_ps_temp; 
    end
end


%normalize and get rid off offset of power spectrum

x_ps = 20*log10(x_ps);

%get rid of offese
scaleMean = zeros(L/2,1);
for i = 1:L/2
    tempXs = x_ps(i,:);
    scaleMean(i) = mean(tempXs);
    x_ps(i,:) = tempXs - mean(tempXs);
end

%normalized
for i = 1:L/2
    x_ps(i,:) = x_ps(i,:)./sqrt( sum(x_ps(i,:).*x_ps(i,:)) );
end



%% Find 10 peaks as characteristics 
%% useless
% psSum = zeros(L/2,1);
% for i = 1:N
%     psSum = psSum + x_ps(:,N);
% end
% 
% loc = peakfinder(psSum);
% x_ps_peak = [x_ps(loc,1) x_ps(loc,2) x_ps(loc,3) x_ps(loc,4) x_ps(loc,5) x_ps(loc,6) x_ps(loc,7) x_ps(loc,8) x_ps(loc,9) x_ps(loc,10) ];


%PCA Analysis
%[pcs, trans, evs] = princomp(x_ps');
%[pcs, trans, evs] = princomp(x_ps_peak');
%[trans,pcs,evs] = pca1(x_ps_peak);
% trans = real(trans);
% Below are SVD transpose
pcsN = 10;
[U , S, V] = svd(x_ps');
%trans = (pcs'*x_ps)';
trans = (V'*x_ps)';


plot(trans(1,1:pcsN))
hold on
plot(trans(2,1:pcsN))
plot(trans(3,1:pcsN))
plot(trans(4,1:pcsN))
plot(trans(5,1:pcsN))
plot(trans(6,1:pcsN))
plot(trans(7,1:pcsN))
plot(trans(8,1:pcsN))
plot(trans(9,1:pcsN))
plot(trans(10,1:pcsN))
plot(trans(11,1:pcsN))
plot(trans(12,1:pcsN))
plot(trans(13,1:pcsN))
plot(trans(14,1:pcsN))
plot(trans(15,1:pcsN))
plot(trans(16,1:pcsN))
plot(trans(17,1:pcsN))
plot(trans(18,1:pcsN))
plot(trans(19,1:pcsN))
plot(trans(20,1:pcsN))
plot(trans(21,1:pcsN), 'r')
plot(trans(22,1:pcsN), 'r')
plot(trans(23,1:pcsN), 'r')
plot(trans(24,1:pcsN), 'r')
plot(trans(25,1:pcsN), 'r')
plot(trans(26,1:pcsN), 'r')
plot(trans(27,1:pcsN), 'r')
plot(trans(28,1:pcsN), 'r')
plot(trans(29,1:pcsN), 'r')
plot(trans(30,1:pcsN), 'r')
plot(trans(31,1:pcsN), 'r')
plot(trans(32,1:pcsN), 'r')
plot(trans(33,1:pcsN), 'r')
plot(trans(34,1:pcsN), 'r')
plot(trans(35,1:pcsN), 'r')
plot(trans(36,1:pcsN), 'r')
plot(trans(37,1:pcsN), 'r')
plot(trans(38,1:pcsN), 'r')
plot(trans(39,1:pcsN), 'r')
plot(trans(40,1:pcsN), 'r')

%k mean cluster
X = trans(:,1:3);
[idx, c1] = kmeans(X,2);
 %c1 = [-3.7662   -0.3641; 5.6494    0.5462];
% t = 1:51200/L:51200/2;

% %% SVD Test
% for i = 1:5
%     meanV(:,i) = mean(x(:,i))*ones(L,1);  
% end
% xz = x - meanV;
% %[U , S, V] = svds(xz,SensorN);
% [U , S, V] = svds(xz,1);
% xr = U*S*V';
% xr_fs = abs(fft(xr(:,1)));
% xr_fs = xr_fs(1:length(xr_fs)/2);



