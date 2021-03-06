% clc;
% Dajan Test 
% Apply PCA to data set of the spectrum with normalization of each freq bin, and get new "projection".
% Select 30 wave files, include 10 pass, 10 fail and 10 unkonw to build the
% pca model.
% The result can be used as basics of classification.

% After model is build by the m-file. Run "group3_transfer" to classify
% wave file. 

% Refer to
% http://research.med.helsinki.fi/corefacilities/proteinchem/pca_introduction_basics.pdf
% and
% https://learnche.org/pid/latent-variable-modelling/principal-component-analysis/pca-example-analysis-of-spectral-data
% Albert Hsu  2017/11/01
% �ڭC �ҦV�L��

clear all;

L_r = 4096*30;
N = 30; % 60 number of oberservation
group = 3;
FrameN = 30; %length(x1)/L;
L = L_r/FrameN;

OKDir = 'C:\MATLAB701\work\Dajan 3000 pcs\Dajan 3000 pcs Test\ThreeGroup\passModel\';
temp0=dir([OKDir,'*.wav']);
num_temp0=length(temp0);



%files below are fail for sure
failDir = 'C:\MATLAB701\work\Dajan 3000 pcs\Dajan 3000 pcs Test\ThreeGroup\failModel\';
temp1=dir([failDir,'*.wav']);
num_temp1=length(temp1);


% x21 = wavread('1_fail.wav'); x21 = x21(1:L_r);
% x22 = wavread('26_fail.wav'); x22 = x22(1:L_r);
% x23 = wavread('3_fail.wav'); x23 = x23(1:L_r);
% x24 = wavread('27_fail.wav'); x24 = x24(1:L_r);
% x25 = wavread('5_fail.wav'); x25 = x25(1:L_r);
% x26 = wavread('6_fail.wav'); x26 = x26(1:L_r);
% x27 = wavread('7_fail.wav'); x27 = x27(1:L_r);
% x28 = wavread('8_fail.wav'); x28 = x28(1:L_r);
% x29 = wavread('30_fail.wav'); x29 = x29(1:L_r);
% x30 = wavread('11_fail.wav'); x30 = x30(1:L_r);
% x21 = wavread('1_fail.wav'); x21 = x21(1:L_r);
% x22 = wavread('2_fail.wav'); x22 = x22(1:L_r);
% x23 = wavread('3_fail.wav'); x23 = x23(1:L_r);
% x24 = wavread('4_fail.wav'); x24 = x24(1:L_r);
% x25 = wavread('5_fail.wav'); x25 = x25(1:L_r);
% x26 = wavread('6_fail.wav'); x26 = x26(1:L_r);
% x27 = wavread('7_fail.wav'); x27 = x27(1:L_r);
% x28 = wavread('8_fail.wav'); x28 = x28(1:L_r);
% x29 = wavread('26_fail.wav'); x29 = x29(1:L_r);
% x30 = wavread('11_fail.wav'); x30 = x30(1:L_r);
% x31 = wavread('27_fail.wav'); x31 = x31(1:L_r);
% x32 = wavread('13_fail.wav'); x32 = x32(1:L_r);
% x33 = wavread('14_fail.wav'); x33 = x33(1:L_r);
% x34 = wavread('17_fail.wav'); x34 = x34(1:L_r);
% x35 = wavread('19_fail.wav'); x35 = x35(1:L_r);
% x36 = wavread('21_fail.wav'); x36 = x36(1:L_r);
% x37 = wavread('22_fail.wav'); x37 = x37(1:L_r);
% x38 = wavread('23_fail.wav'); x38 = x38(1:L_r);
% x39 = wavread('24_fail.wav'); x39 = x39(1:L_r);
% x40 = wavread('25_fail.wav'); x40 = x40(1:L_r);

%files below are hard to determine
minorfailDir = 'C:\MATLAB701\work\Dajan 3000 pcs\Dajan 3000 pcs Test\ThreeGroup\minorModel\';
temp2=dir([minorfailDir,'*.wav']);
num_temp2=length(temp2);
% x42 = wavread('12_fail.wav'); x42 = x42(1:L_r);
% x49 = wavread('80_fail.wav'); x49 = x49(1:L_r);
% x52 = wavread('96_fail.wav'); x52 = x52(1:L_r);
% x53 = wavread('97_fail.wav'); x53 = x53(1:L_r);
% x54 = wavread('98_fail.wav'); x54 = x54(1:L_r);
% x55 = wavread('99_fail.wav'); x55 = x55(1:L_r);
% x57 = wavread('109_fail.wav'); x57 = x57(1:L_r);
% x58 = wavread('110_fail.wav'); x58 = x58(1:L_r);
% x59 = wavread('111_fail.wav'); x59 = x59(1:L_r);
% x60 = wavread('112_fail.wav'); x60 = x60(1:L_r);

N = num_temp0 + num_temp1 + num_temp2; 
x = zeros(L_r , N);

for i = 1:num_temp0
    filename=[OKDir,temp0(i).name];
    x(:, i) = wavread(filename, L_r);
end

for i = 1:num_temp1
    filename=[failDir,temp1(i).name];
    x(:, num_temp0+i) = wavread(filename, L_r);
end

for i = 1:num_temp2
    filename=[minorfailDir,temp2(i).name];
    x(:, num_temp0+ num_temp1 +i) = wavread(filename, L_r);
end
% x41 = wavread('9_fail.wav'); x41 = x41(1:L_r);
% x42 = wavread('12_fail.wav'); x42 = x42(1:L_r);
% x43 = wavread('15_fail.wav'); x43 = x43(1:L_r);
% x44 = wavread('18_fail.wav'); x44 = x44(1:L_r);
% x45 = wavread('20_fail.wav'); x45 = x45(1:L_r);
% x46 = wavread('28_fail.wav'); x46 = x46(1:L_r);
% x47 = wavread('31_fail.wav'); x47 = x47(1:L_r);
% x48 = wavread('32_fail.wav'); x48 = x48(1:L_r);
% x49 = wavread('80_fail.wav'); x49 = x49(1:L_r);
% x50 = wavread('89_fail.wav'); x50 = x50(1:L_r);
% x51 = wavread('92_fail.wav'); x51 = x51(1:L_r);
% x52 = wavread('96_fail.wav'); x52 = x52(1:L_r);
% x53 = wavread('97_fail.wav'); x53 = x53(1:L_r);
% x54 = wavread('98_fail.wav'); x54 = x54(1:L_r);
% x55 = wavread('99_fail.wav'); x55 = x55(1:L_r);
% x56 = wavread('108_fail.wav'); x56 = x56(1:L_r);
% x57 = wavread('109_fail.wav'); x57 = x57(1:L_r);
% x58 = wavread('110_fail.wav'); x58 = x58(1:L_r);
% x59 = wavread('111_fail.wav'); x59 = x59(1:L_r);
% x60 = wavread('112_fail.wav'); x60 = x60(1:L_r);


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


% normalized std for each spectrum data
for i = 1:L/2
    x_ps(i,:) = x_ps(i,:)./ sqrt(sum(x_ps(i,:).*x_ps(i,:))); 
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

%test: normalized projection on each pc
% for i = 1:pcsN
% trans(:,i) = trans(:,i)./S(i,i);
% end



hold on
for i = 1:num_temp0
    plot(trans(i,1:pcsN))
end

for i = 1:num_temp1
    plot(trans(num_temp0+i,1:pcsN), 'r')
end

for i = 1:num_temp2
    plot(trans(num_temp0+num_temp1+i, 1:pcsN), 'g')
end



%k mean cluster
cntPC = 2;
X = trans(:,1:cntPC);

[idx,c1] = kmeans(X,group)

% Find my cluster centers for group1 2 and 3
% Group1: pass. Group2:fail. 
% Group3: seems fail. Pick files 12 80 96 97 98 109 110 111 112. These files seems like fail 
c2 = zeros(group,cntPC);
c2(1,1:cntPC) = sum(trans(1:num_temp0, 1:cntPC)) / num_temp0;
c2(2,1:cntPC) = sum(trans(num_temp0+1 : num_temp0 + num_temp1, 1:cntPC)) / num_temp1;  
c2(3,1:cntPC) = sum(trans(num_temp0+num_temp1+1 : num_temp0 + num_temp1+num_temp2, 1:cntPC)) / num_temp2; 
 
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



