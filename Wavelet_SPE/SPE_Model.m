% clc;
% New method (SPE, or said Q-static), PCA + dwt + SPE for modeling 
% Apply PCA to data set of the spectrum with normalization of each freq bin, and get new "projection".
% Select "good" wave files only, include 10 pass to build the pca model.

% After model is build by the m-file. Run "transfer.m" to classify
% DUT. 

% Refer to
% https://www.exeley.com/exeley/journals/journal_of_theoretical_and_applied_mechanics/54/2/pdf/10.15632_jtam-pl.54.2.659.pdf
% Albert Hsu 2018/05/10

clear all;

L_r = 8192;
N = 30; % 60 number of oberservation
FrameN = 1; %length(x1)/L;
L = L_r/FrameN;
dwtDepth = 5;
base = 'db5'; %bior6.8 rbio6.8 db6 coif5


OKDir = 'C:\MATLAB701\work\Dajan 3000 pcs\Dajan 3000 pcs Test\Wavelet_SPE\Pass_forModel\';
%OKDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\passModel\';
temp0=dir([OKDir,'*.wav']);
num_temp0=length(temp0);



N = num_temp0; 
x = zeros(L_r , N);

%%----------------discrete packle wavele transform-----------------
for i = 1:num_temp0
    filename=[OKDir,temp0(i).name];
    x(:, i) = wavread(filename, L_r);
end



%h = hann(L/(2^dwtDepth));
h = hann(L);
for i = 1:N
    x(:,i) = x(:,i).*h;
end


if dwtDepth == 3
    [A,D] = dwt(x(:,1), base);
    L1 = length(A);
    dwtLev1 = zeros(L1,N*2);
    [AA,DA] = dwt(A(:,1), base);
    L2 = length(AA);
    dwtLev2 = zeros(L2,N*2*2);
    [AAA,DAA] = dwt(AA(:,1), base);
    L3 = length(AAA);
    dwtLev3 = zeros(L3, N*2*2*2);
    %x_ps_temp = zeros(L3, N*2^dwtDepth);
    x_psD1 = zeros(L1/2, N);
    x_psD2 = zeros(L2/2, N);
    x_psD3 = zeros(L3/2, N);
    x_psLA = zeros(L3/2, N);
    for i = 1:N                                                       %A  %D
        [dwtLev1(:,i),dwtLev1(:,N+i)] = dwt(x(:,i), base);   %dwtLev1(1:N N+1:2N): dwtLev1(:,1:N) is A. dwtLev1(:,N+1:2N) is D (D1) 
    end
    
    for i = 1:N*2                                                             %AA   AD       DA      DD
        [dwtLev2(:,i), dwtLev2(:,2*N+i)] = dwt(dwtLev1(:,i), base);  %dwtLev2(1:N , N+1:2N , 2N+1:3N, :)   dwtLev2(:,1:N) is AA. dwtLev2(:,2N+1:3N) is DA (D2)
    end
    
    for i = 1:N*4                                                            %AAA   AAD    ADA        ADD    DAA
        [dwtLev3(:,i), dwtLev3(:,4*N+i)] = dwt(dwtLev2(:,i), base); %dwtLev3(1:N,  N+1:2N, 2N+1:3N, 3N+1:4N, 4N+1:5N, ). dwtLev3(4N+1:5N) is D3
    end
    
    x_ps_temp = abs(fft(dwtLev1(:,N+1:2*N)));
    x_psD1 = x_ps_temp(1:L1/2,:); 
    x_ps_temp = abs(fft(dwtLev2(:,2*N+1:3*N)));
    x_psD2 = x_ps_temp(1:L2/2,:);
    x_ps_temp = abs(fft(dwtLev3(:,4*N+1:5*N)));
    x_psD3 = x_ps_temp(1:L3/2,:);
    x_ps_temp = abs(fft(dwtLev3(:,1:N)));
    x_psLA = x_ps_temp(1:L3/2,:);
end

if dwtDepth == 4
    [A,D] = dwt(x(:,1), base);
    L1 = length(A);
    dwtLev1 = zeros(L1,N*2);
    [AA,DA] = dwt(A(:,1), base);
    L2 = length(AA);
    dwtLev2 = zeros(L2,N*2*2);
    [AAA,DAA] = dwt(AA(:,1), base);
    L3 = length(AAA);
    dwtLev3 = zeros(L3, N*2*2*2);
    [AAAA, DAAA] = dwt(AAA(:,1), base);
    L4 = length(AAAA);
    dwtLev4 = zeros(L4, N*2*2*2*2);
    %x_ps_temp = zeros(L4, N*2^dwtDepth);
    x_psD1 = zeros(L1/2, N);
    x_psD2 = zeros(L2/2, N);
    x_psD3 = zeros(L3/2, N);
    x_psD4 = zeros(L4/2, N);
    x_psLA = zeros(L4/2, N);
    
    for i = 1:N
        [dwtLev1(:,i),dwtLev1(:,N+i)] = dwt(x(:,i), base);
    end
    
    for i = 1:N*2
        [dwtLev2(:,i), dwtLev2(:,2*N+i)] = dwt(dwtLev1(:,i), base);
    end
    
    for i = 1:N*4
        [dwtLev3(:,i), dwtLev3(:,4*N+i)] = dwt(dwtLev2(:,i), base);
    end
    
    for i = 1:N*8
        [dwtLev4(:,i), dwtLev4(:,8*N+i)] = dwt(dwtLev3(:,i), base);
    end
    
    x_ps_temp = abs(fft(dwtLev1(:,N+1:2*N)));
    x_psD1 = x_ps_temp(1:L1/2,:); 
    x_ps_temp = abs(fft(dwtLev2(:,2*N+1:3*N)));
    x_psD2 = x_ps_temp(1:L2/2,:);
    x_ps_temp = abs(fft(dwtLev3(:,4*N+1:5*N)));
    x_psD3 = x_ps_temp(1:L3/2,:);
    x_ps_temp = abs(fft(dwtLev4(:,8*N+1:9*N)));
    x_psD4 = x_ps_temp(1:L4/2,:);
    x_ps_temp = abs(fft(dwtLev4(:,1:N)));
    x_psLA = x_ps_temp(1:L4/2,:);
end

if dwtDepth == 5
    [A,D] = dwt(x(:,1), base);
    L1 = length(A);
    dwtLev1 = zeros(L1,N*2);
    [AA,DA] = dwt(A(:,1), base);
    L2 = length(AA);
    dwtLev2 = zeros(L2,N*2*2);
    [AAA,DAA] = dwt(AA(:,1), base);
    L3 = length(AAA);
    dwtLev3 = zeros(L3, N*2*2*2);
    [AAAA, DAAA] = dwt(AAA(:,1), base);
    L4 = length(AAAA);
    dwtLev4 = zeros(L4, N*2*2*2*2);
    [AAAAA, DAAAA] = dwt(AAAA(:,1), base);
    L5 = length(AAAAA);
    dwtLev5 = zeros(L5, N*2*2*2*2*2);
    x_psD1 = zeros(L1/2, N);
    x_psD2 = zeros(L2/2, N);
    x_psD3 = zeros(L3/2, N);
    x_psD4 = zeros(L4/2, N);
    x_psD5 = zeros(L5/2, N);
    x_psLA = zeros(L5/2, N);
    for i = 1:N
        [dwtLev1(:,i),dwtLev1(:,N+i)] = dwt(x(:,i), base);
    end
    
    for i = 1:N*2
        [dwtLev2(:,i), dwtLev2(:,2*N+i)] = dwt(dwtLev1(:,i), base);
    end
    
    for i = 1:N*4
        [dwtLev3(:,i), dwtLev3(:,4*N+i)] = dwt(dwtLev2(:,i), base);
    end
    
    for i = 1:N*8
        [dwtLev4(:,i), dwtLev4(:,8*N+i)] = dwt(dwtLev3(:,i), base);
    end
    
    for i = 1:N*16
        [dwtLev5(:,i), dwtLev5(:,16*N+i)] = dwt(dwtLev4(:,i), base);
    end
    
    x_ps_temp = abs(fft(dwtLev1(:,N+1:2*N)));
    x_psD1 = x_ps_temp(1:L1/2,:); 
    x_ps_temp = abs(fft(dwtLev2(:,2*N+1:3*N)));
    x_psD2 = x_ps_temp(1:L2/2,:);
    x_ps_temp = abs(fft(dwtLev3(:,4*N+1:5*N)));
    x_psD3 = x_ps_temp(1:L3/2,:);
    x_ps_temp = abs(fft(dwtLev4(:,8*N+1:9*N)));
    x_psD4 = x_ps_temp(1:L4/2,:);
    x_ps_temp = abs(fft(dwtLev5(:,16*N+1:17*N)));
    x_psD5 = x_ps_temp(1:L5/2,:);
    x_ps_temp = abs(fft(dwtLev5(:,1:N)));
    x_psLA = x_ps_temp(1:L5/2,:);
end

if dwtDepth == 6
    [A,D] = dwt(x(:,1), base);
    L1 = length(A);
    dwtLev1 = zeros(L1,N*2);
    [AA,DA] = dwt(A(:,1), base);
    L2 = length(AA);
    dwtLev2 = zeros(L2,N*2*2);
    [AAA,DAA] = dwt(AA(:,1), base);
    L3 = length(AAA);
    dwtLev3 = zeros(L3, N*2*2*2);
    [AAAA, DAAA] = dwt(AAA(:,1), base);
    L4 = length(AAAA);
    dwtLev4 = zeros(L4, N*2*2*2*2);
    [AAAAA, DAAAA] = dwt(AAAA(:,1), base);
    L5 = length(AAAAA);
    dwtLev5 = zeros(L5, N*2*2*2*2*2);
    [AAAAAA, DAAAAA] = dwt(AAAAA(:,1), base);
    L6 = length(AAAAAA);
    dwtLev6 = zeros(L6, N*2*2*2*2*2*2);
    for i = 1:N
        [dwtLev1(:,i),dwtLev1(:,N+i)] = dwt(x(:,i), base);
    end
    
    for i = 1:N*2
        [dwtLev2(:,i), dwtLev2(:,2*N+i)] = dwt(dwtLev1(:,i), base);
    end
    
    for i = 1:N*4
        [dwtLev3(:,i), dwtLev3(:,4*N+i)] = dwt(dwtLev2(:,i), base);
    end
    
    for i = 1:N*8
        [dwtLev4(:,i), dwtLev4(:,8*N+i)] = dwt(dwtLev3(:,i), base);
    end
    
    for i = 1:N*16
        [dwtLev5(:,i), dwtLev5(:,16*N+i)] = dwt(dwtLev4(:,i), base);
    end
    
    for i = 1:N*32
        [dwtLev6(:,i), dwtLev6(:,32*N+i)] = dwt(dwtLev5(:,i), base);
    end
    
%     PCA_feat = zeros(features*2^dwtDepth, N);
%     for i = 1:N
%         for j = 1:2^dwtDepth
%             PCA_feat(features*(j-1)+1,i) = var(dwtLev6(:,N*(j-1)+i));
%             PCA_feat(features*(j-1)+2,i) = skewness(dwtLev6(:,N*(j-1)+i));
%             PCA_feat(features*(j-1)+3,i) = kurtosis(dwtLev6(:,N*(j-1)+i));
%             PCA_feat(features*(j-1)+4,i) = sqrt(sum((dwtLev6(:,N*(j-1)+i)).*dwtLev6(:,N*(j-1)+i))/length(dwtLev6(:,1)));
%             PCA_feat(features*(j-1)+5,i) = mean(dwtLev6(:,N*(j-1)+i));
%             PCA_feat(features*(j-1)+6,i) = (max(dwtLev6(:,N*(j-1)+i))-min(dwtLev6(:,N*(j-1)+i))) / sqrt(mean([dwtLev6(:,N*(j-1)+i)].^2)); %peak 2 rms. aka crest factor
%         end
%     end
end

x_psD1 = 20*log10(x_psD1);
x_psD2 = 20*log10(x_psD2);
x_psD3 = 20*log10(x_psD3);
if dwtDepth >= 4
    x_psD4 = 20*log10(x_psD4);
end

if dwtDepth >= 5
    x_psD5 = 20*log10(x_psD5);
end
x_psLA = 20*log10(x_psLA);

%normalize and get rid off offset of power spectrum
scaleMeanD1 = zeros(length(x_psD1(:,1)),1);
for i = 1:length(x_psD1(:,1))
    tempXs = x_psD1(i,:);
    scaleMeanD1(i) = mean(tempXs);
    x_psD1(i,:) = tempXs - mean(tempXs);
end

scaleMeanD2 = zeros(length(x_psD2(:,1)),1);
for i = 1:length(x_psD2(:,1))
    tempXs = x_psD2(i,:);
    scaleMeanD2(i) = mean(tempXs);
    x_psD2(i,:) = tempXs - mean(tempXs);
end

scaleMeanD3 = zeros(length(x_psD3(:,1)),1);
for i = 1:length(x_psD3(:,1))
    tempXs = x_psD3(i,:);
    scaleMeanD3(i) = mean(tempXs);
    x_psD3(i,:) = tempXs - mean(tempXs);
end

if dwtDepth >= 4
    scaleMeanD4 = zeros(length(x_psD4(:,1)),1);
    for i = 1:length(x_psD4(:,1))
        tempXs = x_psD4(i,:);
        scaleMeanD4(i) = mean(tempXs);
        x_psD4(i,:) = tempXs - mean(tempXs);
    end
end

if dwtDepth >= 5
    scaleMeanD5 = zeros(length(x_psD5(:,1)),1);
    for i = 1:length(x_psD5(:,1))
        tempXs = x_psD5(i,:);
        scaleMeanD5(i) = mean(tempXs);
        x_psD5(i,:) = tempXs - mean(tempXs);
    end
end

scaleMeanLA = zeros(length(x_psLA(:,1)),1);
for i = 1:length(x_psLA(:,1))
    tempXs = x_psLA(i,:);
    scaleMeanLA(i) = mean(tempXs);
    x_psLA(i,:) = tempXs - mean(tempXs);
end


% normalized std for each spectrum data
nor_xpsD1 = zeros(length(x_psD1(:,1)),1);
for i = 1:length(x_psD1(:,1))
    nor_xpsD1(i) = sqrt(sum(x_psD1(i,:).*x_psD1(i,:))/N);
    x_psD1(i,:) = x_psD1(i,:)./ nor_xpsD1(i); 
end

nor_xpsD2 = zeros(length(x_psD2(:,1)),1);
for i = 1:length(x_psD2(:,1))
    nor_xpsD2(i) = sqrt(sum(x_psD2(i,:).*x_psD2(i,:))/N);
    x_psD2(i,:) = x_psD2(i,:)./ nor_xpsD2(i); 
end

nor_xpsD3 = zeros(length(x_psD3(:,1)),1);
for i = 1:length(x_psD3(:,1))
    nor_xpsD3(i) = sqrt(sum(x_psD3(i,:).*x_psD3(i,:))/N);
    x_psD3(i,:) = x_psD3(i,:)./ nor_xpsD3(i); 
end

if dwtDepth >= 4
    nor_xpsD4 = zeros(length(x_psD4(:,1)),1);
    for i = 1:length(x_psD4(:,1))
        nor_xpsD4(i) = sqrt(sum(x_psD4(i,:).*x_psD4(i,:))/N);
        x_psD4(i,:) = x_psD4(i,:)./ nor_xpsD4(i);
    end
end

if dwtDepth >= 5
    nor_xpsD5 = zeros(length(x_psD5(:,1)),1);
    for i = 1:length(x_psD5(:,1))
        nor_xpsD5(i) = sqrt(sum(x_psD5(i,:).*x_psD5(i,:))/N);
        x_psD5(i,:) = x_psD5(i,:)./ nor_xpsD5(i);
    end
end

nor_xpsLA = zeros(length(x_psLA(:,1)),1);
for i = 1:length(x_psLA(:,1))
    nor_xpsLA(i) = sqrt(sum(x_psLA(i,:).*x_psLA(i,:))/N); 
    x_psLA(i,:) = x_psLA(i,:)./ nor_xpsLA(i); 
end

%PCA modeling
if dwtDepth == 3
    PCAIn = [x_psD1, x_psD2, x_psD3, x_psLA];
end

if dwtDepth == 4
    PCAIn = [x_psD1, x_psD2, x_psD3, x_psD4, x_psLA];
end

if dwtDepth == 5
     PCAIn = [x_psD1', x_psD2', x_psD3', x_psD4', x_psD5', x_psLA'];
end

if dwtDepth == 6
    PCAIn = [x_psD1, x_psD2, x_psD3, x_psD4, x_psD5, x_psD6, x_psLA];
end
%pcsN = 10;

[U , S, V] = svd(PCAIn, 0);

%Threshold Q-statistic:
pcsN = 80;
thi1 = 0;
thi2 = 0;
thi3 = 0;
for j = pcsN+1:N
    thi1 = thi1 + S(j,j)^(1/2);
    thi2 = thi2 + S(j,j)^(2/2);
    thi3 = thi3 + S(j,j)^(3/2);
end
h0 = 1- 2*thi1*thi3/(3*thi2^2);
criticalValue  = 1.96; %95%
Qstatic = thi1*(criticalValue*h0*sqrt(2*thi2)/thi1 + 1 + thi2*h0*(h0-1)/thi1^2)
%SPE caculation
% test = V(:,1:N)*V(:,1:N)';
% e = eye(4116, 4116) - test;
% Qstatic = PCAIn(1,:)*e;

