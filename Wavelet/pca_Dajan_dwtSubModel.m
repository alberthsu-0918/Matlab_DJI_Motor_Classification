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
% Albert Hsu  2018/04/20
% This is wavelet transform test

clear all;

L_r = 4096;
N = 30; % 60 number of oberservation
group = 2;
FrameN = 1; %length(x1)/L;
L = L_r/FrameN;
dwtDepth = 5;
base = 'db12'; %bior6.8 rbio6.8 db6 coif5
features = 6;

OKDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\Wavelet\passModel\';
%OKDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\passModel\';
temp0=dir([OKDir,'*.wav']);
num_temp0=length(temp0);


%files below are hard to determine
minorfailDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\Wavelet\minorModel\';
%minorfailDir = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\ThreeGroup\minorModel\';
temp1=dir([minorfailDir,'*.wav']);
num_temp1=length(temp1);


N = num_temp0 + num_temp1; 
x = zeros(L_r , N);

%%----------------discrete packle wavele transform-----------------
for i = 1:num_temp0
    filename=[OKDir,temp0(i).name];
    x(:, i) = wavread(filename, L_r);
end

for i = 1:num_temp1
    filename=[minorfailDir,temp1(i).name];
    x(:, num_temp0+i) = wavread(filename, L_r);
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
    
    PCA_feat = zeros(features*2^dwtDepth, N);
    for i = 1:N
        for j = 1:2^dwtDepth
            PCA_feat(features*(j-1)+1,i) = var(dwtLev6(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+2,i) = skewness(dwtLev6(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+3,i) = kurtosis(dwtLev6(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+4,i) = sqrt(sum((dwtLev6(:,N*(j-1)+i)).*dwtLev6(:,N*(j-1)+i))/length(dwtLev6(:,1)));
            PCA_feat(features*(j-1)+5,i) = mean(dwtLev6(:,N*(j-1)+i));
            PCA_feat(features*(j-1)+6,i) = (max(dwtLev6(:,N*(j-1)+i))-min(dwtLev6(:,N*(j-1)+i))) / sqrt(mean([dwtLev6(:,N*(j-1)+i)].^2)); %peak 2 rms. aka crest factor
        end
    end
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
pcsN = 10;
[U , S, VD1] = svd(x_psD1', 0);
transD1 = (VD1(:,1:pcsN)'*x_psD1)';

[U , S, VD2] = svd(x_psD2', 0);
transD2 = (VD2(:,1:pcsN)'*x_psD2)';

[U , S, VD3] = svd(x_psD3', 0);
transD3 = (VD3(:,1:pcsN)'*x_psD3)';

if dwtDepth >= 4
    [U , S, VD4] = svd(x_psD4', 0);
    transD4 = (VD4(:,1:pcsN)'*x_psD4)';
end

if dwtDepth >= 5
    [U , S, VD5] = svd(x_psD5', 0);
    transD5 = (VD5(:,1:pcsN)'*x_psD5)';
end

[U , S, VLA] = svd(x_psLA', 0);
transLA = (VLA(:,1:pcsN)'*x_psLA)';

good = zeros(num_temp0,dwtDepth+1);
bad  = zeros(num_temp1,dwtDepth+1);


if dwtDepth == 3
    for i = 1:num_temp0
        good(i,1) = transD1(i,1);
        good(i,2) = transD2(i,1);
        good(i,3) = transD3(i,1);
        good(i,4) = transLA(i,1);
    end

    for i = 1:num_temp1
        bad(i,1) = transD1(num_temp0+i,1);
        bad(i,2) = transD2(num_temp0+i,1);
        bad(i,3) = transD3(num_temp0+i,1);
        bad(i,4) = transLA(num_temp0+i,1);
    end
end

if dwtDepth == 4
    for i = 1:num_temp0
        good(i,1) = transD1(i,1);
        good(i,2) = transD2(i,1);
        good(i,3) = transD3(i,1);
        good(i,4) = transD4(i,1);
        good(i,5) = transLA(i,1);
    end

    for i = 1:num_temp1
        bad(i,1) = transD1(num_temp0+i,1);
        bad(i,2) = transD2(num_temp0+i,1);
        bad(i,3) = transD3(num_temp0+i,1);
        bad(i,4) = transD4(num_temp0+i,1);
        bad(i,5) = transLA(num_temp0+i,1);
    end
end

if dwtDepth == 5
    for i = 1:num_temp0
        good(i,1) = transD1(i,1);
        good(i,2) = transD2(i,1);
        good(i,3) = transD3(i,1);
        good(i,4) = transD4(i,1);
        good(i,5) = transD5(i,1);
        good(i,6) = transLA(i,1);
    end

    for i = 1:num_temp1
        bad(i,1) = transD1(num_temp0+i,1);
        bad(i,2) = transD2(num_temp0+i,1);
        bad(i,3) = transD3(num_temp0+i,1);
        bad(i,4) = transD4(num_temp0+i,1);
        bad(i,5) = transD5(num_temp0+i,1);
        bad(i,6) = transLA(num_temp0+i,1);
    end
end

%PCA Analysis
%[pcs, trans, evs] = princomp(x_ps');
%[pcs, trans, evs] = princomp(x_ps_peak');
%[trans,pcs,evs] = pca1(x_ps_peak);
% trans = real(trans);
% % -----------------Below are SVD transpose in freq
% pcsN = 10;
% [U , S, V] = svd(x_ps');
% %trans = (pcs'*x_ps)';
% trans = (V'*x_ps)';

%test: normalized projection on each pc
% for i = 1:pcsN
% trans(:,i) = trans(:,i)./S(i,i);
% end


figure
hold on
for i = 1:num_temp0
    plot(good(i,:))
end

for i = 1:num_temp1
    plot(bad(i,:), 'r')
end




%k mean cluster
%cntPC = 2^dwtDepth;
cntPC = dwtDepth+1 % 
X = zeros(num_temp0 + num_temp1,cntPC);
X(1:num_temp0,1:cntPC) = good;
X(num_temp0+1:num_temp0 + num_temp1,1:cntPC) = bad;


[idx,c1] = kmeans(X,group)

% Find my cluster centers for group1 2 and 3
% Group1: pass. Group2:fail. 
% Group3: seems fail. Pick files 12 80 96 97 98 109 110 111 112. These files seems like fail 
c2 = zeros(group,cntPC);
c2(1,1:cntPC) = sum(good(:,1:1:cntPC)) / num_temp0;
c2(2,1:cntPC) = sum(bad(:,1:1:cntPC)) / num_temp1;  
 


% New standard: using threshold (top & bottom)
stdT = 4; %3
limit1T = mean(good(:,1)) + stdT*sqrt(var(good(:,1)));
limit1B = mean(good(:,1)) - stdT*sqrt(var(good(:,1)));
limit2T = mean(good(:,2)) + stdT*sqrt(var(good(:,2)));
limit2B = mean(good(:,2)) - stdT*sqrt(var(good(:,2)));
limit3T = mean(good(:,3)) + stdT*sqrt(var(good(:,3)));
limit3B = mean(good(:,3)) - stdT*sqrt(var(good(:,3)));
limit4T = mean(good(:,4)) + stdT*sqrt(var(good(:,4)));
limit4B = mean(good(:,4)) - stdT*sqrt(var(good(:,4)));

if dwtDepth >=4
    limit5T = mean(good(:,5)) + stdT*sqrt(var(good(:,5)));
    limit5B = mean(good(:,5)) - stdT*sqrt(var(good(:,5)));
end

if dwtDepth >=5
    limit6T = mean(good(:,6)) + stdT*sqrt(var(good(:,6)));
    limit6B = mean(good(:,6)) - stdT*sqrt(var(good(:,6)));
end


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



