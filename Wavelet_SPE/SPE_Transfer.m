% clc;
% New method (SPE, or said Q-static), PCA + dwt + SPE for transpose
% Apply PCA to data set of the spectrum with normalization of each freq bin, and get new "projection".
% using V which is previously build in SPE_Model to caculate SPE, and see if it is under threshold.

% Refer to
% Refer to
% https://www.exeley.com/exeley/journals/journal_of_theoretical_and_applied_mechanics/54/2/pdf/10.15632_jtam-pl.54.2.659.pdf
% Albert Hsu 2018/05/10

L_r = 8192;
N = 30; % 60 number of oberservation
FrameN = 1; %length(x1)/L;
L = L_r/FrameN;
dwtDepth = 5;
base = 'db5'; %bior6.8 rbio6.8


num_temp0 = 0;
num_temp1 = 0;
num_temp2 = 0;
num_temp3 = 0;
num_temp4 = 0;


Dir0 = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\Wavelet_SPE\DUT_pass\';
temp0=dir([Dir0,'*.wav']);
num_temp0=length(temp0);

Dir1 = 'C:\Matlab\work\PCA\Dajan 3000 pcs Test\Wavelet_SPE\DUT_fail\';
temp1=dir([Dir1,'*.wav']);
num_temp1=length(temp1);

% Dir2 = 'C:\Matlab\work\PCA\Payton\DUT_normal2\';
% temp2=dir([Dir2,'*.csv']);
% num_temp2=length(temp2);
% 
% Dir3 = 'C:\Matlab\work\PCA\Payton\DUT_normal3\';
% temp3=dir([Dir3,'*.csv']);
% num_temp3=length(temp3);
% 
% Dir4 = 'C:\Matlab\work\PCA\Payton\DUT_normal4\';
% temp4=dir([Dir4,'*.csv']);
% num_temp4=length(temp4);



N = num_temp0 + num_temp1 + num_temp2 + num_temp3 + num_temp4; 
x = zeros(L_r , N);

%%----------------discrete packle wavele transform-----------------
for i = 1:num_temp0
    filename=[Dir0,temp0(i).name];
    x(:, i) = wavread(filename, L_r);
end

for i = 1:num_temp1
    filename=[Dir1,temp1(i).name];
    x(:, num_temp0+i) = wavread(filename, L_r);
end

for i = 1:num_temp2
    filename=[Dir2,temp2(i).name];
    x(:, num_temp0 + num_temp1 +i) = wavread(filename, L_r);
end

for i = 1:num_temp3
    filename=[Dir3,temp3(i).name];
    x(:, num_temp0 + num_temp1 + num_temp2 +i) = wavread(filename, L_r);
end

for i = 1:num_temp4
    filename=[Dir4,temp4(i).name];
    x(:, num_temp0 + num_temp1 + num_temp2 + num_temp3+i) = wavread(filename, L_r);
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
offsetD1 = repmat(scaleMeanD1, 1, N);
offsetD2 = repmat(scaleMeanD2, 1, N);
offsetD3 = repmat(scaleMeanD3, 1, N);
if dwtDepth>=4
    offsetD4 = repmat(scaleMeanD4, 1, N);
end
if dwtDepth>=5
    offsetD5 = repmat(scaleMeanD5, 1, N);
end
offsetLA = repmat(scaleMeanLA, 1, N);

x_psD1 = x_psD1 - offsetD1;
x_psD2 = x_psD2 - offsetD2;
x_psD3 = x_psD3 - offsetD3;
if dwtDepth>=4
    x_psD4 = x_psD4 - offsetD4;
end
if dwtDepth>=5
    x_psD5 = x_psD5 - offsetD5;
end
x_psLA = x_psLA - offsetLA;





% normalized std for each spectrum data
for i = 1:length(x_psD1(:,1))
    x_psD1(i,:) = x_psD1(i,:)./ nor_xpsD1(i); 
end

for i = 1:length(x_psD2(:,1))
    x_psD2(i,:) = x_psD2(i,:)./ nor_xpsD2(i); 
end

for i = 1:length(x_psD3(:,1))
    x_psD3(i,:) = x_psD3(i,:)./ nor_xpsD3(i); 
end

if dwtDepth >= 4
    for i = 1:length(x_psD4(:,1))
        x_psD4(i,:) = x_psD4(i,:)./ nor_xpsD4(i);
    end
end

if dwtDepth >= 5
    for i = 1:length(x_psD5(:,1))
        x_psD5(i,:) = x_psD5(i,:)./ nor_xpsD5(i);
    end
end

for i = 1:length(x_psLA(:,1))
    x_psLA(i,:) = x_psLA(i,:)./ nor_xpsLA(i); 
end

%PCA transfer
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

%SPE caculation
pcsN = 80;
e = zeros(1,N);
I = eye(length(V(:,1)), length(V(1,:)));
P = V(:,1:pcsN)*V(:,1:pcsN)';
Q = PCAIn*(I-P);

for i =1:N
    e(i) = norm(Q(i,:));
end

% if e(i)>Qstatic, it is fail
for i = 1:num_temp0
    if e(i)>Qstatic
        temp0(i).name
    end
end

for i = 1:num_temp1
    if e(i+num_temp0) < Qstatic
        temp1(i).name
    end
end


% pcsN = 10;
% transD1 = (VD1(:,1:pcsN)'*x_psD1)';
% transD2 = (VD2(:,1:pcsN)'*x_psD2)';
% transD3 = (VD3(:,1:pcsN)'*x_psD3)';
% 
% if dwtDepth >= 4
%     transD4 = (VD4(:,1:pcsN)'*x_psD4)';
% end
% 
% if dwtDepth >= 5
%     transD5 = (VD5(:,1:pcsN)'*x_psD5)';
% end
% 
% transLA = (VLA(:,1:pcsN)'*x_psLA)';





