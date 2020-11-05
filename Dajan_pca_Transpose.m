% clc;
% Dajan Test.  
% After run 'pcaTest_Dajan_fDomain.m' to get pcs, the transpose matrix
% Run this m file to see if the pcs can project other wavefile to a better
% 3d axis.
% Further clustering is needed
%clear all;
L_r = 128000;
L = 4096; %load 4096 samples of wavefile
N = 40; %number of oberservation


x1 = audioread('21_pass.wav'); x1 = x1(1:L_r);
x2 = audioread('22_pass.wav'); x2 = x2(1:L_r);
x3 = audioread('23_pass.wav'); x3 = x3(1:L_r);
x4 = audioread('24_pass.wav'); x4 = x4(1:L_r);
x5 = audioread('25_pass.wav'); x5 = x5(1:L_r);
x6 = audioread('26_pass.wav'); x6 = x6(1:L_r);
x7 = audioread('27_pass.wav'); x7 = x7(1:L_r);
x8 = audioread('28_pass.wav'); x8 = x8(1:L_r);
x9 = audioread('29_pass.wav'); x9 = x9(1:L_r);
x10 = audioread('30_pass.wav'); x10 = x10(1:L_r);
x11 = audioread('31_pass.wav'); x11 = x11(1:L_r);
x12 = audioread('32_pass.wav'); x12 = x12(1:L_r);
x13 = audioread('33_pass.wav'); x13 = x13(1:L_r);
x14 = audioread('34_pass.wav'); x14 = x14(1:L_r);
x15 = audioread('35_pass.wav'); x15 = x15(1:L_r);
x16 = audioread('36_pass.wav'); x16 = x16(1:L_r);
x17 = audioread('37_pass.wav'); x17 = x17(1:L_r);
x18 = audioread('38_pass.wav'); x18 = x18(1:L_r);
x19 = audioread('39_pass.wav'); x19 = x19(1:L_r);
x20 = audioread('40_pass.wav'); x20 = x20(1:L_r);

% x1 = audioread('41_pass.wav'); x1 = x1(1:L_r);
% x2 = audioread('42_pass.wav'); x2 = x2(1:L_r);
% x3 = audioread('43_pass.wav'); x3 = x3(1:L_r);
% x4 = audioread('44_pass.wav'); x4 = x4(1:L_r);
% x5 = audioread('45_pass.wav'); x5 = x5(1:L_r);
% x6 = audioread('46_pass.wav'); x6 = x6(1:L_r);
% x7 = audioread('47_pass.wav'); x7 = x7(1:L_r);
% x8 = audioread('48_pass.wav'); x8 = x8(1:L_r);
% x9 = audioread('49_pass.wav'); x9 = x9(1:L_r);
% x10 = audioread('50_pass.wav'); x10 = x10(1:L_r);
% x11 = audioread('51_pass.wav'); x11 = x11(1:L_r);
% x12 = audioread('52_pass.wav'); x12 = x12(1:L_r);
% x13 = audioread('53_pass.wav'); x13 = x13(1:L_r);
% x14 = audioread('54_pass.wav'); x14 = x14(1:L_r);
% x15 = audioread('55_pass.wav'); x15 = x15(1:L_r);
% x16 = audioread('56_pass.wav'); x16 = x16(1:L_r);
% x17 = audioread('57_pass.wav'); x17 = x17(1:L_r);
% x18 = audioread('58_pass.wav'); x18 = x18(1:L_r);
% x19 = audioread('59_pass.wav'); x19 = x19(1:L_r);
% x20 = audioread('60_pass.wav'); x20 = x20(1:L_r);

% x1 = audioread('61_pass.wav'); x1 = x1(1:L_r);
% x2 = audioread('62_pass.wav'); x2 = x2(1:L_r);
% x3 = audioread('63_pass.wav'); x3 = x3(1:L_r);
% x4 = audioread('64_pass.wav'); x4 = x4(1:L_r);
% x5 = audioread('65_pass.wav'); x5 = x5(1:L_r);
% x6 = audioread('66_pass.wav'); x6 = x6(1:L_r);
% x7 = audioread('67_pass.wav'); x7 = x7(1:L_r);
% x8 = audioread('68_pass.wav'); x8 = x8(1:L_r);
% x9 = audioread('69_pass.wav'); x9 = x9(1:L_r);
% x10 = audioread('70_pass.wav'); x10 = x10(1:L_r);
% x11 = audioread('71_pass.wav'); x11 = x11(1:L_r);
% x12 = audioread('72_pass.wav'); x12 = x12(1:L_r);
% x13 = audioread('73_pass.wav'); x13 = x13(1:L_r);
% x14 = audioread('74_pass.wav'); x14 = x14(1:L_r);
% x15 = audioread('75_pass.wav'); x15 = x15(1:L_r);
% x16 = audioread('76_pass.wav'); x16 = x16(1:L_r);
% x17 = audioread('77_pass.wav'); x17 = x17(1:L_r);
% x18 = audioread('78_pass.wav'); x18 = x18(1:L_r);
% x19 = audioread('79_pass.wav'); x19 = x19(1:L_r);
% x20 = audioread('80_pass.wav'); x20 = x20(1:L_r);

% x1 = audioread('81_pass.wav'); x1 = x1(1:L_r);
% x2 = audioread('82_pass.wav'); x2 = x2(1:L_r);
% x3 = audioread('83_pass.wav'); x3 = x3(1:L_r);
% x4 = audioread('84_pass.wav'); x4 = x4(1:L_r);
% x5 = audioread('85_pass.wav'); x5 = x5(1:L_r);
% x6 = audioread('86_pass.wav'); x6 = x6(1:L_r);
% x7 = audioread('87_pass.wav'); x7 = x7(1:L_r);
% x8 = audioread('88_pass.wav'); x8 = x8(1:L_r);
% x9 = audioread('89_pass.wav'); x9 = x9(1:L_r);
% x10 = audioread('90_pass.wav'); x10 = x10(1:L_r);
% x11 = audioread('91_pass.wav'); x11 = x11(1:L_r);
% x12 = audioread('92_pass.wav'); x12 = x12(1:L_r);
% x13 = audioread('93_pass.wav'); x13 = x13(1:L_r);
% x14 = audioread('94_pass.wav'); x14 = x14(1:L_r);
% x15 = audioread('95_pass.wav'); x15 = x15(1:L_r);
% x16 = audioread('96_pass.wav'); x16 = x16(1:L_r);
% x17 = audioread('97_pass.wav'); x17 = x17(1:L_r);
% x18 = audioread('98_pass.wav'); x18 = x18(1:L_r);
% x19 = audioread('99_pass.wav'); x19 = x19(1:L_r);
% x20 = audioread('100_pass.wav'); x20 = x20(1:L_r);

% x1 = audioread('101_pass.wav'); x1 = x1(1:L_r);
% x2 = audioread('102_pass.wav'); x2 = x2(1:L_r);
% x3 = audioread('103_pass.wav'); x3 = x3(1:L_r);
% x4 = audioread('104_pass.wav'); x4 = x4(1:L_r);
% x5 = audioread('105_pass.wav'); x5 = x5(1:L_r);
% x6 = audioread('106_pass.wav'); x6 = x6(1:L_r);
% x7 = audioread('107_pass.wav'); x7 = x7(1:L_r);
% x8 = audioread('108_pass.wav'); x8 = x8(1:L_r);
% x9 = audioread('109_pass.wav'); x9 = x9(1:L_r);
% x10 = audioread('110_pass.wav'); x10 = x10(1:L_r);
% x11 = audioread('111_pass.wav'); x11 = x11(1:L_r);
% x12 = audioread('112_pass.wav'); x12 = x12(1:L_r);
% x13 = audioread('113_pass.wav'); x13 = x13(1:L_r);
% x14 = audioread('114_pass.wav'); x14 = x14(1:L_r);
% x15 = audioread('115_pass.wav'); x15 = x15(1:L_r);
% x16 = audioread('116_pass.wav'); x16 = x16(1:L_r);
% x17 = audioread('117_pass.wav'); x17 = x17(1:L_r);
% x18 = audioread('118_pass.wav'); x18 = x18(1:L_r);
% x19 = audioread('119_pass.wav'); x19 = x19(1:L_r);
% x20 = audioread('120_pass.wav'); x20 = x20(1:L_r);

% x1 = audioread('121_pass.wav'); x1 = x1(1:L_r);
% x2 = audioread('122_pass.wav'); x2 = x2(1:L_r);
% x3 = audioread('123_pass.wav'); x3 = x3(1:L_r);
% x4 = audioread('124_pass.wav'); x4 = x4(1:L_r);
% x5 = audioread('125_pass.wav'); x5 = x5(1:L_r);
% x6 = audioread('126_pass.wav'); x6 = x6(1:L_r);
% x7 = audioread('127_pass.wav'); x7 = x7(1:L_r);
% x8 = audioread('128_pass.wav'); x8 = x8(1:L_r);
% x9 = audioread('129_pass.wav'); x9 = x9(1:L_r);
% x10 = audioread('130_pass.wav'); x10 = x10(1:L_r);
% x11 = audioread('131_pass.wav'); x11 = x11(1:L_r);
% x12 = audioread('132_pass.wav'); x12 = x12(1:L_r);
% x13 = audioread('133_pass.wav'); x13 = x13(1:L_r);
% x14 = audioread('134_pass.wav'); x14 = x14(1:L_r);
% x15 = audioread('135_pass.wav'); x15 = x15(1:L_r);
% x16 = audioread('136_pass.wav'); x16 = x16(1:L_r);
% x17 = audioread('137_pass.wav'); x17 = x17(1:L_r);
% x18 = audioread('138_pass.wav'); x18 = x18(1:L_r);
% x19 = audioread('139_pass.wav'); x19 = x19(1:L_r);
% x20 = audioread('140_pass.wav'); x20 = x20(1:L_r);

% x1 = audioread('141_pass.wav'); x1 = x1(1:L_r);
% x2 = audioread('142_pass.wav'); x2 = x2(1:L_r);
% x3 = audioread('143_pass.wav'); x3 = x3(1:L_r);
% x4 = audioread('144_pass.wav'); x4 = x4(1:L_r);
% x5 = audioread('145_pass.wav'); x5 = x5(1:L_r);
% x6 = audioread('146_pass.wav'); x6 = x6(1:L_r);
% x7 = audioread('147_pass.wav'); x7 = x7(1:L_r);
% x8 = audioread('148_pass.wav'); x8 = x8(1:L_r);
% x9 = audioread('149_pass.wav'); x9 = x9(1:L_r);
% x10 = audioread('150_pass.wav'); x10 = x10(1:L_r);
% x11 = audioread('151_pass.wav'); x11 = x11(1:L_r);
% x12 = audioread('152_pass.wav'); x12 = x12(1:L_r);
% x13 = audioread('153_pass.wav'); x13 = x13(1:L_r);
% x14 = audioread('154_pass.wav'); x14 = x14(1:L_r);
% x15 = audioread('155_pass.wav'); x15 = x15(1:L_r);
% x16 = audioread('156_pass.wav'); x16 = x16(1:L_r);
% x17 = audioread('157_pass.wav'); x17 = x17(1:L_r);
% x18 = audioread('158_pass.wav'); x18 = x18(1:L_r);
% x19 = audioread('159_pass.wav'); x19 = x19(1:L_r);
% x20 = audioread('160_pass.wav'); x20 = x20(1:L_r);

% x1 = audioread('161_pass.wav'); x1 = x1(1:L_r);
% x2 = audioread('162_pass.wav'); x2 = x2(1:L_r);
% x3 = audioread('163_pass.wav'); x3 = x3(1:L_r);
% x4 = audioread('164_pass.wav'); x4 = x4(1:L_r);
% x5 = audioread('165_pass.wav'); x5 = x5(1:L_r);
% x6 = audioread('166_pass.wav'); x6 = x6(1:L_r);
% x7 = audioread('167_pass.wav'); x7 = x7(1:L_r);
% x8 = audioread('168_pass.wav'); x8 = x8(1:L_r);
% x9 = audioread('169_pass.wav'); x9 = x9(1:L_r);
% x10 = audioread('170_pass.wav'); x10 = x10(1:L_r);
% x11 = audioread('171_pass.wav'); x11 = x11(1:L_r);
% x12 = audioread('172_pass.wav'); x12 = x12(1:L_r);
% x13 = audioread('173_pass.wav'); x13 = x13(1:L_r);
% x14 = audioread('174_pass.wav'); x14 = x14(1:L_r);
% x15 = audioread('175_pass.wav'); x15 = x15(1:L_r);
% x16 = audioread('176_pass.wav'); x16 = x16(1:L_r);
% x17 = audioread('177_pass.wav'); x17 = x17(1:L_r);
% x18 = audioread('178_pass.wav'); x18 = x18(1:L_r);
% x19 = audioread('179_pass.wav'); x19 = x19(1:L_r);
% x20 = audioread('180_pass.wav'); x20 = x20(1:L_r);

% x1 = audioread('201_pass.wav'); x1 = x1(1:L_r);
% x2 = audioread('202_pass.wav'); x2 = x2(1:L_r);
% x3 = audioread('203_pass.wav'); x3 = x3(1:L_r);
% x4 = audioread('204_pass.wav'); x4 = x4(1:L_r);
% x5 = audioread('205_pass.wav'); x5 = x5(1:L_r);
% x6 = audioread('206_pass.wav'); x6 = x6(1:L_r);
% x7 = audioread('207_pass.wav'); x7 = x7(1:L_r);
% x8 = audioread('208_pass.wav'); x8 = x8(1:L_r);
% x9 = audioread('209_pass.wav'); x9 = x9(1:L_r);
% x10 = audioread('210_pass.wav'); x10 = x10(1:L_r);
% x11 = audioread('211_pass.wav'); x11 = x11(1:L_r);
% x12 = audioread('212_pass.wav'); x12 = x12(1:L_r);
% x13 = audioread('213_pass.wav'); x13 = x13(1:L_r);
% x14 = audioread('214_pass.wav'); x14 = x14(1:L_r);
% x15 = audioread('215_pass.wav'); x15 = x15(1:L_r);
% x16 = audioread('216_pass.wav'); x16 = x16(1:L_r);
% x17 = audioread('217_pass.wav'); x17 = x17(1:L_r);
% x18 = audioread('218_pass.wav'); x18 = x18(1:L_r);
% x19 = audioread('219_pass.wav'); x19 = x19(1:L_r);
% x20 = audioread('220_pass.wav'); x20 = x20(1:L_r);



% Fail
% x21 = audioread('1_fail.wav'); x21 = x21(1:L_r);
% x22 = audioread('2_fail.wav'); x22 = x22(1:L_r);
% x23 = audioread('3_fail.wav'); x23 = x23(1:L_r);
% x24 = audioread('4_fail.wav'); x24 = x24(1:L_r);
% x25 = audioread('5_fail.wav'); x25 = x25(1:L_r);
% x26 = audioread('6_fail.wav'); x26 = x26(1:L_r);
% x27 = audioread('7_fail.wav'); x27 = x27(1:L_r);
% x28 = audioread('8_fail.wav'); x28 = x28(1:L_r);
% x29 = audioread('9_fail.wav'); x29 = x29(1:L_r);
% x30 = audioread('11_fail.wav');  x30 = x30(1:L_r);
% x31 = audioread('12_fail.wav'); x31 = x31(1:L_r);
% x32 = audioread('13_fail.wav'); x32 = x32(1:L_r);
% x33 = audioread('14_fail.wav'); x33 = x33(1:L_r);
% x34 = audioread('15_fail.wav'); x34 = x34(1:L_r);
% x35 = audioread('16_fail.wav'); x35 = x35(1:L_r);
% x36 = audioread('17_fail.wav'); x36 = x36(1:L_r);
% x37 = audioread('18_fail.wav'); x37 = x37(1:L_r);
% x38 = audioread('19_fail.wav'); x38 = x38(1:L_r);
% x39 = audioread('20_fail.wav'); x39 = x39(1:L_r);
% x40 = audioread('21_fail.wav'); x40 = x40(1:L_r);

x21 = audioread('21_fail.wav'); x21 = x21(1:L_r);
x22 = audioread('22_fail.wav'); x22 = x22(1:L_r);
x23 = audioread('23_fail.wav'); x23 = x23(1:L_r);
x24 = audioread('24_fail.wav'); x24 = x24(1:L_r);
x25 = audioread('25_fail.wav'); x25 = x25(1:L_r);
x26 = audioread('26_fail.wav'); x26 = x26(1:L_r);
x27 = audioread('27_fail.wav'); x27 = x27(1:L_r);
x28 = audioread('28_fail.wav'); x28 = x28(1:L_r);
x29 = audioread('29_fail.wav'); x29 = x29(1:L_r);
x30 = audioread('30_fail.wav');  x30 = x30(1:L_r);
x31 = audioread('31_fail.wav'); x31 = x31(1:L_r);
x32 = audioread('32_fail.wav'); x32 = x32(1:L_r);
x33 = audioread('68_fail.wav'); x33 = x33(1:L_r);
x34 = audioread('72_fail.wav'); x34 = x34(1:L_r);
x35 = audioread('73_fail.wav'); x35 = x35(1:L_r);
x36 = audioread('74_fail.wav'); x36 = x36(1:L_r);
x37 = audioread('75_fail.wav'); x37 = x37(1:L_r);
x38 = audioread('76_fail.wav'); x38 = x38(1:L_r);
x39 = audioread('77_fail.wav'); x39 = x39(1:L_r);
x40 = audioread('78_fail.wav'); x40 = x40(1:L_r);


% x21 = audioread('79_fail.wav'); x21 = x21(1:L_r);
% x22 = audioread('80_fail.wav'); x22 = x22(1:L_r);
% x23 = audioread('81_fail.wav'); x23 = x23(1:L_r);
% x24 = audioread('82_fail.wav'); x24 = x24(1:L_r);
% x25 = audioread('83_fail.wav'); x25 = x25(1:L_r);
% x26 = audioread('84_fail.wav'); x26 = x26(1:L_r);
% x27 = audioread('85_fail.wav'); x27 = x27(1:L_r);
% x28 = audioread('86_fail.wav'); x28 = x28(1:L_r);
% x29 = audioread('87_fail.wav'); x29 = x29(1:L_r);
% x30 = audioread('88_fail.wav');  x30 = x30(1:L_r);
% x31 = audioread('89_fail.wav'); x31 = x31(1:L_r);
% x32 = audioread('90_fail.wav'); x32 = x32(1:L_r);
% x33 = audioread('91_fail.wav'); x33 = x33(1:L_r);
% x34 = audioread('92_fail.wav'); x34 = x34(1:L_r);
% x35 = audioread('93_fail.wav'); x35 = x35(1:L_r);
% x36 = audioread('101_fail.wav'); x36 = x36(1:L_r);
% x37 = audioread('96_fail.wav'); x37 = x37(1:L_r);
% x38 = audioread('97_fail.wav'); x38 = x38(1:L_r);
% x39 = audioread('98_fail.wav'); x39 = x39(1:L_r);
% x40 = audioread('99_fail.wav'); x40 = x40(1:L_r);

% x21 = audioread('79_fail.wav'); x21 = x21(1:L_r);
% x22 = audioread('80_fail.wav'); x22 = x22(1:L_r);
% x23 = audioread('101_fail.wav'); x23 = x23(1:L_r);
% x24 = audioread('102_fail.wav'); x24 = x24(1:L_r);
% x25 = audioread('103_fail.wav'); x25 = x25(1:L_r);
% x26 = audioread('104_fail.wav'); x26 = x26(1:L_r);
% x27 = audioread('105_fail.wav'); x27 = x27(1:L_r);
% x28 = audioread('106_fail.wav'); x28 = x28(1:L_r);
% x29 = audioread('107_fail.wav'); x29 = x29(1:L_r);
% x30 = audioread('108_fail.wav');  x30 = x30(1:L_r);
% x31 = audioread('109_fail.wav'); x31 = x31(1:L_r);
% x32 = audioread('110_fail.wav'); x32 = x32(1:L_r);
% x33 = audioread('111_fail.wav'); x33 = x33(1:L_r);
% x34 = audioread('112_fail.wav'); x34 = x34(1:L_r);
% x35 = audioread('113_fail.wav'); x35 = x35(1:L_r);
% x36 = audioread('114_fail.wav'); x36 = x36(1:L_r);
% x37 = audioread('116_fail.wav'); x37 = x37(1:L_r);
% x38 = audioread('117_fail.wav'); x38 = x38(1:L_r);
% x39 = audioread('118_fail.wav'); x39 = x39(1:L_r);
% x40 = audioread('119_fail.wav'); x40 = x40(1:L_r);

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
FrameN = 3; %length(x1)/L;


h = hann(L);
%h = 1;
x1_ps = zeros(L/2,1); 
x2_ps = zeros(L/2,1);
x3_ps = zeros(L/2,1); 
x4_ps = zeros(L/2,1);
x5_ps = zeros(L/2,1); 
x6_ps = zeros(L/2,1);
x7_ps = zeros(L/2,1); 
x8_ps = zeros(L/2,1);
x9_ps = zeros(L/2,1); 
x10_ps = zeros(L/2,1);
x11_ps = zeros(L/2,1); 
x12_ps = zeros(L/2,1);
x13_ps = zeros(L/2,1); 
x14_ps = zeros(L/2,1);
x15_ps = zeros(L/2,1); 
x16_ps = zeros(L/2,1);
x17_ps = zeros(L/2,1); 
x18_ps = zeros(L/2,1);
x19_ps = zeros(L/2,1); 
x20_ps = zeros(L/2,1);
x21_ps = zeros(L/2,1); 
x22_ps = zeros(L/2,1);
x23_ps = zeros(L/2,1); 
x24_ps = zeros(L/2,1);
x25_ps = zeros(L/2,1); 
x26_ps = zeros(L/2,1);
x27_ps = zeros(L/2,1); 
x28_ps = zeros(L/2,1);
x29_ps = zeros(L/2,1); 
x30_ps = zeros(L/2,1);
x31_ps = zeros(L/2,1); 
x32_ps = zeros(L/2,1);
x33_ps = zeros(L/2,1); 
x34_ps = zeros(L/2,1);
x35_ps = zeros(L/2,1); 
x36_ps = zeros(L/2,1);
x37_ps = zeros(L/2,1); 
x38_ps = zeros(L/2,1);
x39_ps = zeros(L/2,1); 
x40_ps = zeros(L/2,1);

for i = 1:FrameN
   x1h = x1( (i-1)*L+1 : L*i );
   x1h = x1h.*h;
   x1_ps_temp = abs(fft(x1h));
   x1_ps_temp = x1_ps_temp(1:L/2);
   x1_ps = x1_ps + x1_ps_temp;
   %
   x2h = x2( (i-1)*L+1 : L*i );
   x2h = x2h.*h;
   x2_ps_temp = abs(fft(x2h));
   x2_ps_temp = x2_ps_temp(1:L/2);
   x2_ps = x2_ps + x2_ps_temp;
   %
   x3h = x3( (i-1)*L+1 : L*i );
   x3h = x3h.*h;
   x3_ps_temp = abs(fft(x3h));
   x3_ps_temp = x3_ps_temp(1:L/2);
   x3_ps = x3_ps + x3_ps_temp;
   %
   x4h = x4( (i-1)*L+1 : L*i );
   x4h = x4h.*h;
   x4_ps_temp = abs(fft(x4h));
   x4_ps_temp = x4_ps_temp(1:L/2);
   x4_ps = x4_ps + x4_ps_temp;
   %
   x5h = x5( (i-1)*L+1 : L*i );
   x5h = x5h.*h;
   x5_ps_temp = abs(fft(x5h));
   x5_ps_temp = x5_ps_temp(1:L/2);
   x5_ps = x5_ps + x5_ps_temp;
   %
   x6h = x6( (i-1)*L+1 : L*i );
   x6h = x6h.*h;
   x6_ps_temp = abs(fft(x6h));
   x6_ps_temp = x6_ps_temp(1:L/2);
   x6_ps = x6_ps + x6_ps_temp;
   %
   x7h = x7( (i-1)*L+1 : L*i );
   x7h = x7h.*h;
   x7_ps_temp = abs(fft(x7h));
   x7_ps_temp = x7_ps_temp(1:L/2);
   x7_ps = x7_ps + x7_ps_temp;
   %
   x8h = x8( (i-1)*L+1 : L*i );
   x8h = x8h.*h;
   x8_ps_temp = abs(fft(x8h));
   x8_ps_temp = x8_ps_temp(1:L/2);
   x8_ps = x8_ps + x8_ps_temp;
   %
   x9h = x9( (i-1)*L+1 : L*i );
   x9h = x9h.*h;
   x9_ps_temp = abs(fft(x9h));
   x9_ps_temp = x9_ps_temp(1:L/2);
   x9_ps = x9_ps + x9_ps_temp;
   %
   x10h = x10( (i-1)*L+1 : L*i );
   x10h = x10h.*h;
   x10_ps_temp = abs(fft(x10h));
   x10_ps_temp = x10_ps_temp(1:L/2);
   x10_ps = x10_ps + x10_ps_temp;
   %
   x11h = x11( (i-1)*L+1 : L*i );
   x11h = x11h.*h;
   x11_ps_temp = abs(fft(x11h));
   x11_ps_temp = x11_ps_temp(1:L/2);
   x11_ps = x11_ps + x11_ps_temp;
   %
   x12h = x12( (i-1)*L+1 : L*i );
   x12h = x12h.*h;
   x12_ps_temp = abs(fft(x12h));
   x12_ps_temp = x12_ps_temp(1:L/2);
   x12_ps = x12_ps + x12_ps_temp;
   %
   x13h = x13( (i-1)*L+1 : L*i );
   x13h = x13h.*h;
   x13_ps_temp = abs(fft(x13h));
   x13_ps_temp = x13_ps_temp(1:L/2);
   x13_ps = x13_ps + x13_ps_temp;
   %
   x14h = x14( (i-1)*L+1 : L*i );
   x14h = x14h.*h;
   x14_ps_temp = abs(fft(x14h));
   x14_ps_temp = x14_ps_temp(1:L/2);
   x14_ps = x14_ps + x14_ps_temp;
   %
   x15h = x15( (i-1)*L+1 : L*i );
   x15h = x15h.*h;
   x15_ps_temp = abs(fft(x15h));
   x15_ps_temp = x15_ps_temp(1:L/2);
   x15_ps = x15_ps + x15_ps_temp;
   %
   x16h = x16( (i-1)*L+1 : L*i );
   x16h = x16h.*h;
   x16_ps_temp = abs(fft(x16h));
   x16_ps_temp = x16_ps_temp(1:L/2);
   x16_ps = x16_ps + x16_ps_temp;
   %
   x17h = x17( (i-1)*L+1 : L*i );
   x17h = x17h.*h;
   x17_ps_temp = abs(fft(x17h));
   x17_ps_temp = x17_ps_temp(1:L/2);
   x17_ps = x17_ps + x17_ps_temp;
   %
   x18h = x18( (i-1)*L+1 : L*i );
   x18h = x18h.*h;
   x18_ps_temp = abs(fft(x18h));
   x18_ps_temp = x18_ps_temp(1:L/2);
   x18_ps = x18_ps + x18_ps_temp;
   %
   x19h = x19( (i-1)*L+1 : L*i );
   x19h = x19h.*h;
   x19_ps_temp = abs(fft(x19h));
   x19_ps_temp = x19_ps_temp(1:L/2);
   x19_ps = x19_ps + x19_ps_temp;
   %
   x20h = x20( (i-1)*L+1 : L*i );
   x20h = x20h.*h;
   x20_ps_temp = abs(fft(x20h));
   x20_ps_temp = x20_ps_temp(1:L/2);
   x20_ps = x20_ps + x20_ps_temp;
   %
   x21h = x21( (i-1)*L+1 : L*i );
   x21h = x21h.*h;
   x21_ps_temp = abs(fft(x21h));
   x21_ps_temp = x21_ps_temp(1:L/2);
   x21_ps = x21_ps + x21_ps_temp;
   %
   x22h = x22( (i-1)*L+1 : L*i );
   x22h = x22h.*h;
   x22_ps_temp = abs(fft(x22h));
   x22_ps_temp = x22_ps_temp(1:L/2);
   x22_ps = x22_ps + x22_ps_temp;
   %
   x23h = x23( (i-1)*L+1 : L*i );
   x23h = x23h.*h;
   x23_ps_temp = abs(fft(x23h));
   x23_ps_temp = x23_ps_temp(1:L/2);
   x23_ps = x23_ps + x23_ps_temp;
   %
   x24h = x24( (i-1)*L+1 : L*i );
   x24h = x24h.*h;
   x24_ps_temp = abs(fft(x24h));
   x24_ps_temp = x24_ps_temp(1:L/2);
   x24_ps = x24_ps + x24_ps_temp;
   %
   x25h = x25( (i-1)*L+1 : L*i );
   x25h = x25h.*h;
   x25_ps_temp = abs(fft(x25h));
   x25_ps_temp = x25_ps_temp(1:L/2);
   x25_ps = x25_ps + x25_ps_temp;
   %
   x26h = x26( (i-1)*L+1 : L*i );
   x26h = x26h.*h;
   x26_ps_temp = abs(fft(x26h));
   x26_ps_temp = x26_ps_temp(1:L/2);
   x26_ps = x26_ps + x26_ps_temp;
   %
   x27h = x27( (i-1)*L+1 : L*i );
   x27h = x27h.*h;
   x27_ps_temp = abs(fft(x27h));
   x27_ps_temp = x27_ps_temp(1:L/2);
   x27_ps = x27_ps + x27_ps_temp;
   %
   x28h = x28( (i-1)*L+1 : L*i );
   x28h = x28h.*h;
   x28_ps_temp = abs(fft(x28h));
   x28_ps_temp = x28_ps_temp(1:L/2);
   x28_ps = x28_ps + x28_ps_temp;
   %
   x29h = x29( (i-1)*L+1 : L*i );
   x29h = x29h.*h;
   x29_ps_temp = abs(fft(x29h));
   x29_ps_temp = x29_ps_temp(1:L/2);
   x29_ps = x29_ps + x29_ps_temp;
   %
   x30h = x30( (i-1)*L+1 : L*i );
   x30h = x30h.*h;
   x30_ps_temp = abs(fft(x30h));
   x30_ps_temp = x30_ps_temp(1:L/2);
   x30_ps = x30_ps + x30_ps_temp;
   %
   x31h = x31( (i-1)*L+1 : L*i );
   x31h = x31h.*h;
   x31_ps_temp = abs(fft(x31h));
   x31_ps_temp = x31_ps_temp(1:L/2);
   x31_ps = x31_ps + x31_ps_temp;
   %
   x32h = x32( (i-1)*L+1 : L*i );
   x32h = x32h.*h;
   x32_ps_temp = abs(fft(x32h));
   x32_ps_temp = x32_ps_temp(1:L/2);
   x32_ps = x32_ps + x32_ps_temp;
   %
   x33h = x33( (i-1)*L+1 : L*i );
   x33h = x33h.*h;
   x33_ps_temp = abs(fft(x33h));
   x33_ps_temp = x33_ps_temp(1:L/2);
   x33_ps = x33_ps + x33_ps_temp;
   %
   x34h = x34( (i-1)*L+1 : L*i );
   x34h = x34h.*h;
   x34_ps_temp = abs(fft(x34h));
   x34_ps_temp = x34_ps_temp(1:L/2);
   x34_ps = x34_ps + x34_ps_temp;
   %
   x35h = x35( (i-1)*L+1 : L*i );
   x35h = x35h.*h;
   x35_ps_temp = abs(fft(x35h));
   x35_ps_temp = x35_ps_temp(1:L/2);
   x35_ps = x35_ps + x35_ps_temp;
   %
   x36h = x36( (i-1)*L+1 : L*i );
   x36h = x36h.*h;
   x36_ps_temp = abs(fft(x36h));
   x36_ps_temp = x36_ps_temp(1:L/2);
   x36_ps = x36_ps + x36_ps_temp;
   %
   x37h = x37( (i-1)*L+1 : L*i );
   x37h = x37h.*h;
   x37_ps_temp = abs(fft(x37h));
   x37_ps_temp = x37_ps_temp(1:L/2);
   x37_ps = x37_ps + x37_ps_temp;
   %
   x38h = x38( (i-1)*L+1 : L*i );
   x38h = x38h.*h;
   x38_ps_temp = abs(fft(x38h));
   x38_ps_temp = x38_ps_temp(1:L/2);
   x38_ps = x38_ps + x38_ps_temp;
   %
   x39h = x39( (i-1)*L+1 : L*i );
   x39h = x39h.*h;
   x39_ps_temp = abs(fft(x39h));
   x39_ps_temp = x39_ps_temp(1:L/2);
   x39_ps = x39_ps + x39_ps_temp;
   %
   x40h = x40( (i-1)*L+1 : L*i );
   x40h = x40h.*h;
   x40_ps_temp = abs(fft(x40h));
   x40_ps_temp = x40_ps_temp(1:L/2);
   x40_ps = x40_ps + x40_ps_temp;
end

%normalize and get rid off offset of power spectrum
x_ps = [x1_ps x2_ps x3_ps x4_ps x5_ps x6_ps x7_ps x8_ps x9_ps x10_ps x11_ps x12_ps x13_ps x14_ps x15_ps x16_ps x17_ps x18_ps x19_ps x20_ps x21_ps x22_ps x23_ps x24_ps x25_ps x26_ps x27_ps x28_ps x29_ps x30_ps x31_ps x32_ps x33_ps x34_ps x35_ps x36_ps x37_ps x38_ps x39_ps x40_ps];


% using pre-built pca transfer model's axis center
% Therefore new data need to minus pre-built offset, scaleMean.
x_ps = 20*log10(x_ps);
offset = repmat(scaleMean, 1, N);
x_ps = x_ps - offset;

%normalized
for i = 1:L/2
    x_ps(i,:) = x_ps(i,:)./sqrt( sum(x_ps(i,:).*x_ps(i,:)) );
end


%PCA transpose
pcsN = 10;
%transT = pcs'*x_ps;
transT = (V'*x_ps);
for i = 1:N/2
hold on
plot(transT(1:pcsN,i))
end
for i = N/2+1:N
plot(transT(1:pcsN,i), 'r')
end

%clustering
X = transT(:,1:2);
if(c1(1,1)>c1(2,1))
    fail = c1(1,:);
    pass = c1(2,:);
end

if(c1(1,1)<c1(2,1))
    fail = c1(2,:);
    pass = c1(1,:);
end
cluster = zeros(N,1);
for i =1:N
    transTr = transT(1:3,i)';
    distP = transTr - pass;
    distP = sum(distP.*distP);
    distF = transTr - fail;
    distF = sum(distF.*distF);
    if(distF>distP)
        cluster(i) = 1;
    else
        cluster(i) = 2;
    end
end
cluster'
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



