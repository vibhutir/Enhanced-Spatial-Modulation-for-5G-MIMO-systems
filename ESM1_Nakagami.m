close all;
clear all;
clc;
Nt =4;
Nr= 8;
Nu= 2;
M = 16;
N = 8;
SNRdB = [0:2:30];
NUM = 10^4;
SER_ESM1 = zeros(size(SNRdB));
m = 1;
lambda = 0.8;
% SER_ESM11 = zeros(size(SNR));


sigConstP = qammod((0:M-1).', M);
% scatterplot(sigConstP);
%sigConstS = [4+6j 4-6j -4+6j -4-6j].';% 2+4j 2-4j -2+4j -2-4j 6+4j 6-4j -6+4j -6-4j 4+2j 4-2j -4+2j -4-2j].';    %assuming m = 16
%sigConstS=[2+2i 2-2i -2+2i -2-2i];% these points are exactly interpolated for 16. 
sigConstS = [2 -2 2j -2j 2+2j 2-2j -2-2j -2+2j];%for n = 8

spatialbits = 3;
%2^(spatialbits+log2(N))
eta = log2(M) + log2(N) + spatialbits;

all_ant_comb = nchoosek(1:Nt, Nu);
all_ant_comb([3 4], :) = [];
all_ant_comb = repmat(all_ant_comb, 2,1);

poss_ant_comb = all_ant_comb(reshape(repmat(1:2^spatialbits, M*N, 1),2^eta,1) ,:);

poss_sig_symP = sort(repmat((1:M)', N, 1));
%poss_sig_symP = sort(poss_sig_symP);
poss_sig_symS = repmat((1:N)', 2^(spatialbits), 1);

PossibleAntInd = zeros(Nu, 2^eta);

for i = 1:Nu
    PossibleAntInd(i, :) = sub2ind([Nt 2^eta], poss_ant_comb(:,i).', 1:2^eta);
end

ESMconsdia = zeros(Nt, 2^eta);
ESMconsdia(PossibleAntInd(1, 1:2^(eta-1))) = repmat(sigConstP(poss_sig_symP), 1,2^(eta-1)/size(poss_sig_symP,1));
ESMconsdia(PossibleAntInd(1, 2^(eta-1)+1:2^eta)) = repmat(sigConstS(poss_sig_symS), 1, 2^(eta-1)/size(poss_sig_symS,1));
ESMconsdia(PossibleAntInd(2, 1:2^(eta-1))) = repmat(sigConstS(poss_sig_symS), 1, 2^(eta-1)/size(poss_sig_symS,1));
ESMconsdia(PossibleAntInd(2, 2^(eta-1)+1:2^eta)) = repmat(sigConstP(poss_sig_symP), 1, 2^(eta-1)/size(poss_sig_symP,1));

for xx = 1:length(SNRdB)
    for i = 1:NUM
        xx;
        pregendata = randi([0 (2^eta - 1)]);
        
        x_t = ESMconsdia(:,pregendata+1);
        
        %% Nakagami H_cap (NR x Nt)
        H_cap = sqrt(sum(abs(randn(Nr,Nt,m)/sqrt(2*m)).^2,3))+1j.*sqrt(sum(abs(randn(Nr,Nt,m)/sqrt(2*m)).^ 2 ,3 )) ;
        
        %% correlation matrices sigma_r and sigma_t using exponential decay model
        for i = 1:Nr
          for j = 1:Nr
              sigma_r(i, j) = lambda^(abs(i-j));
          end
        end
        
        for i = 1:Nt
          for j = 1:Nt
              sigma_t(i, j) = lambda^(abs(i-j));
          end
        end
        
        %% Final H
        H = sigma_r^(0.5)*H_cap*sigma_t^(0.5);
        
        snr = power(10,(SNRdB(xx)/10));
        sig_pwr = sum((abs(x_t)).^2)/max([mean(abs(sigConstP)) mean(abs(sigConstS))]);                                   % this is very important...
        noise_pwr = sig_pwr./snr;
        std_dev = power(noise_pwr,0.5);
        
        %%
        %Calculation of sigma in line 52 is wrong. firstly calculate the
        %abs power then divide the signal power by the max of abs
        %esmconsdia, then equal power will be transmitted otherwise unequal
        %power leads to better BER that is not correct, means it is not
        %averaged out. so calculation what i have shown here should be
        %followed to everything from now on.
%         noise = sqrt(sigma)*(randn(Nr,1)+1j*randn(Nr,1))/sqrt(2);
        noise = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        
        y = H*x_t + noise;  %there is no spatial correlation
        
%         e = sqrt(sigma)*(randn(Nr,1)+1j*randn(Nr,1))/sqrt(2);
%         H_CSE = H-e;
%         y_CSE = H_CSE*x_t + noise;
        
        siz = size(ESMconsdia, 2);
        
        [~,Idx_min_Error]= min(sum(abs((repmat(y,1,siz)-(H*ESMconsdia))).^2,1)) ;
        ML_Binary_Results = dec2bin(Idx_min_Error-1,log2(siz));
        BER_SMT_ML =(sum(dec2bin(pregendata,eta)~=ML_Binary_Results))/eta ;
        
%         [~,Idx_min_Error1]= min(sum(abs((repmat(y_CSE,1,siz)-(H_CSE*ESMconsdia))).^2,1)) ;
%         ML_Binary_Results1 = dec2bin(Idx_min_Error1-1,log2(siz));
%         BER_SMT_ML1 =(sum(dec2bin(pregendata,eta)~=ML_Binary_Results1))/eta ;
        
        SER_ESM1(xx) = SER_ESM1(xx) + BER_SMT_ML;
%         SER_ESM11(xx) = SER_ESM11(xx) + BER_SMT_ML1;
    end
end

##figure(1)
##grid on 
##hold on 
SER_ESM1 = SER_ESM1/NUM;
##% SER_ESM11 = SER_ESM11/NUM;
##semilogy(SNRdB, SER_ESM1, 'mo-', 'LineWidth', 1);
##
##% semilogy(SNR, SER_ESM11, 'mo-', 'LineWidth', 2);
##
##ylabel('Symbol Error Rate') 
##xlabel('Average Eb/No,dB');
##%title('SER_ESM1 vs SNR of ESM-1 (Nakagami Channel, m=1)');
##axis([0 30 10^(-7) 1]);
        
        
        
        
        
        

