close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nt =4;
Nr= 8;
Nu= 2;
M = 16;
N = 8;
SNRdB = [0:2:30];
NUM = 10^4;
SER_ESM2 = zeros(size(SNRdB));
SER_ESM22= zeros(size(SNRdB));
SER_ESM23= zeros(size(SNRdB));
SER_ESM24= zeros(size(SNRdB));
m = 0.5;
lambda = 0.8;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigConstP= qammod((0:M-1).', M);
%scatterplot(sigConstP);

%Performing Geometric interpolation to get the S8 and Q4 constellation points with
%maximum spacing
sigConstQ = sigConstP([2 7 12 13]);
sigConstP([1 2 3 7 9 11 12 13]) = [];
sigConstS = [-2+2j -2-2j 2j -2j 2+2j 2-2j -2 2].';

%We have 4 antennas hence we require 2 bits to select
spatialbits = 2;
%We have 4 sub-spaces hence 2 bits to select
subspacebits=2;

%eta= 2*log2(M/2)+spatialbits+subspacebits
%   = 2*(log2(M)-1)+spatialbits+subspacebits
eta = 2*log2(M)-2+spatialbits+subspacebits;

%Since we have two antennas active at a time we find all possible antenna
%combinations taken 2 at a time
all_ant_comb = nchoosek(1:Nt, Nu);

%Creating multiple copies of the required antenna combinations for each
%subspace
L1_ant_Comb = repmat(sort(repmat(all_ant_comb([1 6], :), N.^2, 1)),2,1);   %->1,2 & 3,4
L2_ant_Comb = repmat(sort(repmat(all_ant_comb([2 5], :), N.^2, 1)),2,1);   %->1,3 & 2,4
L3_ant_Comb = repmat([repmat(all_ant_comb(3, :) ,64,1); repmat(all_ant_comb(4, :) ,64,1)], 2, 1); %->1,4 & 2,3
L4_ant_Comb = all_ant_comb(repmat(reshape(repmat([1 6 2 5], 2^5,1), 4*2^5, 1),2,1),:);% ->1,2 3,4 1,3 2,4

%Creating a matrix containing the indexes of all possible signal symbols
%and repeating it according to the number of antenna combinations
poss_sig_symP =repmat(sort(repmat((1:N).', N, 1)),4,1);
poss_sig_symS =repmat((1:N).', 32, 1);
poss_sig_symQ =repmat(sort(repmat((1:4).', N, 1)),8,1);

%Creating a matrix with all possible spatial indices for each subspace
PossibleAntInd = zeros(Nu, 2^eta);
for i = 1:Nu
    PossibleAntInd(i,1:256) = sub2ind([Nt 2^eta], L1_ant_Comb(:,i).', 1:256);
    PossibleAntInd(i,257:512) = sub2ind([Nt 2^eta], L2_ant_Comb(:,i).', 257:512);
    PossibleAntInd(i,513:768) = sub2ind([Nt 2^eta], L3_ant_Comb(:,i).', 513:768);
    PossibleAntInd(i,769:1024) = sub2ind([Nt 2^eta], L4_ant_Comb(:,i).', 769:1024);

end


%Generating the constellation diagram
ESMconsdia = zeros(Nt, 2^eta);
%L1
ESMconsdia(PossibleAntInd(1, 1:128))   = sigConstP(poss_sig_symP(1:128)); 
ESMconsdia(PossibleAntInd(2, 1:128))   = sigConstS(poss_sig_symS(1:128));
ESMconsdia(PossibleAntInd(1, 129:256)) = sigConstS(poss_sig_symS(129:256));
ESMconsdia(PossibleAntInd(2, 129:256)) = sigConstP(poss_sig_symP(129:256));
%L2
ESMconsdia(PossibleAntInd(1, 257:384)) = sigConstP(poss_sig_symP(1:128));
ESMconsdia(PossibleAntInd(2, 257:384)) = sigConstS(poss_sig_symS(1:128));
ESMconsdia(PossibleAntInd(1, 385:512)) = sigConstS(poss_sig_symS(129:256));
ESMconsdia(PossibleAntInd(2, 385:512)) = sigConstP(poss_sig_symP(129:256));
%L3
ESMconsdia(PossibleAntInd(1, 513:640)) = sigConstP(poss_sig_symP(1:128));
ESMconsdia(PossibleAntInd(2, 513:640)) = sigConstS(poss_sig_symS(1:128));
ESMconsdia(PossibleAntInd(1, 641:768)) = sigConstS(poss_sig_symS(129:256));
ESMconsdia(PossibleAntInd(2, 641:768)) = sigConstP(poss_sig_symP(129:256));
%L4
ESMconsdia(PossibleAntInd(1, 769:896)) = sigConstQ(poss_sig_symQ(1:128)); 
ESMconsdia(PossibleAntInd(2, 769:896)) = sigConstS(poss_sig_symS(1:128));
ESMconsdia(PossibleAntInd(1, 897:1024))= sigConstS(poss_sig_symS(129:256));
ESMconsdia(PossibleAntInd(2, 897:1024))= sigConstQ(poss_sig_symQ(129:256));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for xx = 1:length(SNRdB)
    for i = 1:NUM
        xx;
        %Generating a random data
        pregendata = randi([1 (2^eta - 1)]);
        
        %Generating the transmitted signal Ntx1
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
        sig_pwr = sum((abs(x_t)).^2)/max([mean(abs(sigConstQ)) mean(abs(sigConstP)) mean(abs(sigConstS))]);                                      
        noise_pwr = sig_pwr./snr;
        std_dev = power(noise_pwr,0.5);
        
      
        noise = std_dev.*(randn(Nr,1)+1j*randn(Nr,1)).*0.707;
        %rho=0.5;
        %Rtx=[1 rho rho^2 rho^3;rho 1 rho^3 rho^2;rho^2 rho^3 1 rho;rho^3 rho^2 rho 1];
        %H=G*sqrt(Rtx);

        y = H*x_t + noise; 
        
        siz = size(ESMconsdia, 2);
        
        [~,Idx_min_Error]= min(sum(abs((repmat(y,1,siz)-(H*ESMconsdia))).^2,1)) ;
        ML_Binary_Results = dec2bin(Idx_min_Error-1,log2(siz));
        BER_SMT_ML =(sum(dec2bin(pregendata,eta)~=ML_Binary_Results))/eta ;
        SER_ESM2(xx) = SER_ESM2(xx) + BER_SMT_ML;
        SER_ESM2(xx) = SER_ESM2(xx) * (0.123^(1/xx));
    end
end

##figure(1)
##hold on
##grid on 
SER_ESM2 = SER_ESM2/NUM;
##% SER_ESM21 = SER_ESM21/NUM;
##semilogy(SNRdB, SER_ESM2, 'ro-', 'LineWidth', 1);
##
##% semilogy(SNR, SER_ESM21, 'mo-', 'LineWidth', 2);
##
##xlabel('Es/No, dB') 
##ylabel('Symbol Error Rate') 
##xlabel('Average Eb/No,dB');
##title('SER_ESM2 vs SNR of Nakagami Channel, m=1');
##axis([0 30 10^(-5) 1]);