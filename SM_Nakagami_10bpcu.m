hold on;
clc;


SNRdB = [0:2:50];
SER_SM= zeros(size(SNRdB));

M= 256;  %size of signal constellation diagram

NUM=100;
Nt=8;          %no of transmit antennas 
Nr=4;          %no or recieving antennas


m = 1;
%lambda = 0.8;

%%% CONSTRUCTING AN M- QAM MODULATOR AND GENERATING M- QAM SIGNAL CONSTELLATION DIAGRAM %%%
S=qammod([0:M-1], M)'; 

%%% CALCULATING AVG. POWER OF SIGNAL CONSTELLATION DIAGRAM %%%
S_avg=S'*S/M;

%%% NORMALISING THE POWER OF THE SIGNAL CONSTELLATION DIAGRAM %%%
S_diag=S/sqrt(S_avg);

%%% NO. OF BITS TRANSMITTED IN SPATIAL DOMAIN %%%
SpatialBits = floor ( log2 ( nchoosek (Nt , 1 ) ) ) ; 
%nchoosek :number of combinations of Nt antennas taken nu at a time

%%% CALCULATING SPECTRAL EFFICIENCY %%%
no=log2(M)+ SpatialBits;

%%% GENERATING ALL POSSIBLE ACTIVE ANTENNA COMBINATIONS TAKING nu ANTENNAS AT A TIME %%%
All_Ant_Comb = flipud ( nchoosek ( 1 :Nt , 1 ) ) ;
%C = nchoosek(v,k) returns a matrix containing all possible combinations of the elements of vector v taken k at a time. 
%Matrix C has k columns and m!/((m–k)! k!) rows, where m is length(v).

%%% Spectral efficiency is 4 bits -> 16 number of symbol and antenna combinations %%%
%%% Considering 4 symbols each combination, we consider 4 combinations of antenna %%%

%%% TAKING A SUBSET OF ALL ANTENNAS WHICH IS POSSIBLE ANTENNAS %%%
Possible_Ant_Comb = All_Ant_Comb ( reshape ( repmat ( 1:2^ SpatialBits ,M, 1 ) , 2 ^ no , 1 ) , : ) ;

%%% GENERATING GSM CONSTELLATION DIAGRAM %%%
% Creating a matrix containing the indexes of all possible signal symbols (1,2,3,4)
% repeating 4 times for each antenna combination
PossibleSigSymb = repmat ( ( 1 :M) .' , 2 ^ SpatialBits , 1 ) ;

% Creating a matrix containing all possible spatial indexes
PossibleAnt_Ind = zeros (1 ,2^ no ) ;


  PossibleAnt_Ind  = sub2ind ( [Nt 2^ no ] , Possible_Ant_Comb .' , 1:2^ no ) ;


%%% GSM CONSTELLATION DIAGRAM %%%
GSM_Cons_Diagram = zeros (Nt ,2^ no) ;
GSM_Cons_Diagram( PossibleAnt_Ind ) = repmat ( S_diag( PossibleSigSymb ) .' , 1 , 1 ) / sqrt ( 1 );

R_tx = eye(Nt(1));  %Identity Matrix 4x4
 R_rx = eye(Nr);    %4x4
 R_s = kron(R_rx,R_tx);     %16x16 : Multiplication of every possible pair
 LH = R_s;
eta=log2(size(GSM_Cons_Diagram,2));   %log2M
for xx=1:length(SNRdB)      % 1: 16
    for i=1:NUM           % 1: 10^5
    xx;
    pregendata= randi([1 (2^eta-1)]);       %pregenerated data
    %%Uniformly distributed pseudorandom integers from any one of the M constellation symbols
    
    x_t = GSM_Cons_Diagram(:,pregendata+1);    %4 bit transmitted data corresponding to a symbol
    %% Nakagami H_cap (NR x Nt)
        H_cap = sqrt(sum(abs(randn(Nr,Nt,m)/sqrt(2*m)).^2,3))+1j.*sqrt(sum(abs(randn(Nr,Nt,m)/sqrt(2*m)).^ 2 ,3 )) ;
        
##        %% correlation matrices sigma_r and sigma_t using exponential decay model
##        for i = 1:Nr
##          for j = 1:Nr
##              sigma_r(i, j) = lambda^(abs(i-j));
##          end
##        end
##        
##        for i = 1:Nt
##          for j = 1:Nt
##              sigma_t(i, j) = lambda^(abs(i-j));
##          end
##        end
        
        %% Final H
##        H = sigma_r^(0.5)*H_cap*sigma_t^(0.5);

H = H_cap;
        
    sigma = 1/10^(SNRdB(xx)/10);  %%variance
    sigmae=sigma;
    noise = sqrt(sigma)*(randn(Nr,1)+1j*randn(Nr,1))/sqrt(2);  %%AWGN
    %y = H_SC*x_t + noise ;
    y = H*x_t + noise ;%%without sc
    e= sqrt(sigmae)*(randn(Nr,1)+1j*randn(Nr,1))/sqrt(2);
    %H_SC_CSE=H_SC-0;
    H_SC_CSE=H-e;%%without sc
    siz= size(GSM_Cons_Diagram,2);
    % Ca l cul a t ing the Eucl idean di s t anc e s of a l l pos s ibl e SMT transmitted
    % vector s , and then f inding the index of the SMT symbol with the
    % minimum Euc l ide an di s t anc e , i . e . e r r o r .
    [~,Idx_min_Error]= min(sum(abs((repmat(y,1,siz)-(H_SC_CSE*GSM_Cons_Diagram))).^2,1)) ;
    % Converting the index to binary , to r e t r i e v e the transmitted binary b i t s
    ML_Binary_Results = dec2bin(Idx_min_Error-1,log2(siz));
    BER_SMT_ML =(sum(dec2bin(pregendata,eta)~=ML_Binary_Results))/eta; 
  
       SER_SM(xx)=SER_SM(xx)+BER_SMT_ML;
    
    end
end

figure(2)
grid on; 
hold on; 
SER_SM= SER_SM/(NUM); 


semilogy(SNRdB,SER_SM,'go-','LineWidth',2);

ylabel('Symbol Error Rate') 
xlabel('Average Eb/No,dB');
##%title('SER vs SNR of GSM Nakagami Channel, m=1');
axis([0 50 10^(-5) 1]);





