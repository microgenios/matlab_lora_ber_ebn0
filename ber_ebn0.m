clear
clc
txtD = [];
writematrix(txtD,'C:\Users\microgenios\Documents\dados.txt')

BW = 125e3 ;
fc = 915e6 ;
Power = 20 ;
SF = 12:1:12; % Spreading Factor from 7 to 12
Fs = 10e6 ;
%Fs = BW;
Fc = 921.5e6 ;%921.500 000
CR = 1;
SNR_dB = -25:1:25;           % SNR in DB
Total_iterations = 2000;

%snrr = zeros(1,Total_iterations);
%bersnr = zeros(Total_iterations, Total_iterations);



%Avg_BER = zeros(length(1:1:length(SF)), length(bit_error));

a = 17; b = 17; 
r = (b-a).*rand(1000,1) + a;
info_word_length = floor(r(1));

bit_error = zeros(1,Total_iterations);
ber = zeros(length(SF),length(SNR_dB)*length(SF)*Total_iterations);
ebn0 = zeros(length(SF),length(SNR_dB)*length(SF)*Total_iterations);
snr = zeros(length(SF),length(SNR_dB)*length(SF)*Total_iterations);
count = 0;

for sfreq = 1:1:length(SF)
    count = 0;
    bit_error = zeros(1,Total_iterations);
    ber = zeros(length(SF),length(SNR_dB)*length(SF)*Total_iterations);
    ebn0 = zeros(length(SF),length(SNR_dB)*length(SF)*Total_iterations);
    snr = zeros(length(SF),length(SNR_dB)*length(SF)*Total_iterations);

    for snr_r = 1:1:length(SNR_dB)
        for ite = 1:1:Total_iterations   
        
			%BW = 125e3 ;
			%fc = 915e6 ;
			%Power = 14 ;

			%Fs = 10e6 ;
			%Fc = 921.5e6 ;
            
			message_in = randstr([1 info_word_length],2,'useWildCards',false);
            
            %message_in = 'pxiM2BfpitNossU8iyHAS3ncpfMoVlUjFX6rys3h';
            %info_word_length = length(message_in);
			
            [signalIQ, signal,pld_swp, packet_hdr] = LoRa_Tx(message_in,BW,SF(sfreq),Power,Fs,Fc - fc);
			[message_out] = LoRa_Rx(signalIQ,BW,SF(sfreq),2,Fs,Fc - fc,SNR_dB(snr_r)) ;
%message_out
			%fprintf('SF = %d\n', SF(sfreq));
			%fprintf('SNR = %d\n', SNR_dB(snr_r));
			%fprintf('Message In = %s\n', message_in);
			%fprintf('Message Out = %s\n', char(message_out) );

			%mess_in= uint8(convertStringsToChars(message_in));
			%mess_out = message_out;

			bit_count = 0;
			if message_out(isnan(message_out)==0)
				
				m_in = cell2mat( cellstr( reshape( string(dec2bin(message_in,8)),size(message_in,1),[]) ) ) - '0';
				m_out = cell2mat( cellstr( reshape( string(dec2bin(message_out,8)),size(message_out,1),[]) ) ) - '0';

				if length(m_in) == length(m_out)
					bit_count=bit_count+length(find(m_out~=m_in));
				elseif length(m_in) < length(m_out)
					bit_count=bit_count+length(find(m_out(1:length(m_in))~=m_in));
					bit_count=bit_count+(length(m_out) - length(m_in));
				elseif length(m_in) > length(m_out)
					bit_count=bit_count+length(find(m_out~=m_in(1:length(m_out))));
					bit_count=bit_count+(length(m_in) - length(m_out));
				end
			else 
				bit_count=bit_count+(info_word_length*8); 
			end
            bit_error(ite) = bit_count;
            %bit_count
            %info_word_length
            %info_word_length
            %bit_error(ite)
            clear message_in;
            clear signalIQ;
            clear signal; 
            clear pld_swp;
            clear packet_hdr;
            clear message_out;
            clear m_in;
            clear m_out;
        end

        
        RB = SF(sfreq) * ((4.0/(4.0+CR))/( pow2(SF(sfreq))/BW ));
        proc_gain = 10 * log10(BW/RB);
        eb_n0 = SNR_dB(snr_r) + proc_gain;

        count = count + 1

        total = sum(bit_error)/((info_word_length*8)*Total_iterations);
        ebn0(sfreq,count) = eb_n0;
        ber(sfreq,count) = total;

        snr(sfreq,count) = SNR_dB(snr_r);

        txtD = [SF(sfreq) SNR_dB(snr_r) eb_n0 total];
        writematrix(txtD,'C:\Users\microgenios\Documents\dados.txt','WriteMode','append');
        %save('C:\Users\microgenios\Documents\pqfile.txt','eb_n0','total','-ascii')
       


         %bit_error
         %SF(sfreq)
         %snr_r
         %mean(bit_error)
         %eb_n0

        %Rb = Fs*BW/(2^SF(sfreq));
        %eb_n0 = SNR_dB(snr_r) + 10 * log10(BW/Rb);
        %EBN0(ite, snr) = eb_n0;

        %berr2(1,snr_r) = mean(berr);
        %%snrr(ite) = SNR_dB(snr_r);
    end



%berr2
	%Avg_BER(sfreq,:) = berr2; 
	%resut = mean(berr);
	%semilogy(SNR_dB,Avg_BER);
    figure(1)
    semilogy(snr,ber);
	hold on

    figure(2)
    semilogy(ebn0,ber);
	hold on  


end

%% Plotting 
% Plotting the BER vs SNR curve
% SNR_dB
% Avg_BER

%semilogy(SNR_dB,Avg_BER);
figure(1)
title('BER vs SNR');
xlabel('SNR');
ylabel('BER (Bit Error Rate)');
%legend('SF7','SF8','SF9','SF10','SF11','SF12');
grid on;
hold off

figure(2)
title('BER vs EbN0');
xlabel('EbN0 (Energy per bit to Noise Power Spectral Density Ratio)');
ylabel('BER (Bit Error Rate)');
%legend('SF7','SF8','SF9','SF10','SF11','SF12');
grid on;
hold off
