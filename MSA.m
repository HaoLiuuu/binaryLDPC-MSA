clear 
clc


H = [1 1 1 1 0 0 0 0 0 0;...
     1 0 0 0 1 1 1 0 0 0;... 
     1 1 0 0 1 0 0 1 1 0;... 
     0 0 1 1 0 1 0 1 0 1;... 
     0 0 0 1 0 0 1 0 1 1];



%%----------------------------- 相关参数
code_length = size(H,2);
[hRowIndex,...
 hColIndex,...
 prefix_index,...
 find_out_llr_index,...
 zero_index,...
 mink_tail,...
 m,n] = Get_All_Index(H);

R = 0.6;

threads = 24;

codeword_bits       =       zeros(1,code_length);
transmit_signal     =       codeword_bits * (-2) + 1;
receive_signal      =       zeros(threads,code_length);
out_llr             =       zeros(threads,code_length);
hard_decision       =       zeros(threads,code_length);
syndrome            =       zeros(threads,1);


alpha               =       zeros(size(hRowIndex,1),size(hRowIndex,2));
in_llr              =       zeros(size(hRowIndex,1),size(hRowIndex,2));
alpha_min           =       zeros(size(hRowIndex,1),size(hRowIndex,2));
gamma               =       zeros(size(hRowIndex,1),size(hRowIndex,2));
gamma_min           =       zeros(size(hRowIndex,1),size(hRowIndex,2));
alpha_gamma         =       zeros(size(hRowIndex,1),size(hRowIndex,2));
gamma_mink          =       zeros(size(hRowIndex,1),size(hRowIndex,2));

%%----------------------------- 信噪比参数
EBN0 = 1.0:1:8;
BER_LOG = zeros(1,length(EBN0));
BLER_LOG = zeros(1,length(EBN0));
%%----------------------------- 仿真
for k = 1 : length(EBN0)
    
    BER = 0;
    BLER = 0;
    cw_num = 0;
    SNR = EBN0(k) - 10 * log10(1 / (2 * R));
    tic
    %----------------------------------------------------------------------
    while BER < 200 && BLER < 4
        for l = 1 : threads
            
            receive_signal(l,:)     =       awgn(transmit_signal,SNR);
            % 加inf是为了保证hRowIndex有0的位置取到inf(非规则的情况)
            temp_receive_signal     =       [inf,receive_signal(l,:)];
            in_llr                  =       temp_receive_signal(hRowIndex);
            for i = 1 : 10 
                alpha               =       sign(in_llr);
                alpha_min           =       prod(alpha,2) .* alpha;
                %----------------------------------------------------------------------
                gamma                       =       abs(in_llr);
                gamma_mink                  =       [mink(gamma,2,2),mink_tail];
                index                       =       prefix_index + ((gamma == gamma_mink(:,1)) + 1);
                gamma_mink                  =       gamma_mink';
                gamma_min                   =       gamma_mink(index);
                alpha_gamma                 =       0.725 * (alpha_min .* gamma_min);
                alpha_gamma(zero_index)     =       0;
                temp_alpha_gamma            =       [zeros(m,1),alpha_gamma];
                %----------------------------------------------------------------------
                out_llr(l,:)        =       receive_signal(l,:) + sum(temp_alpha_gamma(find_out_llr_index));
                temp_out_llr        =       [inf,out_llr(l,:)];
                in_llr              =       temp_out_llr(hRowIndex) - alpha_gamma;           
                %----------------------------------------------------------------------
                hard_decision(l,:)      =   ~(out_llr(l,:)  > 0);
                syndrome(l,1)           =   sum(mod(sum(hard_decision(l,:) .* H,2),2));
                if syndrome(l,1)  == 0
                    break;
                end
            end
        end
        cw_num = cw_num + threads;
        total_error = sum(sum(hard_decision ~= codeword_bits));
        BER = BER + total_error;
        BLER = BLER + sum(sum(hard_decision,2) > 0);
    end
    BER_LOG(k) = BER / (cw_num * code_length);
    BLER_LOG(k) = BLER / cw_num;
   disp(['仿真码字数量:',num2str(cw_num),'   ',...
          'Eb/N0: ',num2str(EBN0(k)),' dB','   ',...
          'BER:  ',num2str(BER_LOG(k)),'   ', ...
          'BLER: ', num2str( BLER_LOG(k))]);
    toc
end


semilogy(EBN0,BER_LOG,EBN0,BLER_LOG);
grid on;
xlabel('EbN0 (dB)');
ylabel('BER / BLER');
























