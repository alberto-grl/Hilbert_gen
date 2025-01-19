% Hilbert filter generator and visuals for Verilog
% Copyright (C) 2025 Alberto Garlassi I4NZX
% Part of the code comes from
% https://dsp.stackexchange.com/questions/50015/discrepancies-between-fft-based-hilbert-transform-and-fir-filter-results

pkg load signal;
n_samples = 1024
f_sample = 62500/4;
f_sig =1500;
f_lo = 1000;
N = 39;      %7 filter length of FIR Hilbert transformer is 2*N+1
bits_hilbert = 12; % n bits of integer representation of Hilbert coefficients
SIZE = 12; % n bits of integer representation of signal
idx = (0:n_samples-1)';
sig = sin((2*pi *  f_sig / f_sample) .* idx);
lo_sin = sin((2*pi *  f_lo / f_sample) .* idx);
lo_cos = cos((2*pi *  f_lo / f_sample) .* idx);
x0_sin= sig .* lo_sin;
x0_cos = sig .* lo_cos;
Y = fft(x0_sin);

P2 = abs(Y/n_samples);
P1 = P2(1:n_samples/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = f_sample*(0:(n_samples/2))/n_samples;

figure(1)
plot(f,P1)
title("Input signals");

h = hilbert(x0_sin);
nn = (1:N)';
g = 2*(sin(pi*nn/2)).^2./ (pi*nn);
g = [-g(N:-1:1);0;g];
hw = hamming(2*N+1);
g = g.*hw(:);    % windowed coefficients
I = x0_cos;
Q = filter(g,1,x0_cos);
h2 = I(1:end-N) + 1i*Q(N+1:end);
A1 = abs(h);
A2 = abs(h2);

I_sin = x0_sin;
idx2 = idx + round ((f_sample /f_sig )/4);
SSB =  I_sin(1:end-N)+ Q(N+1:end);


Y_SSB = fft(SSB);

P2 = abs(Y_SSB/(n_samples-N));
P1 = P2(1:(n_samples-N)/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = f_sample*(0:((n_samples-N)/2))/(n_samples-N);

figure(3);
plot(f,P1);
title('FFT SSB');


g_int = round(g * 2^bits_hilbert);    % windowed coefficients
I = x0_sin;
Q_int = filter(g_int,1,x0_sin)/2^bits_hilbert;

figure (4);
idx2 = idx + round ((f_sample /f_sig )/4);
plot (idx2(1:end-N), I(1:end-N), idx(1:end-N), Q(N+1:end ));
title('I and Q, real');

sig_int = round(sin((2*pi *  f_sig / f_sample) .* idx) * 2^SIZE);
lo_int = round(sin((2*pi *  f_lo / f_sample) .* idx) * 2^SIZE);
x0_sin_int = sig_int .* lo_int;
I = x0_sin_int / 2^SIZE;
Q_int = filter(g_int,1,x0_sin_int)/ 2^bits_hilbert / 2^SIZE ;

figure (5);
idx2 = idx + round ((f_sample /f_sig )/4);
plot (idx2(1:end-N), I(1:end-N), idx(1:end-N), Q_int(N+1:end ));
title('I and Q, integer');

f_id = fopen('C:\Users\alberto\Lattice\hilbert.v', 'w');
fprintf(f_id, "/*\n   Hilbert filter for Verilog\n");
fprintf(f_id,     "   Copyright (C) 2025 Alberto Garlassi I4NZX\n \
  GPL v3\n*/\n\n");

fprintf(f_id, "module Hilbert (\n \
input  clk  ,\n \
input  data_ready  ,\n \
input  [31:0] data_in,\n \
output reg signed [31:0] Q_out\n \
);\n\n\
parameter SIZE = %d;\n\n", SIZE);

N = floor(length(g)/2);
for idx = 1:2:N
    fprintf(f_id, "parameter k_%d = %d ;\n", idx-1, round((2 ^ bits_hilbert) * -g(idx)));
end
fprintf(f_id, "\n");
for idx = 0:2:N-1
    fprintf(f_id, "reg signed [SIZE-1:0] sum_%d ;\n", idx);
end
fprintf(f_id, "\n");
for idx = 0:2:N-1
    fprintf(f_id, "reg signed [SIZE-1:0] mul_%d ;\n", idx);
end
fprintf(f_id, "\n");
for idx = 0:1:2*N-1
    fprintf(f_id, "reg signed [SIZE-1:0] v_%d;\n", idx);
end
fprintf(f_id, "\nalways @(posedge clk)\n \
begin\n");
fprintf(f_id, "   if (data_ready)\n \
      begin\n");
fprintf(f_id, "        sum_0 <= (v_0 - data_in);\n");
for idx = 2:2:N
    fprintf(f_id, "        sum_%d <= (v_%d - v_%d) ;\n", idx, idx, N*2 - idx);
end
for idx = 0:2:N-1
    fprintf(f_id, "        mul_%d <= sum_%d * k_%d; \n", idx, idx, idx);
end
%fprintf(f_id, "        Q_out <= sum_0 * k_0 +  \n");
%fprintf(f_id, "        0;\n");
for idx = 0:1:2*N-2
    fprintf(f_id, "        v_%d <= v_%d;\n", idx, idx+1);
end
fprintf(f_id, "        v_%d <= data_in;\n", 2*N-1);
fprintf(f_id, "      end\n");
fprintf(f_id, "   end\n");
fprintf(f_id, "endmodule\n");

fclose(f_id)

