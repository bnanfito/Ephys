function [X,f,P,Nt,phaseAngle]=myfft(x,dtORfs,plotSpec)
% [X,f,P,Nt]=myfft(x,dtORfs,plotSpec)
%  
% Inputs
% x       : time domain sigmal
% dtORfs  : sampling period or frequency. Currently assumes that dtORfs <=1
%           is a sampling period, whereas dtORfs >1 is a sampling frequency. 
% plotSpec: =1 => show plot; not specified, or =0 => don't plot
%
% Outputs
% X : fft(x)
% f : frequency range (Hz)
% P : power spectrum of x
% Nt: time stamp for x (based on your input dtORfs).
%
% Example: 
% dt=0.01; t=0:dt:5;
% x=7*sin(2*pi*10*t) + 19*sin(2*pi*23.5*t);
% [X,f,P,Nt]=myfft(x,dtORfs,1);
% 
% Note in this example that the multiplier of t inside the sin function (i.e,. 2*pi*10)
% is frequency in radians/sec (omega). Recall that f = 2*pi*omega. f is in Hz.
% So in the example above, the two frequencies are 10Hz or 2*pi*10 rad/s;
% and 23.5Hz or 2*pi*23.5 rad/s
%
% --------------------
% Shreesh P. Mysore
% shreesh@stanford.edu
% 2.10.07
% --------------------

xhat=[];
if nargin<3, plotSpec=0; end

%-- dt or fs?
if dtORfs > 1, Fs = dtORfs; dt = 1/Fs; else dt = dtORfs; Fs =1/dt; end 
% min sampling frequency is 1Hz just as a way to identify what is being inputted.

%-- fft
N=length(x); Nt = dt*(1:N);
Nfft = 2^nextpow2(N); %N=2^round(log2(N));
%X=fft(x,Nfft)/N;Pxx=2*abs(X);P=Pxx(1:Nfft/2);

X=fft(x,N)/N;Pxx=2*abs(X);
%<==> 
%X = fft(x);Pxx = 2*sqrt(X.*conj(X))/N; %Pxx = X.*conj(X)/N;

%-- power spectrum %magnitude of complex Fourier coefficients)
P=Pxx(1:floor(N/2)+1);
%P=0.5*Pxx; % all points in Pxx

% %-- phase of coefficients
phaseAngle = unwrap (angle (X(1:floor(N/2)+1))); %in radians
phaseAngle = phaseAngle*180/pi; %in degrees


%-- frequency range
f = Fs/2*linspace(0,1,Nfft/2);
f = Fs/2*linspace(0,1,floor(N/2)+1);
f = ((0: 1/(floor(N/2)+1): 1-1/(floor(N/2)+1))*Fs/2).'; 
%[length(f) N]
% pause
%f = ((0: 1/N: 1-1/N)*Fs).'; %all points in Pxx 


%-- plotting
if plotSpec==1
    figure
    %  ind = find(f>0 & f~=10); P(ind) = zeros(size(ind));
    %  X(ind) = zeros(size(ind));
    %  X = [X(1:Nfft/2) fliplr(X(1:Nfft/2))];
    %xhat = ifft(N*X);
    
    subplot(3,1,1), plot(Nt,x);xh=xlabel('time (seconds)');
    set(gca,'fontsize',12)
    yh=ylabel('x'); set(xh,'fontsize',15);set(yh,'fontsize',15)
    axis tight;
    
    subplot(3,1,2), plot(f,P,'.-'), 
    set(gca,'fontsize',12)
    xh=xlabel('frequency (Hz)'); yh=ylabel('| fft(x) |'); %th = title('Frequency spectrum of X');
    set(xh,'fontsize',15);set(yh,'fontsize',15); %set(th,'fontsize',20);
    %subplot(3,1,3), plot(Nt,xhat(1:N));title('x'), xlabel('time (seconds)')
    set(gcf,'color','w')
    
%     subplot(3,1,3), 
%     plot(f,phaseAngle,'-'), 
%     set(gca,'fontsize',12); set(gca,'ylim',[-180 180])
%     xh=xlabel('frequency (Hz)'); yh=ylabel('phase (deg)');
%     set(xh,'fontsize',15);set(yh,'fontsize',15); 
%     %subplot(3,1,3), plot(Nt,xhat(1:N));title('x'), xlabel('time (seconds)')
%     set(gcf,'color','w')
   
end

