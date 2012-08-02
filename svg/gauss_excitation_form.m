%
% this script evaluates the same gaussian excitation function, as openEMS does
%

clear
close all
clc

f0 = 10e9;
fc = 6e9;
dT = 1e-12; % sample time-step


sigma = 1/sqrt(8/9)/pi/fc;
t0 = sqrt(18)/sqrt(8/9)/pi/fc;

len = 2 * 9/(2*pi*fc) / dT; % gauss length

for n=1:len
    t_(n) = (n-1)*dT;
    ex(n) = cos(2*pi*f0*((n-1)*dT - 9/(2*pi*fc))) .* exp(-((t_(n)-t0)/sigma)^2/2);
end

%%  time domain
plot(t_/1e-9,ex)
title('Gaussian pulse in time domain','fontsize',14,'fontweight','b','color','k');
xlabel( 'time (ns)','fontsize',12,'fontweight','b','color','k' );
ylabel( 'amplitude','fontsize',12,'fontweight','b','color','k' );
 hold on;
 
 % t0 line
 x_posi=t0/1e-9;
 yLim=get(gca,'YLim');
                y_posi=yLim(2);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               hold on;
               y_posi=yLim(1);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
disp( ['Amplitude at t=0: ' num2str(20*log10(abs(ex(1))/1)) ' dB'] );
text(x_posi,y_posi*0.9,'t_0','fontsize',12,'fontweight','b','color','k');
grid on
%% freq domain
val = DFT_time2freq( t_, ex, [f0-fc f0 f0+fc] );
disp( ['Amplitude at f=f0-fc: ' num2str(20*log10(abs(val(1))/abs(val(2)))) ' dB'] );
disp( ['Amplitude at f=f0+fc: ' num2str(20*log10(abs(val(3))/abs(val(2)))) ' dB'] );

% calculate frequency domain via slow DFT
freq = linspace(f0-fc,f0+fc,1000);
val = DFT_time2freq( t_, ex, freq );
figure
plot( freq/1e9, abs(val) )

hold on

x_posi=f0/1e9;
 yLim=get(gca,'YLim');
                y_posi=yLim(2);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               text(x_posi,y_posi*0.02,'f_0','fontsize',12,'fontweight','b','color','k');
               
               x_posi=(f0-fc/2)/1e9;
               y_posi=yLim(1);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               
               hold on;
               y_posi=yLim(2);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               text(x_posi,y_posi*0.02,'f_0-f_c/2','fontsize',12,'fontweight','b','color','k');
               hold on;
               
 x_posi=(f0+fc/2)/1e9;
               y_posi=yLim(1);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               
               hold on;
               y_posi=yLim(2);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               text(x_posi,y_posi*0.02,'f_0+f_c/2','fontsize',12,'fontweight','b','color','k');
               hold on;
% % overlay the FFT result
% [f,val_fft] = FFT_time2freq( t_, ex );
% val_fft = val_fft((f0-fc<=f) & (f<=f0+fc));
% f = f((f0-fc<=f) & (f<=f0+fc));
% hold on
% plot( f/1e9, abs(val_fft), 'r' )
title('Gaussian pulse in frequency domain','fontsize',14,'fontweight','b','color','k');
              
if (f0==0)
    Fw = sigma*sqrt(2*pi)*exp(-0.5*(sigma*2*pi*f).^2);
    plot( f/1e9, 2*abs(Fw), 'g--' )
%     legend('dft','fft','analytic')
else
%     legend('dft','fft')
end

xlim([0 1.1*max(freq)/1e9])

xlabel( 'frequency (GHz)','fontsize',12,'fontweight','b','color','k' );
ylabel( 'amplitude','fontsize',12,'fontweight','b','color','k' );
grid on

%% dB
figure
val = val(freq>=0);
freq = freq(freq>=0);
plot( freq/1e9, 20*log10(abs(val)/max(abs(val))))
xlabel( 'frequency (GHz)','fontsize',12,'fontweight','b','color','k' );
ylabel( 'amplitude (dB)','fontsize',12,'fontweight','b','color','k' );
title('Amplitude (dB) of Gaussian pulse in frequency domain','fontsize',14,'fontweight','b','color','k');
xlim([0 1.1*max(freq)/1e9]);
ylim([min(20*log10(abs(val)/max(abs(val)))) max(20*log10(abs(val)/max(abs(val))))+1]);
hold on;

x_posi=f0/1e9;
 yLim=get(gca,'YLim');
               y_posi=yLim(1);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               text(x_posi,y_posi*0.96,'f_0','fontsize',12,'fontweight','b','color','k');
               hold on;
               y_posi=yLim(2);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               hold on;

x_posi=(f0-fc/2)/1e9;
               y_posi=yLim(1);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               text(x_posi,y_posi*0.96,'f_0-f_c/2','fontsize',12,'fontweight','b','color','k');
               hold on;
               y_posi=yLim(2);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               hold on;
               
 x_posi=(f0+fc/2)/1e9;
               y_posi=yLim(1);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               text(x_posi,y_posi*0.96,'f_0+f_c/2','fontsize',12,'fontweight','b','color','k');
               hold on;
               y_posi=yLim(2);
               stem(x_posi,y_posi,'--.','Color',[0.7,0.7,0.7]);
               hold on;
   plot([0 1.1*max(freq)/1e9], [-5 -5],'--.','Color',[0.7,0.7,0.7])            
               grid on
               text(-0.7,-5,'-5','fontsize',12,'fontweight','b','color','k');