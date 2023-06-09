function EigenvalueSpectrum(freq,ev)
figure('Name','Eigenvalue Spectrum','WindowState','maximized');
subplot(3,2,1);
plot(freq,real(ev),'k.'); set(gca, 'FontSize', 13);
xlim([freq(1),freq(end)]);ylim([-15,15]);grid on;
xlabel('Frequency(Hz)');ylabel('Re(\lambda_n)');
subplot(3,2,2);
plot(freq,imag(ev),'k.'); set(gca, 'FontSize', 13);
xlim([freq(1),freq(end)]);ylim([-15,15]);grid on;
xlabel('Frequency(Hz)');ylabel('Im(\lambda_n)');
subplot(3,2,3);
plot(freq,abs(real(ev)),'k.'); set(gca, 'YScale', 'log', 'FontSize', 13);
xlim([freq(1),freq(end)]);ylim([1e-2,1e4]);grid on;
xlabel('Frequency(Hz)');ylabel('|Re(\lambda_n)|');
subplot(3,2,4);
plot(freq,abs(imag(ev)),'k.'); set(gca, 'YScale', 'log', 'FontSize', 13);
xlim([freq(1),freq(end)]);ylim([1e-2,1e4]);grid on;
xlabel('Frequency(Hz)');ylabel('|Im(\lambda_n)|');
subplot(3,2,5);
plot(freq,1./abs(1+1i*ev),'k.'); set(gca, 'FontSize', 14);
xlim([freq(1),freq(end)]);ylim([0,1]);grid on;
xlabel('Frequency(Hz)');ylabel('MS_n');
subplot(3,2,6);
plot(freq,180 - rad2deg(atan(real(ev))),'k.'); set(gca, 'FontSize', 13);
xlim([freq(1),freq(end)]);ylim([90,270]);grid on;
xlabel('Frequency(Hz)');ylabel('\alpha_n');
end
