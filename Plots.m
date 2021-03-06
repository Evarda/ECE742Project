figure
semilogy((1:nmax)*dt,abs(Error(2:6,:)-Error(1,:))/abs(max(Error(1,:))));
legend('\sigma_{max}=10^{-1}','\sigma_{max}=10^{0}','\sigma_{max}=10^{1}','\sigma_{max}=10^{2}','\sigma_{max}=10^{3}','location','southeast');
hold off
grid minor
grid on
title('Relative Error over Time for Varying Values of \sigma_{max}, m = 2, \kappa_{max} = 3')
xlabel('Time (s)');
ylabel('Relative Error')
axis([0, 2.1e-9, 1e-10, 10])

figure
semilogy((1:nmax)*dt,abs(Error(7:12,:)-Error(1,:))/abs(max(Error(1,:))));
legend('m = 0','m = 1','m = 2','m = 3','m = 4','m = 5','location','southeast');
hold off
grid minor
grid on
title('Relative Error over Time for Varying Values of m, \kappa_{max} = 1, \sigma_{max} = Optimal')
xlabel('Time (s)');
ylabel('Relative Error')
axis([0, 2.1e-9, 1e-10, 10])

figure

semilogy((1:nmax)*dt,abs(Error(13:17,:)-Error(1,:))/abs(max(Error(1,:))));
legend('\sigma_{max}=10^{-1}','\sigma_{max}=10^{0}','\sigma_{max}=10^{1}','\sigma_{max}=10^{2}','\sigma_{max}=10^{3}','location','southeast');
grid minor;
grid on;
title('Relative Error over Time for Varying Values of \sigma_{max}, m = 3, \kappa_{max} = 1')
xlabel('Time (s)');
ylabel('Relative Error');
axis([0, 2.1e-9, 1e-10, 10])