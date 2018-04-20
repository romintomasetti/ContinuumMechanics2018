%% Compare the analytical and numerical results:

anal = XY;
lambda_anal = lambda_vec;

configureFigure(figure);
subplot(1,2,1);
hold on;
plot(anal(1:end-9,1),lambda_anal(1:end-9),'b-');
plot(UX2_SPH,LAMBDA_SPH,'r--')
xlabel('$u$')
ylabel('$\lambda$')
legend('N.-R.','SAL')
subplot(1,2,2);
hold on;
plot(anal(1:end-9,2),lambda_anal(1:end-9),'b-');
plot(spherical.uy2,spherical.lambda,'r--')
xlabel('$v$')
ylabel('$\lambda$')
legend('N.-R.','SAL')

