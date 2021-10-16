%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Modelo de Crescimento Deterministico
%Baseado em Miao (2014) com atualizações para evitar comandos
%descontinuados
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;

%--------------------------------------------------------------------------
%Preambulo
%--------------------------------------------------------------------------

var c n k z;
varexo e;
parameters beta chi delta alpha rho;

%--------------------------------------------------------------------------
%Calibração
%--------------------------------------------------------------------------

alpha=0.33;
beta=0.99;
delta=0.023;
chi=1.75;
rho=0.95;

%--------------------------------------------------------------------------
%modelo
%--------------------------------------------------------------------------

model;
  (1/c) = (beta/c(+1))*(1-delta+alpha*exp(z(+1))*k^(alpha-1)*n(+1)^(1-alpha));
  (chi*c/(1-n)) = (1-alpha)*exp(z)*(k(-1)^(alpha))*(n^(-alpha));
  c+k-(1-delta)*k(-1) = exp(z)*(k(-1)^(alpha))*(n^(1-alpha));
  z = rho*z(-1)+e;
end;

%--------------------------------------------------------------------------
%Estados estacionário e inicial e choques
%--------------------------------------------------------------------------

initval;
  k = 9;
  c = 0.76;
  n = 0.3;
  z = 0;
  e = 0;
end;
steady;
check;

%Choque temporário nos 10 primeiros períodos
shocks;
  var e ; periods 1:10;
  values 0.01;
end;

%Choque permanente de 0.01 começando no período 10
%endval;
% k = 9;
%  c = 0.76;
%  n = 0.3;
%  z = 0;
%  e = 0.01;
%end;
%steady;

%shocks;
%  var e ; periods 1:9;
%  values 0;
%end;

%choque temporário de 0.015 dos períodos 5 a 10
% shocks;
%  var e ; periods 5:10;
%  values 0.015;
%end;

%choque permanente de 0.015
%endval;
% k = 9;
%  c = 0.76;
%  n = 0.3;
%  z = 0;
%  e = 0.015;
%end;
%steady;

%--------------------------------------------------------------------------
%Simulação
%--------------------------------------------------------------------------

perfect_foresight_setup(periods=200);
perfect_foresight_solver;

%--------------------------------------------------------------------------
%Gráficos disponíveis no original(consumo, trabalho, capital e tecnologia)
%--------------------------------------------------------------------------

tt=1:202';

subplot(2,2,1);
plot(tt,oo_.endo_simul(1,:),tt,oo_.steady_state(1)*ones(1,202),'LineWidth',2);
title('C')

%O código disponível em Miao erroneamente inverte K e N nos gráficos
subplot(2,2,2);
plot(tt,oo_.endo_simul(2,:),tt,oo_.steady_state(2)*ones(1,202),'LineWidth',2);

title('N')

subplot(2,2,3);
plot(tt,oo_.endo_simul(3,:),tt,oo_.steady_state(3)*ones(1,202),'LineWidth',2);
title('K');


subplot(2,2,4);
plot(tt,oo_.endo_simul(4,:),'LineWidth',2);
title('z')

%--------------------------------------------------------------------------
%Gráficos adicionais (produto, investimento, salário e aluguel do capital)
%--------------------------------------------------------------------------

y = exp(z(2:202)).*(n(2:202).^(1-alpha)).*(k(1:201).^(alpha));
y_bar = exp(oo_.steady_state(4))*(oo_.steady_state(3)^(alpha))*(oo_.steady_state(2)^(1-alpha));

I = y-c(2:202);
I_bar = y_bar-oo_.steady_state(1);

w = (1-alpha)*exp(z(2:202)).*(n(2:202).^(-alpha)).*(k(1:201).^(alpha));
w_bar = (1-alpha)*exp(oo_.steady_state(4))*(oo_.steady_state(3)^(alpha))*(oo_.steady_state(2)^(-alpha));

r = alpha*exp(z(2:202)).*(n(2:202).^(1-alpha)).*(k(1:201).^(alpha-1));
r_bar = alpha*exp(oo_.steady_state(4))*(oo_.steady_state(3)^(alpha-1))*(oo_.steady_state(2)^(1-alpha));

figure();
tt2=1:201';

subplot(2,2,1);
plot(tt2,y,tt2,y_bar*ones(1,201),'LineWidth',2);
title('Y')

subplot(2,2,2);
plot(tt2,I,tt2,I_bar*ones(1,201),'LineWidth',2);

title('I')

subplot(2,2,3);
plot(tt2,w,tt2,w_bar*ones(1,201),'LineWidth',2);
title('W');

subplot(2,2,4);
plot(tt2,r,tt2,r_bar*ones(1,201),'LineWidth',2);
title('R')