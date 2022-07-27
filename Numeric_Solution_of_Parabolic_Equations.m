clear all
clc
tic
%% Numeric Solutions of Parabolic Equations
%% a) FTCS Explicit Method

% We can show the FTCS Explicit Method as ...
...(u_i^(n+1)-u_i^n)/delta_T= alpha(u_(i+1)^n-2u_i^n+u_(i-1)^n)/(delta_x)^2 ...

length = 0.04; % m
alpha = 0.000217; % m/sec^2
delta_t = 0.002; % sec
delta_x = 0.001; % m
Final_t = 1.08; % sec

Stability_Control = alpha*(delta_t)/(delta_x)^2; % Check whether the...
... solution is stabile or not.
if Stability_Control<=0.5
    fprintf('The Solution is stabile.\n')
else
    fprintf('The Solution is not stabile.\n')
end

n_x = length/delta_x +1; % The numbers of the point of delta_x
n_t = Final_t/delta_t +1; % The numbers of the point of delta_t
    
u(1:n_x,1:n_t) = zeros; %The solutions of FTCS, the matrix name is U.
u(1,1:n_t) = 40; % U, L = 0 m and initial conditions for t = 0-1.08 sec
u(n_x,1:n_t) = 0; % U, L = 0.04 m and initial conditions for t = 0-1.08 sec
u(2:(n_x-1),1) = 0; % U, each delta_x interval between 0<L<0.04 m...
...and initial conditions for t = 0
    
i = 2; %Initial value of matrix value
    
n=1; % Initial value for time interval

while n>=1 && n<=n_t-1
    i=2;
    while i>=2 && i<=n_x-1
    u(i,n+1) = u(i,n) + Stability_Control*(u(i+1,n) - 2*u(i,n) + u(i-1,n));
    i=i+1;
    end
n=n+1;
end
%% b) Laasonen Implicit Method

% We can show the Laasonen Implicit Method as...
...(u_i^(n+1)-u_i^n)/delta_T = ...
    ...alpha*(u_(i+1)^(n+1)-2u_i^(n+1)+u_(i-1)^(n+1))/(delta_x)^2 ...
    
delta_t_laasonen=0.01;
n_t_laasonen = Final_t/delta_t_laasonen + 1;
Stability_Control_Laasonen = alpha*(delta_t_laasonen)/(delta_x)^2;
C1=-Stability_Control_Laasonen;
    
C2=1+2*Stability_Control_Laasonen;

m(1:n_x,1:n_t_laasonen) = zeros; %The solutions of Laasonen...
...the matrix name is m.
m(1,1:n_t_laasonen) = 40; % U,Boundary conditions for L = 0 and t = 0-1.08 sec
m(n_x,1:n_t_laasonen) = 0; % U,Boundary conditions for L = 0.04 m and t = 0-1.08 sec
m(2:(n_x-1),1) = 0; % U, each delta_x interval between 0<L<0.04 m and...
...initial conditions for t = 0

k=2; %Initial value of matrix value
while k>=2 && k<=n_t_laasonen
    A = full(gallery('tridiag',n_x - 2,C1,C2,C1)); %diagonal matris
    b = repmat(m(2:(n_x-1),k-1),1); %b matrix
    b(1) = b(1)-C1*m(1,1); %Boundary conditions of b matrix for L=0
    b(n_x - 2) = b(n_x - 2) - C1*m(n_x,1); %Boundary conditions of b matrix...
    ... for L=0.04
    m(2:n_x-1,k) = inv(A*A')*A'*b;
    k = k + 1;
end
%% c)Crank-Nicolson Implicit Method

delta_t_CN = 0.01;
n_t_CN = Final_t/delta_t_CN + 1;
SC_CN = alpha*(delta_t_CN)/(delta_x)^2;
Cns1 = -SC_CN/2;
Cns2 = 1 + SC_CN;

z(1:n_x,1:n_t_CN) = zeros; %The solutions of Crank-Nicholson...
...the matrix name is z.
z(1,1:n_t_CN) = 40; % U,Boundary conditions for y = 0 and t = 0-1.08 sec
z(n_x,1:n_t_CN) = 0; % U,Boundary conditions for L = 0.04 m and t = 0-1.08 sec
z(2:(n_x-1),1) = 0; % U,each delta_x interval between 0<L<0.04 m and...
...initial conditions for t = 0

o=2;
while o>=2 && o<=n_t_CN
    AA = full(gallery('tridiag',n_x - 2,Cns1,Cns2,Cns1)); %diagonal matris
    bb=0;
    p=2;
    while p>=2 && p<=n_x-1
        pp(p-1,o-1) = -Cns1*(z(p-1,o-1) + z(p+1,o-1))+(1-SC_CN)*z(p,o-1);
        p = p+1;
    end
    bb = repmat(pp(1:(n_x-2),o-1),1); %bb matrix
    bb(1) = bb(1)-Cns1*z(1,1); % Boundary conditions of bb matrix for L=0
    bb(n_x - 2) = bb(n_x - 2) - Cns1*z(n_x,1); % Boundary conditions of...
    ...b matrix for L=0.04m
    z(2:n_x-1,o) = inv(AA*AA')*AA'*bb;
    o = o + 1;
end
%% Writing the data (Excel)
excel_deltax={num2cell(0:0.001:0.04)'}; %delta_x column
excel_deltat={num2cell(0:0.18:1.08)}; %delta_t row
excel_u=num2cell(u(1:41,1:90:541)); %FTCS value
ZZ=['delta_x/delta_t',excel_deltat{1};excel_deltax{1},excel_u];
xlswrite('FTCS_Part(a)',ZZ);

excel_m=num2cell(m(1:41,1:18:109));
excel_deltat_L={num2cell(0:0.18:1.08)}; %Laasonen values
TT=['delta_x/delta_t',excel_deltat_L{1};excel_deltax{1},excel_m];
xlswrite('Laasonen_Part(b)',TT);

excel_z=num2cell(z(1:41,1:18:109));
excel_deltat_Z={num2cell(0:0.18:1.08)}; %Crank-Nicholson values
OO=['delta_x/delta_t',excel_deltat_Z{1};excel_deltax{1},excel_z];
xlswrite('Crank-Nicholson_Part(c)',OO);
%% Grafikleri Çizdirme
figure('units','normalized','outerposition',[0 0 1 1]) 

subplot(2,2,1) %part a
plot(u(1:41,1:90:541),0:delta_x:length,'-v')
title('part a) Explicit Method')
xlabel('u(m/sec)')
ylabel('y(m)')
grid on
legend('t=0','t=0.18','t=0.36','t=0.54','t=0.72','t=0.90','t=1.08')

subplot(2,2,2) %part b
plot(m(1:41,1:18:109),0:delta_x:length,'-v')
title('part b) Laasonen Implicit Method')
xlabel('u(m/sec)')
ylabel('y(m)')
grid on
legend('t=0','t=0.18','t=0.36','t=0.54','t=0.72','t=0.90','t=1.08')

subplot(2,2,3) %part c
plot(z(1:41,1:18:109),0:delta_x:length,'-v')
title('part c) Crank-Nicholson Implicit Method')
xlabel('u(m/sec)')
ylabel('y(m)')
grid on
legend('t=0','t=0.18','t=0.36','t=0.54','t=0.72','t=0.90','t=1.08')
toc