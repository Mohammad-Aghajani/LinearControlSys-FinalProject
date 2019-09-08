%Parameters from Table 2.1
l     = 1.2     ; % Length of link
Gamma = 1.2     ; % Gamma = pA
b     = 0.59    ; % b = F1;
EI    = 1.94    ; % = E*I;
w1    = 3       ; % Frequency of delta1
w2    = 19      ; % Frequency of delta2
c1    = 0.4     ; % Constant in phi1 formula that normalizes the eigenfunction
c2    = 4       ; % Constant in phi2 formula that normalizes the eigenfunction
a     = 0.6     ; % The ratio for stabilizing the system
f1    = 2.5   ; % First element of F2 matrix
f2    = 9    ; % Second element of F2 matrix
k1    = 17.4561 ; % First element of K matrix
k2    = 685.5706; % Second element of K matrix

%Defining M matrix from appendix B
M = [ 0.9929   1.0703  -0.0282;
      1.0703   1.6235  -0.4241;
     -0.0282  -0.4241   2.5920];
 
Minv = inv(M);

%Defining KK matrix
KK = [0 0  0;
      0 k1 0;
      0 0  k2];
FF = [b 0  0;
      0 f1 0;
      0 0  f2];
 % Definig A matrix
 A = [zeros(3,3) eye(3,3);
      -Minv*KK   -Minv*FF];
 % Defining B matrix
 P = [1;0;0];
 B = [zeros(3,1);
      Minv*P];
 
 %For calculating C matrix:
 beta1 = sqrt(w1*(sqrt(Gamma*l^4/EI)));
 beta2 = sqrt(w2*(sqrt(Gamma*l^4/EI)));
 phi1  = @(x) c1*((sin(beta1*x)-sinh(beta1*x)) - (sin(beta1*l)+sinh(beta1*l))*(cos(beta1*x)-cosh(beta1*x))/(cos(beta1*l)+cosh(beta1*l)));
 phi2  = @(x) c2*((sin(beta2*x)-sinh(beta2*x)) - (sin(beta2*l)+sinh(beta2*l))*(cos(beta2*x)-cosh(beta2*x))/(cos(beta2*l)+cosh(beta2*l)));
 
 phi1L = phi1(l);
 phi2L = phi2(l);
 %Defining C matrix
 C = [1 a*phi1L/l a*phi2L/l 0 0 0];
 D = 0;
 %Convert from state space to transfer function
 [numCo,denumCo] = ss2tf(A,B,C,D);
 G = tf(numCo,denumCo);
 %sisotool(G)


 
 
 
 
 
 
 
 
 