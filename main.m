%This Matlab file is done based on the Refereance (Polyakov, Andrey, Denis Efimov, and Wilfrid Perruquetti.
%"Finite-time and fixed-time stabilization: Implicit Lyapunov function approach." Automatica 51 (2015): 332-340).
A=[3 1 2;4 2 1;-2 -1 -1];%System matrix
B=[0;2;1];%input Matrix
rank([B,A*B,A*A*B])%Check that the system is controllable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Decomposition method: Step 1: Find (G)
%inialization
A0=A;
B0=B;
T0=eye(3);
% itieration 1
B0_orth=null(B0')';
B0_tild=null(B0_orth)';
A1=B0_orth*A0*B0_orth';
B1=B0_orth*A0*B0_tild';
T1=[B0_orth;B0_tild];
rank(B1)<size(A1,1)%check the condition to end the algorithm
% itieration 2
B1_orth=null(B1')';
B1_tild=null(B1_orth)';
A2=B1_orth*A1*B1_orth';
B2=B1_orth*A1*B1_tild';
T2=[B1_orth;B1_tild];
rank(B2)<size(A2,1)%check the condition to end the algorithm
%%%Calcuation of transformation matrix G
G=blkdiag(T2,1)*T1;%Transformation matrix
Astep1=G*A*G';%System matrix after trasfomation is applied
Bstep1=G*B;%Input matrix after trasfomation is applied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Decomposition method: Step 2: Find (PHI)
A11=Astep1(1,1);
A21=Astep1(2,1);
A22=Astep1(2,2);
A12=Astep1(1,2);
A23=Astep1(2,3);
A34=Bstep1(3);

A12_plus=A12'*inv(A12*A12');
A23_plus=A23'*inv(A23*A23');
A34_plus=A34'*inv(A34*A34');

PHI=eye(3);
PHI(2,1)=A12_plus*A11;
PHI(3,1)=A23_plus*(A21+A12_plus*A11*A11);
PHI(3,2)=A23_plus*(A22+A12_plus*A11*A12);
%System matrix after step 1 and step 2
Anew=PHI*G*A*G'*inv(PHI);
%Input matrix after step 1 and step 2
Bnew=PHI*G*B;
%Calculation for the gain to transfer the system to a cascaded integrator
%form
Kline=Anew(3,:);
Anew-Bnew*A34_plus*Anew(3,:)%Check that the system will be in cascaded integartor form
%% Calculation for the finite time stabilization control parameters (K,P)
% NOTE: The Library "YALMIP" is needed to evaluate this section
% Here we condiser only the cascaded integrator system after applyiing
% using the gain Kline.
clear F
A_int=zeros(3,3);
A_int(1,2)=Anew(1,2);
A_int(2,3)=Anew(2,3);
B_int=[0;0;1];
%The value chosen to solve the LMIs
mui=0.5;
alpha=1;
gamma=4.5;
X=sdpvar(3,3);
y=sdpvar(1,3);
H=[-1-2*mui 0 0;0 -1-1*mui 0;0 0 -1];
F=[X>=0];
F=F+[(gamma*X+X*H+H*X)>=0];
F=F+[(X*H+H*X)<=0];
F=F+[(A_int*X+X*A_int'+B_int*y+y'*B_int'+alpha*X)<=0];
solvesdp(F);

X=double(X);
y=double(y);
%Check that the solution is correct
min(eig(X))%should be positive
min(eig(gamma*X+X*H+H*X))%should be positive
max(eig(X*H+H*X))%should be negative
max(eig(A_int*X+X*A_int'+B_int*y+y'*B_int'+alpha*X))%should be negative
% Gain matrix and Matrix P is to define the implicitly defined Lypunov funcition 
K=y/X;
P=inv(X);
%% Find the time expected for stabilization
x0=[-1;5;3];
s_int=PHI*G*x0;%the value of the s corresponds to inial state variable 
Vint=fminsearch(@(V)(s_int'*[V^(-1-2*mui) 0 0;0 V^(-1-1*mui) 0;0 0 V^(-1)]*P*[V^(-1-2*mui) 0 0;0 V^(-1-1*mui) 0;0 0 V^(-1)]*s_int-1)^2,0.1);
T=gamma*Vint^mui/(alpha*mui)%maximum expected time for stabilization
