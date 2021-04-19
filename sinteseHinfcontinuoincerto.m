clear all
e=10;
k1 = 1; k2 = 1; m1 = 1; m2 = 0.5; c0 = 2;
An{1} = [0 0 1 0; 0 0 0 1; -(k1+k2)/m1  k2/m1 -c0/m1 0; k2/m2 -k2/m2 0 -c0/m2];
An{2} = [0 0 1 0; 0 0 0 1; -(k1+k2)/m1  k2/m1 -c0/m1 0; k2/m2  k2/m2 0 -c0/m2];
An{3} = [0 0 1 0; 0 0 0 1; -(k1+k2)/m1 -k2/m1 -c0/m1 0; k2/m2  k2/m2 0 -c0/m2];
An{4} = [0 0 1 0; 0 0 0 1; -(k1+k2)/m1 -k2/m1 -c0/m1 0; k2/m2 -k2/m2 0 -c0/m2];
B{1} = [0 0; 0 0; -1/m1 0; 0 -1/m2];
C{1} = [0 0 1 0; 0 0 0 1];
D{1} = [0 ; 0];
Bw{1} = [0; 0; 1/m1; 0];
H = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 1 0];
E = [0 0 0 0; 0 0 0 0; 0.6/m2 -0.6/m2 0 -1.97/m2; -1.2/m1 0.6/m1 -1.97/m1 0];
mumelhor=200;
N = size(An,2);
nx = size(An{1},1);
q = size(B{1},2);
alfa=94:2:98;
for i=2:N
    B{i}=B{1};
    C{i}=C{1};
    D{i}=D{1};
    Bw{i}=Bw{1};
end
%% Declara��o das vari�veis da LMI
Beta=zeros(1,length(alfa));
Mu=zeros(1,length(alfa));
%% Declara��o das vari�veis da LMI
for n = 1:length(alfa)
    e=10;
    for bet=1:20
        beta=-bet;
        while (e>0.0001)
            a=1;
            
            while a>0
                P = sdpvar(nx,nx,'symmetric');
                F1 = sdpvar(nx,nx,'full');
                W = sdpvar(q,nx,'full');
                %G3 = sdpvar(nx,nx,'full');
                mu = sdpvar(1,1);
                gama=sdpvar(1,1);
                
                if a==1
                    obj = mu;
                else
                    mu = Xf;
                    obj = [];
                end
                %% LMIs programadas
                LMIs = [];
                for i = 1:N
                    J = [F1+F1'                          P-An{i}*F1-B{i}*W+alfa(n)*F1'                     zeros(nx,q)                           -Bw{i};
                        (P-An{i}*F1-B{i}*W+alfa(n)*F1')'     alfa(n)*(-An{i}*F1-B{i}*W-F1'*An{i}'-W'*B{i}')              -F1*C{i}'*beta                            -alfa(n)*Bw{i};
                        zeros(nx,q)'                       (-F1*C{i}'*beta)'                 2*beta*diag(ones(q,1))+diag(ones(q,1))        -beta*D{i};
                        (-Bw{i})'                             (-alfa(n)*Bw{i})'                          (-beta*D{i})'                              -mu*diag(ones(1,1))];
                    M=[zeros(nx,nx) -E*F1 zeros(nx,q) zeros(nx,1)];
                    N1=[H' alfa(n)*H' zeros(nx,q) zeros(nx,1)];
                    LMI1= [J M' N1';
                        M -gama*diag(ones(nx,1)) zeros(nx,nx);
                        N1 zeros(nx,nx) -gama*diag(ones(nx,1))];
                    LMIs = LMIs+(LMI1<=0)+(P>=0)+(abs(gama)<=1);
                end
                % solver
                options = sdpsettings('savesolveroutput',1,'verbose',1,'warning',1,'solver','sedumi','showprogress',0);
                %soluC{i}ao=
                % Tentativa de solu��o
                tempo = cputime;
                solvesdp(LMIs,obj,options);
                tempo = cputime-tempo;
                %% Verifica a necessidade de outro loop
                teste = min(checkset(LMIs));
                if teste <0 && a==1
                    Xf = 1.001*double(mu);
                    a=a+1;
                else
                    break
                end
            end
            if teste > 0
                Resp = 1; % flag de factibilidade
                
            else
                Resp = -1;
                break
            end
            a=1;
            
            while a>0
                G1=sdpvar(1,1);
                P = sdpvar(nx,nx,'symmetric');
                F1 = double(F1);
                W = sdpvar(q,nx,'full');
                mu2 = sdpvar(1,1);
                gama=sdpvar(1,1);
                if a==1
                    obj = mu2;
                else
                    mu2 = Xf;
                    obj = [];
                end
                %% LMIs programadas
                LMIs2 = [];
                J2 = [F1+F1'                         P-An{i}*F1-B{i}*W+alfa(n)*F1'                     zeros(nx,q)                     -Bw{i};
                    (P-An{i}*F1-B{i}*W+alfa(n)*F1')'     alfa(n)*(-An{i}*F1-B{i}*W-F1'*An{i}'-W'*B{i}')                -F1*C{i}'*G1                      -alfa(n)*Bw{i};
                    zeros(nx,q)'                       (-F1*C{i}'*G1)'                2*G1*diag(ones(q,1))+diag(ones(q,1))       -G1*D{i};
                    (-Bw{i})'                               (-alfa(n)*Bw{i})'                            (-G1*D{i})'                      -mu2*diag(ones(1,1))];
                M=[zeros(nx,nx) -E*F1 zeros(nx,q) zeros(nx,1)];
                N1=[H' alfa(n)*H' zeros(nx,q) zeros(nx,1)];
                LMI2= [J2 M' N1';
                    M -gama*diag(ones(nx,1)) zeros(nx,nx);
                    N1 zeros(nx,nx) -gama*diag(ones(nx,1))];
                LMIs2 = LMIs2+(LMI2<=0)+(P>=0)+(abs(gama)<=1);
                % solver
                options = sdpsettings('savesolveroutput',1,'verbose',1,'warning',1,'solver','sedumi','showprogress',0);
                %solucao=
                % Tentativa de solu��o
                tempo = cputime;
                solvesdp(LMIs2,obj,options);
                tempo = cputime-tempo;
                %% Verifica a necessidade de outro loop
                teste = min(checkset(LMIs2));
                if teste <0 && a==1
                    Xf = 1.001*double(mu2);
                    a=a+1;
                else
                    break
                end
            end
            if teste > 0
                Resp = 1; % flag de factibilidade
                
            else
                Resp = -1;
                break
            end
            
            e=abs(double(mu)-double(mu2));
            beta=double(G1);
            Beta(n)=beta;
            Mu(n)=double(mu2);
            if(double(mu2)<mumelhor)
                mumelhor=double(mu2);
                Pmelhor=double(P);
                Wmelhor=double(W);
                F1melhor=F1;
                gamamelhor=double(gama);
                G1melhor=double(G1);
                alfamelhor=alfa(n);
            end
        end
    end
end
%%

w=double(W);
f=double(F1);
k=w/f;
Hinft = sqrt(mumelhor)
%%