clear all
e=12;
An{1} = [0.009 0.027 0.303 0.066; -0.093 -0.033 0.021 0.195; 0.162 -0.048 -0.057 0.057; -0.06 -0.237 0.138 -0.288];
An{2} = [-0.315 0.294 0.252 0.156; -0.108 -0.249 0.255 -0.093; -0.186 -0.321 0.105 0.117; -0.177 -0.216 -0.033 0.129];
An{3} = [-0.048 -0.048 0.114 -0.261; 0.159 0.09 0.21 -0.006; 0.114 0.192 0.102 0.117; 0.111 0.087 0 0.24];
An{4} = [-0.081 0.198 0.009 0.198; 0.069 -0.189 0.087 0.033; 0.117 -0.114 -0.162 0.162; 0.033 0.048 0.168 -0.054];
B{1} = [1; 0; 0; 0];
C{1} = [0 0 0 1];
D{1} = 0;
Bw{1} = B{1};
%
% H = [2 3; 2 4];
% E = [1 1; 1 1];


N = size(An,2);
nx = size(An{1},1);
q = size(B{1},2);
for i=2:N
    B{i}=B{1};
    C{i}=C{1};
    D{i}=D{1};
    Bw{i}=Bw{1};
end
%% Declara��o das vari�veis da LMI

Beta=[];
Mu=[];
rho = 2.77;%:0.0:2.77;
for n2=1:length(rho)
    e=10;
    for bet=1:20
        beta=-bet;
        t=0;
        while (e>0.01)
            t=t+1;
            if t==5
                break
            end
            a=1;
            while a>0
                
                for i = 1:N
                    P{i} = sdpvar(nx,nx,'symmetric');
                end
                F1 = sdpvar(nx,nx,'full');
                W = sdpvar(q,nx,'full');
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
                for l=1:N
                    for i=1:N
                        J = [P{l}+F1+F1' -rho(n2)*An{i}*F1-B{i}*W zeros(nx,1) -Bw{i};
                            (-rho(n2)*An{i}*F1-B{i}*W)' -P{i} -F1'*C{i}'*beta zeros(nx,1);
                            zeros(nx,q)' (-F1'*C{i}'*beta)' 2*beta*diag(ones(q,1))+diag(ones(q,1)) -beta*D{i};
                            (-Bw{i})' zeros(nx,1)' (-beta*D{i})' -mu*diag(ones(q,1))];
                        %         M=[zeros(nx,nx) -E*F1 zeros(nx,1) zeros(nx,1)];
                        %         N1=[H' zeros(nx,nx) zeros(nx,1) zeros(nx,1)];
                        %         LMI1= [J M' N1';
                        %                M -gama*diag(ones(nx,1)) zeros(nx,nx);
                        %                N1 zeros(nx,nx) -gama*diag(ones(nx,1))];
                        LMIs = LMIs+(J<=0)+(P{l}>=0);
                    end
                end
                % solver
                options = sdpsettings('savesolveroutput',1,'verbose',1,'warning',1,'solver','sedumi','showprogress',0);
                %solucao=
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
            if teste >= 0
                Resp = 1; % flag de factibilidade
                
            else
                Resp = -1;
                break
            end
            a=1;
            
            while a>0
                for i = 1:N
                    P{i} = sdpvar(nx,nx,'symmetric');
                end
                F1=double(F1);
                G1=sdpvar(q,q);
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
                for l=1:N
                    for i=1:N
                        J2 = [P{l}+F1+F1' -rho(n2)*An{i}*F1-B{i}*W zeros(nx,1) -Bw{i};
                            (-rho(n2)*An{i}*F1-B{i}*W)' -P{i} -F1'*C{i}'*G1 zeros(nx,1);
                            zeros(nx,q)' (-F1'*C{i}'*G1)' 2*G1*diag(ones(q,1))+diag(ones(q,1)) -G1*D{i};
                            (-Bw{i})' zeros(nx,1)' (-G1*D{i})' -mu2*diag(ones(q,1))];
                        %         M=[zeros(nx,nx) -E*F1 zeros(nx,1) zeros(nx,1)];
                        %         N1=[H' zeros(nx,nx) zeros(nx,1) zeros(nx,1)];
                        %         LMI2= [J2 M' N1';
                        %                M -gama*diag(ones(nx,1)) zeros(nx,nx);
                        %                N1 zeros(nx,nx) -gama*diag(ones(nx,1))];
                        LMIs2 = LMIs2+(J2<=0)+(P{l}>=0);
                    end
                end
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
            if teste >= 0
                Resp = 1; % flag de factibilidade
                e=double(mu)-double(mu2);
                beta=double(G1);
                Beta=beta;
                Mu=double(mu2);
                rhomaior=rho(n2);
                P1 = double(P{1});
                P2 = double(P{2});
                P3 = double(P{3});
                P4 = double(P{4});
                Wf = double(W);
                F1f = double(F1);
            else
                Resp = -1;
                break
            end
            
        end
    end
end
w=double(W);
f=double(F1);
k=w/f;
