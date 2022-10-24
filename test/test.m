clear;
clc;

m = 11.5;
W = 112.8;
B = 114.8;

zg = 0.02;

Ix = 0.16;
Iy = 0.16;
Iz = 0.16;

Xu_dot = -5.5;
Yv_dot = -12.7;
Zw_dot = -14.57;
Kp_dot = -0.12;
Mq_dot = -0.12;
Nr_dot = -0.12;

Xu = -4.03;
Yv = -6.22;
Zw = -5.18;
Kp = -0.07;
Mq = -0.07;
Nr = -0.07;

Xu_abs = -18.18;
Yv_abs = -21.66;
Zw_abs = -36.99;
Kp_abs = -1.55;
Mq_abs = -1.55;
Nr_abs = -1.55;


MRB = [m  0  0  0  m*zg  0;
       0  m  0  -m*zg  0  0;
       0  0  m  0  0  0;
       0  -m*zg  0  Ix  0  0;
       m*zg  0  0  0  Iy  0;
       0  0  0  0  0  Iz];

MA = (-1).*[Xu_dot  0  0  0  0  0;
            0  Yv_dot  0  0  0  0;
            0  0  Zw_dot  0  0  0;
            0  0  0  Kp_dot  0  0;
            0  0  0  0  Mq_dot  0;
            0  0  0  0  0  Nr_dot];

M=MRB+MA;

M_inv = inv(MRB+MA);



N=10000;My=6;Mu=6;
yd=cell(1);
for k=1:N+1
    yd{k}(1)=sin(0.001*k);
    yd{k}(2)=sin(0.001*k);
    yd{k}(3)=-0.01*k;
    yd{k}(4)=0;
    yd{k}(5)=0;
    yd{k}(6)=0;
end

%控制器阶数
Lu=3;
Ly=1;
eta=0.5;
miu=1;
rou=0.5;
lamda=1;
M=50;
epsilon=10^(-5);
%初始值
y=cell(1);u=cell(1);du=cell(1);
y{1}=[0,0,0,0,0,0];
y{2}=y{1};y{3}=y{2};y{4}=y{3};
u{1}=zeros(Mu,1);u{2}=zeros(Mu,1);
u{3}=u{1};u{4}=u{3};
%控制器伪偏导数初始值
Fai=cell(1);
Fai1=cell(1);
Fai2=cell(1);
Fai3=cell(1);
Fai4=cell(1);

Fai1{1}=diag([0.5,0.5,0.5,0.5,0.5,0.5],0);
Fai2{1}=zeros(6);
Fai3{1}=zeros(6);
Fai4{1}=zeros(6);
Fai{1}=[Fai1{1}(:,:) Fai2{1}(:,:) Fai3{1}(:,:) Fai4{1}(:,:)];
Fai{2}=Fai{1};
Fai{3}=Fai{1};
Fai{4}=Fai{1};
zerou=zeros(Mu,1);zeroy=zeros(My,1);
%程序循环
% for k=5:N
%     %计算dH
% %     dy = y{k-1}(:)-y{k-2}(:);
% %     du{1}=u{k-1}(:)-u{k-2}(:);
% %     du{2}=u{k-2}(:)-u{k-3}(:);
% %     du{3}=u{k-3}(:)-u{k-4}(:);
% %     dH=[dy; du{1}(:); du{2}(:); du{3}(:)];
% 
%     if k==2
%         dH=u{k-1}(:)-0;
%     else
%         dH=u{k-1}(:)-u{k-2}(:);
%     end
%     for i=2:Lu
%         if k>i+1
%             dH=[dH;u{k-i}(:)-u{k-i-1}(:)];
%         else
%             dH=[dH;zerou];
%         end
%     end
%     for i=1:Ly
%         if k>i+1
%             dH=[dH;y{k-i}(:)-y{k-i-1}(:)];
%         else
%             dH=[dH;zeroy];
%         end
%     end
% 
%     %计算并更新Fai
%     Fai{k}=Fai{k-1}+eta*(dy-Fai{k-1}*dH)*dH'/(miu+norm(dH,2)^2);
%     Fai1{k} = Fai{k}(1:6,1:6);
%     Fai2{k} = Fai{k}(1:6,7:12);
%     Fai3{k} = Fai{k}(1:6,13:18);
%     Fai4{k} = Fai{k}(1:6,19:24);
%     for i=1:6
%         for j=1:6
%             if i==j
%                 if (Fai1{k}(i,j)<epsilon)|(Fai1{k}(i,j)>M)
%                     Fai1{k}(i,j)=Fai1{1}(i,j);
%                 end
%                 if (Fai2{k}(i,j)<epsilon)|(Fai2{k}(i,j)>M)
%                     Fai2{k}(i,j)=Fai2{1}(i,j);
%                 end
%                 if (Fai3{k}(i,j)<epsilon)|(Fai3{k}(i,j)>M)
%                         Fai3{k}(i,j)=Fai3{1}(i,j);
%                 end
%                 if (Fai4{k}(i,j)<epsilon)|(Fai4{k}(i,j)>M)
%                         Fai4{k}(i,j)=Fai4{1}(i,j);
%                 end
%             else
%                 if (Fai1{k}(i,j)<epsilon)|(Fai1{k}(i,j)>M)
%                     Fai1{k}(i,j)=Fai1{1}(i,j);
%                 end
%                 if (Fai2{k}(i,j)<epsilon)|(Fai2{k}(i,j)>M)
%                     Fai2{k}(i,j)=Fai2{1}(i,j);
%                 end
%                 if (Fai3{k}(i,j)<epsilon)|(Fai3{k}(i,j)>M)
%                         Fai3{k}(i,j)=Fai3{1}(i,j);
%                 end
%                 if (Fai4{k}(i,j)<epsilon)|(Fai4{k}(i,j)>M)
%                         Fai4{k}(i,j)=Fai4{1}(i,j);
%                 end
%             end
%         end
%     end
%     %计算输出
%     u{k}=u{k-1}+rou*Fai2{k}*(yd{k}(:)-y{k-1}(:)-Fai3{k}*du{2}-Fai1{k}*dy-Fai4{k}*du{3})/(lamda+norm(Fai1{k},2)^2);        
% end
for k=3:N
    %improved projection algorithm proposed by Hou (1999)
    %zerou=zeros(Mu,1);zeroy=zeros(My,1);
    if k==2
        dH=u{k-1}(:)-0;
    else
        dH=u{k-1}(:)-u{k-2}(:);
    end
    for i=2:Lu
        if k>i+1
            dH=[dH;u{k-i}(:)-u{k-i-1}(:)];
        else
            dH=[dH;zerou];
        end
    end
    for i=1:Ly
        if k>i+1
            dH=[dH;y{k-i}(:)-y{k-i-1}(:)];
        else
            dH=[dH;zeroy];
        end
    end
    Fai{k}=Fai{k-1}+eta*(y{k}(:)-y{k-1}(:)-Fai{k-1}(:,:)*dH)*dH'/(miu+norm(dH,2)^2);
%%%%%%%%%% Reset algorithm %%%%%%%%%%%%%%%%
    Fai1{k}=Fai{k}(1:My,1:Mu);
    if (Fai1{k}(1,1)<epsilon)|(Fai1{k}(1,1)>M)
        Fai1{k}(1,1)=Fai1{3}(1,1);
    end
    if (Fai1{k}(2,2)<epsilon)|(Fai1{k}(2,2)>M)
        Fai1{k}(2,2)=Fai1{3}(2,2);
    end
    Fai{k}(1:My,1:Mu)=Fai1{k};
%   control law
    if k==2
        tmp=u{k-1}(:)-0;
    else
        tmp=u{k-1}(:)-u{k-2}(:);
    end
    for i=2:Lu-1
        if k>i+1
            tmp=[tmp;u{k-i}(:)-u{k-i-1}(:)];
        else
            tmp=[tmp;zerou];
        end    
    end
    for i=1:Ly
        if k>i+1
            tmp=[tmp;y{k-i}(:)-y{k-i-1}(:)];
        else
            tmp=[tmp;zeroy];
        end
    end

    u{k}=u{k-1}+rou*Fai1{k}(:,:)'*(yd{k+1}(:)-y{k}(:)-Fai{k}(:,Mu+1:My*Ly+Mu*Lu)*tmp)/(lamda+norm(Fai1{k},2)^2);
    
    
end


