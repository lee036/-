
load('h.mat');
ne0=load('ne.mat');
load('ni.mat');
load('T0.mat');
Te0=load('Te.mat');
load('nn.mat');
%����ԭʼ����
t=0:1e-1:1000;N=length(t);I=length(h);k=1.38e-23;me=0.91e-30;He=2000;dt=1e-1;Vn=100;
%ȷ���������߶ȷ�Χ
Ke=zeros(I,N);Le=zeros(I,N);ven=zeros(I,N); vei=zeros(I,N);
k1=zeros(I,N); k2=zeros(I,N);ne=zeros(I,N);Te=zeros(I,N);MV=zeros(I,N);
we=zeros(I,N);
% Ԥ���ռ���м���
ne(:,1)=ne0.ne;Te(:,1)=Te0.Te;

f0=9.12e6;H=230; %Ϊ�����߶�

[wpe,fpe] = plamsafrequence(ne(:,1));   %����ԭʼ���ݣ��������ò�����
[MVo,MVo2,MVno] = ionofrequency(nn,ni,T0);  %��������ԭ����ײƵ��

[Q0,S0,Le(:,1),Ke(:,1),ven(:,1),vei(:,1),k1(:,1),k2(:,1),we(:,1)] =originalvalue(ne,Te,T0,ni,nn,I,He);
Le0=Le(:,1);Ke0=Ke(:,1);we0=we(:,1);ven0=ven(:,1);k10=k1(:,1);k20=k2(:,1);
Qrw=heating(Te(:,1),ne(:,1),h,nn,fpe,f0);%������ʼ������

for n=1:N-1
    %������µ����¶�
    
    Te(2:I-1,n+1)=(-3/2*k*ne(2:I-1,n).*we(2:I-1,n).*diff(Te(2:I,n))/He-k.*ne(2:I-1,n).*Te(2:I-1,n)...
        .*diff(we(1:I-1,n))/He+diff(Ke(1:I-1,n).*diff(Te(:,n)))/(He^2)+S0(2:I-1)+Qrw(2:I-1)+Le(2:I-1,n))./(3/2*k*ne(2:I-1,n))*dt+Te(2:I-1,n);  %�����¶ȼ��ȹ���
    Te(I,n+1)=Te(I,1)+Te(I-1,n+1)-Te(I-1,1); Te(1,n+1)=Te(1,1);%�߽�����
    %���²���
    %plot(h(2:I-1),diff(Ke(1:I-1,1).*diff(Te(:,n)))/(1000^2));
    Le(:,n)=energyloss(ne(:,n),ni,Te(:,n+1),T0,nn,I);   %����������ʧ����
    % plot(h,Le(:,n)-Le0);
    Ke(:,n)=thermalconductance(Te(:,n+1),nn,ne(:,n));  %�ȵ��ʺ���
    %  plot(h,Ke(:,n)-Ke0);
    [ven(:,n),vei(:,n)]=collisionfrequency(Te(:,n+1),ne(:,n),nn); %��ײƵ��
    %  plot(h,ven(:,n),h,vei(:,n))
    fh=@(te)(4.2e-13*(300./te).^0.85);
    k1(:,n)=fh(Te(:,n+1));         %NO+����ϵ��
    % plot(h,k1(:,n)-k10);
    fh1=@(x)(1.6e-13*(300./x).^0.55);
    k2(:,n)=fh1(Te(:,n+1));
    %  plot(h,k2(:,n)-k20);%O2+����ϵ��
    MV(:,n)=ne(:,n)./(MVo+MVo2+MVno);
    %������µ���ܶ�
    we(1:I-1,n)=-1./ne(2:I,n)*1./(me.*ven(2:I,n)+MV(2:I,n)).*diff(ne(:,n)*k.*(Te(:,n+1)+T0))/He+Vn;
    ne(2:I-1,n+1)=-diff(ne(2:I,n).*we(1:I-1,n))/He*dt+Q0(2:I-1)*dt-(k1(2:I-1,n).*ni(2:I-1,3)+k2(2:I-1,n).*ni(2:I-1,2)).*ne(2:I-1,n)*dt+ne(2:I-1,n);
    ne(I,n+1)=ne(I,1)+ne(I-1,n+1)-ne(I-1,1); ne(1,n+1)=ne(1,1);%�߽�����
    %  plot(h,ne(:,n+1));
    %   pause(1);
    %���¼������
    [ven(:,n),vei(:,n)]=collisionfrequency(Te(:,n+1),ne(:,n+1),nn);
    we(1:I-1,n+1)=-1./ne(2:I,n+1)*1./(me.*ven(2:I,n)+MV(2:I,n)).*diff(ne(:,n+1)*k.*(Te(:,n+1)+T0))/He+Vn; %���µ���Ǩ���ٶ�
    
    %   plot(h(1:I-1),we(1:I-1,n+1));
    %   pause(1);
    Le(:,n+1)=energyloss(ne(:,n+1),ni,Te(:,n+1),T0,nn,I);   %����������ʧ����
    %  plot(h,Le(:,n+1)-Le0);
    Ke(:,n+1)=thermalconductance(Te(:,n+1),nn,ne(:,n+1));  %�ȵ��ʺ���
    %  plot(h,Ke(:,n+1)-Ke0);
    [ven(:,n+1),vei(:,n+1)]=collisionfrequency(Te(:,n+1),ne(:,n+1),nn); %��ײƵ��
    % plot(h,ven(:,n+1),h,vei(:,n+1))
end



