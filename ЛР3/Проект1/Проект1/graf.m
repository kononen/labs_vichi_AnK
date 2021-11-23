Pn = fopen ('polynom.txt','rt');
if Pn == -1 
    error('File is not opened'); 
end 

x=0;                % ������������� ���������� 
fpn=0;
f=0;
cnt=1;              % ������������� �������� 
while ~feof(Pn)    % ����, ���� �� ��������� ����� ����� 
    [Vx,Nx] = fscanf(Pn, '%f',1);  %���������� ������ 
    [Vf,Nf] = fscanf(Pn, '%f',1);
    [Vf_true,Nf_ftrue] = fscanf(Pn, '%f',1);
% �������� double (V �������� �������� 
% ��������, N � ����� ��������� ���������) 
    if (Nf_ftrue > 0)        % ���� ������� ��� �������� �������, �� 
            x(cnt)=Vx;   % ��������� ������-������ �� �������� V
            fpn(cnt)=Vf;
            f(cnt)=Vf_true;
            cnt=cnt+1;  % ����������� ������� �� 1     
    end 
end 
fclose(Pn); 

plot(x,fpn,'ro-', x, f,'b')% ������� �������, ����� �������� �������
 legend('�����������','�������� �������');