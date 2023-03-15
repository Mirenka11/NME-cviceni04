% zvolme pocet iteraci
iteraci=100;

% zvolime nahodnou matici, uvidime jeslti bude konvergovat, ke ktere 
% pricteme jednotkovou diagonalni, zvysime tak pravdepodobnost ze vse bude konvergovat
A=rand(3)/2+eye(3)/2

% zvolime pozadovane reseni, ze ktereho napocitame pravou stranu
origReseni=[1;1;1];

% najdeme pravou stranu
y=A*origReseni

% zvolime nahodny uvodni odhad ze ktereho budeme iterovat
uvodniOdhad=rand(3,1);

% zvolime omega pro superrelaxaci, nebude optimalni, nevadi
omega=[0.1, 0.3, 0.7, 1, 1.5, 1.9];

%% GAUSS-SEIDEL

	L=tril(A,-1);
	R=triu(A,1);
	D=A-L-R;

	F=-inv(D+L)*R;
	G=inv(D+L);

	x=uvodniOdhad;
	for i=1:iteraci
		x=F*x+G*y;
		vzdGaussSeidel(i)=norm(x-origReseni);
	end;

    osaX=1:iteraci;

	plot(osaX,vzdGaussSeidel, '--','LineWidth',3);
	% zobrazime legendu
	legend('Gauss-Seidel');
	% nastavime osu Y jako logaritmickou
	set(gca,'YScale','log');
    hold on

%% volim ruzne omega
    for i=1:6
        vzdSuperRelaxace = superrelax(A,omega(i), uvodniOdhad, iteraci, y, origReseni);
        txt = ['omega=',num2str(omega(i))];
        plot(osaX,vzdSuperRelaxace(1,:),'DisplayName',txt);   
    end

    hold off
    legend show
%%
function[vzdSuperRelaxace] = superrelax(A, omega, uvodniOdhad, iteraci,y, origReseni)
    L=tril(A,-1);
	R=triu(A,1);
	D=A-L-R;

	% toto je ze seidela
	F=-inv(D+L)*R;
	G=inv(D+L);

	x=uvodniOdhad;
	for i=1:iteraci
		xkp1=F*x+G*y;
		x=omega*(xkp1-x)+x;
		vzdSuperRelaxace(i)=norm(x-origReseni);
	end;
    
end
