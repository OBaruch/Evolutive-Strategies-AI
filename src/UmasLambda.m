clear all
close all
clc
%% OMAR BARUCH MORON lOPEZ 
%213605572   
%% Funcion a analizar

% f=@(x,y)((x-2).^2)+((y-2).^2);
% [xfuncion,yfuncion]=meshgrid(-20:.5:20,-10:.5:20);
% z=((xfuncion-2).^2)+((yfuncion-2).^2);
% axis([-20 20 -10 20])

f=@(x,y) x.*exp(-x.^2-y.^2);
[xfuncion,yfuncion]=meshgrid(-20:.1:20,-10:.1:20);
z=xfuncion.*exp(-xfuncion.^2-yfuncion.^2);
axis([-4 4 -4 4])


hold on
%% VARIABLES QUE SE PUEDEN MODIFICAR 
G = 15; %Generaciones
mu = 100; %Poblacion
lambda=30; %Hijos que se tendran
%Este algoritmo reproduce lambda hijos y los mixtea entre los padres y elimina a los peores lambda individuos
%%
xl = [-5; -5]; %Limites 
xu = [5; 5];
D = 2; %Genes


%%
x = zeros(D,mu+lambda);%Poblacion con genes
sigma = zeros(D,mu+lambda); %Sigma de la poblacion

fitness = zeros(1,mu+lambda); %Sirve al final para ordenar y dejar al menos optimo al final
%% Inizializar poblacion
for i=1:mu
    x(:,i) = xl+(xu-xl).*rand(D,1); %Inizlializar corrdenadas random dentro de los limites
    sigma(:,i) = rand(D,1); %Iniziliaziar algun sigma random
end
%% Repetir Generaciones
for t=1:G
    
    for t=mu:mu+lambda
    %%Elejir 2 padres aleatorios difrentes entre si 
    r1 = randi([1 mu]);
    r2 = r1;
    
    while r2==r1
        r2 = randi([1 mu]);
    end
    %
    x(:,t)=(x(:,r1)+x(:,r2))/2; %valores de de dos padres aleatorios entre 2 para nuevos genes del hijo(mu+1) 
    sigma(:,t)=(sigma(:,r1)+sigma(:,r2))/2; 
    
    r = normrnd(0,sigma(:,mu+lambda)); %Valores random con la desviacion estandar para mutar al hijo
    x(:,t) = x(:,t)+r;
    end
    %
    for i=1: mu+lambda
        fitness(i) = f(x(1,i),x(2,i));%Evaluar que tan chingones son 
    end
    
    % ordenar para dejar al final al menos chido que depsues no se tomara
    % en cuenta y terminara por ser sustituido por un hijo mutado todo
    % vergas.
    [~,ind] = sort(fitness);
    
    fitness = fitness(ind);
    x = x(:,ind);
    sigma = sigma(:,ind)
    
    cla;
    for q=1:mu+lambda
        if q>=mu+1
        plot3(x(1,q),x(2,q),f(x(1,q),x(2,q)),'ro');
        end
        plot3(x(1,q),x(2,q),f(x(1,q),x(2,q)),'k*');
    end
     surfc(xfuncion,yfuncion,z);
     pause(.5);
end