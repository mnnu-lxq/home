clear all;close all;clc;
G = gsp path(5);
W = [] %Enter weight matrix here
G.W = W;
G = gsp graph default parameters(G);
G.plotting.vertex size = 200;
G.coords = [0,0;2,2;3,-1;5,1;4,3;]
G.plotting.limits = [-1,6,-2,4];
%-- Define a Graph Signal here
f = [2 -1 2 -1 1]';
%---Graph Spectrum
Gf = gsp compute fourier basis(G);
paramplot.colorbar = 1;
%plotting the signal (sticks)
figure;%gplot(G.A,G.coords,'-or');hold on;box off;axis off;%title('ND signal on ER Graph');
gsp plot graph(G);hold on;box off;axis off;title('Graph signal');
for m=1:G.N
plot3([G.coords(m,1),G.coords(m,1)],[G.coords(m,2),G.coords(m,2)],[0,f(m)],'b','linewidth',3);
end
view([25,45]);
%Plotting spectrum of the graph signal
f hat = gsp gft(Gf,f);
figure;
gsp plot signal spectral(Gf,f hat);
title('Spectrum of the graph signal');
set(gca,'fontsize',20);
xlabel('$\lambda \ell$','Interpreter','LaTex','FontSize',20);
ylabel('$\hat{f}(\lambda \ell) $','Interpreter','LaTex','FontSize',20);
%----Filter definition
tau = 5;
h = @(x) 1./(1+tau*x);
figure;
gsp plot filter(Gf,h);
title('Filter h');
% Perform the filtering operation
f2 = gsp filter(Gf,h,f);
figure;%gplot(G.A,G.coords,'-or');hold on;box off;axis off;title('ND signal on ER Graph');
gsp plot graph(G);hold on;box off;axis off;title('Graph signal (sticks)');
for m=1:G.N
plot3([G.coords(m,1),G.coords(m,1)],[G.coords(m,2),G.coords(m,2)],[0,f2(m)],'b','linewidth',3);
end
title('Filtered Signal');
view([25,45]);
