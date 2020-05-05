% Plotting a Graph Signal and its Spectrum

clear all;close all;clc;
%-- Enter weight matrix of the graph
W =

G.W = W;
G = gsp_graph_default_parameters(G);

G.plotting.vertex_size = 100;
G.coords = [0,0;2,2;3,-1;5,1;4,3;]
G.plotting.limits = [-1,6,-2,4];
%-- Define a Graph Signal here
f =

%---Graph Spectrum
Gf = gsp_compute_fourier_basis(G);
paramplot.colorbar = 1;
%plotting the signal (color coded)
figure;
gsp_plot_signal(G,f,paramplot);
set(gca,'fontsize',20)
title('Graph signal (color coded)');

%plotting the signal (sticks)
figure;
gsp_plot_graph(G); hold on; box off; axis off;
title('Graph signal (sticks)');
for m=1:G.N
  plot3([G.coords(m,1),G.coords(m,1)], ...
    [G.coords(m,2),G.coords(m,2)], ...
    [0,f(m)],'b','linewidth',3);
end
view([25,45]);

%Plotting spectrum of the graph signal
f_hat = gsp_gft(Gf,f);
figure;
gsp_plot_signal_spectral(Gf,f_hat);
title('Spectrum of the graph signal');
set(gca,'fontsize',20);
xlabel('$\lambda_\ell$','Interpreter','LaTex','FontSize',20);
ylabel('$\hat{f}(\lambda_\ell) $','Interpreter','LaTex','FontSize',20);
