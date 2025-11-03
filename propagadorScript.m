clearvars;
close all force;
clear;
clc;
clf;

disp("Programa iniciado")

l = 5;

lam = 532E-9;
w0 = 0.2E-3;
k = 2*pi/lam;

N = 2^12;
NV = -N/2:N/2-1;
L = 20*w0;

dx = 2*L/N;
xs = NV*dx;
ys = NV*dx;
[Xs, Ys] = meshgrid(xs, ys);
Rs = sqrt(Xs.^2 + Ys.^2);
theta = atan2(Ys, Xs);

kmax = pi/dx;
kxs = kmax*(2/N)*NV;
kys = kmax*(2/N)*NV;
[KXs, KYs] = meshgrid(kxs, kys);
KTs = sqrt(KXs.^2 + KYs.^2);

zR = k*w0*w0/2;
z = zR / 4;
nz = 25;
dz = z/(nz-1);

disp("Iniciando calculos...")

% Propagacion paraxial
prop = exp(1i*(-0.5*(KTs.^2)/k)*dz);

propagador = Propagador(prop, N, nz);

bloques = {};

%Fuente Gaussiana
disp("Fuente Gaussiana..."); tic
U = exp( - (Rs.^2) / w0^2);

Uz1 = propagador.propagar(U);
bloques{end+1} = Uz1;  

U = Uz1(:, :, end);
toc

% Paso por SPP
disp("Paso por SPP..."); tic
phase = exp(1i * l * theta);
U = U .* phase;

Uz2 = propagador.propagar(U);
bloques{end+1} = Uz2;

U = Uz2(:, :, end);
toc

% Definir rejilla M
disp("Generando rejilla..."); tic
%d = 1E-5;
%d = 2E-5;
d = 4E-5;
%d = 8E-5;
gamma = 2*pi/d;
M1 = exp(1i * rem(gamma * (Xs + L/2), 2*pi));

M_SPP = exp(1i * 1 * theta);
%M_SPP = exp(1i * l * theta);

M = M_SPP .* M1;

M = 0.5 + 0.5*sign(cos(angle(M + pi/4)));
toc

% Paso por rejilla
disp("Paso por rejilla..."); tic
U = U .* M;

Uz3 = propagador.propagar(U);
bloques{end+1} = Uz3;

U = Uz3(:, :, end);
toc

% Paso por lente delgada
disp("Paso por lente delgada..."); tic
f = zR/4;
z = f;
thin_lens = exp(-1i*k/(2*f)*Rs.^2);
U = U.*thin_lens;

Uz4 = propagador.propagar(U);
bloques{end+1} = Uz4;

U = Uz4(:, :, end);
toc

disp("Calculos listos")

disp("Iniciando graficacion...")

%Imax = max(abs(bloques).^2, [], 'all');

video_name = sprintf('l = %d.mp4', l);
v = VideoWriter(video_name,'MPEG-4');
v.FrameRate = max(1, floor(nz/nz));
v.Quality = 100;
open(v);

fig = figure('Color','k','Position',[100 100 1280 720]);
tiledlayout(fig,1,1);
ax1 = nexttile;

I0 = abs(bloques{1}(:,:,1).^2);
P0 = angle(bloques{1}(:,:,1));

hI = imagesc(ax1, xs, xs, I0);
axis(ax1,'image');
colormap(ax1,'hot');
colorbar(ax1);
xlabel(ax1,'x (mm)');
ylabel(ax1,'y (mm)');
title(ax1,'Intensidad');
%clim(ax1,[0 Imax]);

drawnow;
writeVideo(v, getframe(fig));

total_nz = numel(bloques)*nz;

for i = 1:4
    for j = 1:25
        set(hI, 'CData', abs(bloques{i}(:,:,j)).^2);
        title(ax1, sprintf('Intensidad (frame = %d / %d) | l = %d', (i-1)*25+j, total_nz, l));
        drawnow limitrate
        writeVideo(v, getframe(fig));
    end
end
    
close(v);
disp("Animación lista")