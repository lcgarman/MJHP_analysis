function plot_MJHP_HKL(filename, MINY, MAXY,path_direc)

fileID = fopen(filename,'r');
filestring = split(filename,["-", "_", "."]);
compound = filestring(1);
psp = filestring(2);
psp = strcat("(", psp);
psp = strcat(psp, ")");
reflection = filestring(3);
reflection = strcat("[", reflection);
reflection = strcat(reflection, "]");
ftype = filestring(4);
fig_title = strcat(psp, "  ");
fig_title = strcat(fig_title, compound);
fig_title = strcat(fig_title, "  ");
fig_title = strcat(fig_title, reflection);

%Header
str = fscanf(fileID, '%s',1); %file name
str = fscanf(fileID, '%s',1); %astar
astar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %bstar
bstar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %cstar
cstar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %number of steps
totstep = fscanf(fileID, '%d',1);
str = fscanf(fileID, '%s',1); %starting Energy
Elow = fscanf(fileID, '%f',1);
str = fscanf(fileID, '%s',1); %size of step
Emesh = fscanf(fileID, '%f',1);
str = fscanf(fileID, '%s',1); %HKL
HKL = fscanf(fileID,'%d %d %d',[1 3]);
%End of header

HATOEV = 27.211396;
%Read in reflection Energies
size_rflc = [2 Inf];
rflc = fscanf(fileID,'%d %f',size_rflc);
step = rflc(1,:);
energy = rflc(2,:)*HATOEV;
Yval = (step + Elow)/Emesh;

%initialize varibale for gaussian
nsts = size(energy);
nsts = nsts(2);
%     maxy = max(Yval);
%     miny = min(Yval);
miny=  MINY;
maxy = MAXY;
Eval = miny:0.001:maxy;
nEval = size(Eval);
nEval = nEval(2);
for i = 1:nEval
    RFLC = 0;
    for j = 1:nsts
        sigma = 0.02;%standard deviation; controls width
%             sigma = 0.008; 
        Ei = Eval(i);
        density = energy(j);
        EofSt = Yval(j);
        gauss_rflc = density*exp(-(Ei-EofSt)^2/(sigma));
        RFLC = gauss_rflc + RFLC;
    end
    RFLCplot(i) = RFLC;
end
% figure;
figure('Position', [360,198,300,400]);

maxX = max(RFLCplot);
minX = min(RFLCplot);
abs_minX = abs(minX);
sizeX = size(RFLCplot);
sizeX = sizeX(2);
for j = 1:sizeX
    if RFLCplot(j)>0
        norm_X(j) = RFLCplot(j)/maxX;
    else
        norm_X(j) = RFLCplot(j)/abs_minX;
    end
end

%plot reflection energy values
red = [1,0,0];
blue = [0,0,1];
black = [0,0,0];
length = 100;
arr_blue = [linspace(blue(1),black(1),length)', linspace(blue(2),black(2),length)', linspace(blue(3),black(3),length)'];
arr_red = [linspace(black(1),red(1),length)', linspace(black(2),red(2),length)', linspace(black(3),red(3),length)'];
blue_to_red = cat(1, arr_blue, arr_red);
col = norm_X;
fig = patch([RFLCplot nan],[Eval nan],[col nan],[col nan],'edgecolor', 'interp','linewidth',2.8); 
colormap(blue_to_red);
% colorbar;
hold on;

miny = min(Eval);
maxy = max(Eval);
minx = min(RFLCplot);
maxx = max(RFLCplot);
if abs(minx)>maxx
    maxx = abs(minx);
else 
    minx = -maxx;
end

fermi = 0.0;

% % % making the plot pretty
hold on;
fontname(gcf,"Arial")
axis([minx maxx miny maxy]);
x_range = maxx-minx;
x_step = x_range/4;
set(gca,'TickDir','out');
% xticks(minx:x_step:maxx);
% xticks([minx, 0, maxx]);

xtickformat('%.2f')
% title(['\fontsize{18}' fig_title]);
% xlabel('\fontsize{18} Potential Energy Contribution (eV)');
% ylabel('\fontsize{18} Energy (eV)');
yticks({});
set(gca,'Fontsize',13);
plot([0,0],[miny, maxy],':','linewidth',1,'color',[0,0,0]);
plot([minx,maxx],[fermi,fermi],':','linewidth',2,'color',[0,0,0]);
plot([minx, minx],[miny,maxy],'linewidth',1,'color',[0,0,0]);
plot([minx,maxx],[miny,miny],'linewidth',1,'color',[0,0,0]);
plot([minx,maxx],[maxy, maxy],'linewidth',1,'color',[0,0,0]);
plot([maxx,maxx],[miny,maxy],'linewidth',1,'color',[0,0,0]);
%save figure in directory
if path_direc ~= 0
    fileinput = split(filename,["."]);
    filestore = fileinput(1);
    filestore = strcat(filestore, "-MJHP");
    filestore_tif = strcat(filestore, ".tif");
    filestore_fig = strcat(filestore, ".fig");
    saveas(fig,fullfile(path_direc,[filestore_tif]));
    saveas(fig,fullfile(path_direc,[filestore_fig]));
end

fclose('all');


