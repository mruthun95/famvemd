function TIMF_plot(signal,Colour,nslice,imf,name1,name2)

% default plot attributes
set(groot,'defaultaxesfontname','times');
set(groot,'defaultaxesfontsize',12);
set(groot,'defaulttextInterpreter','latex');
set(groot,'defaultLineLineWidth',2);

[Nx, Ny, Nz] = size(signal);

figure
    xslice = linspace(1,Nx,nslice);
    yslice = linspace(1,Ny,nslice);
    zslice = linspace(1,Nz,nslice);
    volume = slice(signal,xslice,yslice,zslice);
    axis equal;
    xlabel('x');
    ylabel('y');
    zlabel('z');
    switch(name1)
        case 'IMF'
            title(sprintf('%s %d %s',name1,imf,name2));
        case 'Signal'
            title(sprintf('%s %s',name1,name2));
        case 'Residue'
            title(sprintf('%s %s',name1,name2));
    end
    colorbar;
%     caxis([c(1) c(2)]);
    set(volume,'EdgeColor','none',...
        'FaceColor','interp',...
        'FaceAlpha','interp')
    alpha('color')
    view(30,30);
    alphamap('rampup')
    alphamap('decrease',.1)
    colormap(Colour);
