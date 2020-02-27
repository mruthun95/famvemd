
function IMF_plot(signal,t,imf,name1,name2)    

    plot(t,signal,'-k');
    axis ([0 6*pi -5 5]); % TODO: Line not in EMD1DNV - why?
    xlabel('t');
    set(gca,'TickLabelInterpreter','tex')
    switch(name1)
        case 'IMF'
            title(sprintf('%s %d %s',name1,imf,name2));
        case 'Channel'
            title(sprintf('%s %s',name1,name2));
        case 'Residue'
            title(sprintf('%s %s',name1,name2));
    end
    
    
end