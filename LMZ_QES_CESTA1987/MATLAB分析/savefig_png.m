function savefig_png(figHandle, fileName)
% Save figure without requiring exportgraphics.
try
    set(figHandle,'PaperPositionMode','auto');
    print(figHandle,fileName,'-dpng','-r180');
catch ME
    warning('Could not save figure %s: %s',fileName,ME.message);
end
end
