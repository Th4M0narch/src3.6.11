function scatter_factor_plot(obs,pred,titleText,xlab,ylab,outFile)
q=isfinite(obs)&isfinite(pred)&obs>0&pred>0;
obs=obs(q); pred=pred(q);
fig=figure('Color','w');
if isempty(obs)
    text(0.1,0.5,'No positive paired data'); axis off;
else
    loglog(obs,pred,'o','MarkerSize',5); hold on;
    lo=min([obs;pred]); hi=max([obs;pred]);
    if lo<=0, lo=min([obs(obs>0);pred(pred>0)]); end
    xx=logspace(log10(lo),log10(hi),100);
    loglog(xx,xx,'k-','LineWidth',1.0);
    loglog(xx,2*xx,'k--',xx,0.5*xx,'k--');
    loglog(xx,5*xx,'k:',xx,0.2*xx,'k:');
    loglog(xx,10*xx,'k:',xx,0.1*xx,'k:');
    xlim([lo hi]); ylim([lo hi]); grid on;
    legend('Model vs experiment','1:1','FAC2 bounds','','FAC5 bounds','','1:10','10:1','Location','best');
end
xlabel(xlab); ylabel(ylab); title(titleText);
savefig_png(fig,outFile);
end
