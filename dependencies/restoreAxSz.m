function restoreAxSz(despos)
    newpos = get(gca,"position");
    posdif = despos(3:4)-newpos(3:4);
    fpos = get(gcf,"position"); fpos(3:4)=round(fpos(3:4)+posdif);
    set(gcf,pos=fpos);
    set(gca,pos=[newpos(1:2) despos(3:4)]);
end