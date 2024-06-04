function figureformat_forsaving

h = gcf;
ht = findall(h,typ="text");
set(ht,'HorizontalAlignment','Center');%;,'VerticalAlignment','middle')
ha = findall(h,typ="axes");
hl = findall(h,typ="Legend");
set(hl,"box","off");
set(ha,'box','off');