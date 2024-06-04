function [f,h]=stackedTracePlot(traces,dt,o)
arguments 
    traces
    dt
    o.spac = 0;
    o.ncols = 1;
    o.cents = [(1:size(traces,1))',(1:size(traces,1))'];
    o.spikes = [];
    o.colors = rand(size(traces,1),3); 
    o.events = []; 
    o.one_t_range = true;
    o.scalebar = false;
    o.plotdur = 60;
end
if ~iscell(o.events), o.events = {o.events}; end
colset = [0 0 0; 1 0 0; 0 0 1];
if length(o.events)>1, evcc = colset(1:length(o.events),:); else evcc =[0  0 0]; end
if isempty(o.spikes), plotspikes = false; else, plotspikes = true; end
h=[];
nt = size(traces,2);
y = sort(o.cents(:,1),'ascend'); 
fracs = round((0:o.ncols)/o.ncols*length(y));
binix = cell(o.ncols,1);
for k = 1:o.ncols
    binix{k} = find(o.cents(:,1)>=y(fracs(k)+1) & o.cents(:,1)<=y(fracs(k+1)));
end
nc = size(traces,1);
scsz = get(0,'screensize');
ht = min(scsz(4)-90,15*nc+100);
f = figure(pos=[5 round(scsz(4)/2-ht/2) 665 ht]);
nt_toplot = floor(nt/o.plotdur);
if o.one_t_range, nt_toplot = 1; end
for tt = 1:nt_toplot
toplot = o.plotdur*(tt-1)+[0 o.plotdur];
toplot(1) = max(dt,toplot(1));
toplot(2) = min(size(traces,2)*dt,toplot(2));
plotix = round(toplot/dt)+1;
plotix = plotix(1):plotix(2);
plotix(plotix<1) = []; plotix(plotix>size(traces,2)) = [];
tplot = linspace(toplot(1),toplot(2),length(plotix));
for j = 1:o.ncols
    tmpmaxprev = 0;
    curcs = binix{j};
%     [~,pltord] = sort(o.cents(curcs,2),'descend');
%     curcs = curcs(pltord);
    subplot(1,o.ncols,sub2ind([1,o.ncols],j,1));%,'sh',0.01,'mt',0.03,'mb',0.05,'ml',0.01','mr',0.01);
    ycoords = [];
    for kk = 1:length(curcs)
        k = curcs(kk);
        tmpplt = traces(k,plotix);
        o.spac = o.spac+tmpmaxprev;
        if kk>1, o.spac = o.spac - min(tmpplt); end
        ycoords(end+1) = o.spac;
        h(end+1)=plot(tplot,tmpplt+o.spac,'color',o.colors(k,:));
        if plotspikes
            cursps = find(o.spikes(kk,plotix));
            hold on;
            if ~isempty(cursps)
            plot(tplot(cursps),prctile(tmpplt+o.spac,95),'.k','markersize',5);
            end
        end
        hold on;
        tmpmaxprev = max(tmpplt);
        tmppltprev = tmpplt;
    end
    axis tight;
%     yticks([]);
    yticks(ycoords);
    yticklabels(1:length(curcs));
    if ~isempty(o.events)
            for k = 1:length(o.events)
                curevents = o.events{k}(ismember(o.events{k},plotix))-plotix(1)+1;
                hold on;
                for iii = 1:length(curevents)
                    plot(repmat(tplot(curevents(iii)),[1 2]),ylim,'--',color=evcc(k,:));
                end
            end
    end
end
xl=xlim;
yl=ylim;
if o.scalebar, plot((xl(2)-xl(1))/2*[1 1]+xl(1),yl(2)-[1 6],'g','linewidth',3); end
xlabel('Time (s)')
drawnow;
if ~o.one_t_range
s = input('next/[q]uit? ','s');
if strcmp(s,'q'), break; end
clf;
end
h=[];
end