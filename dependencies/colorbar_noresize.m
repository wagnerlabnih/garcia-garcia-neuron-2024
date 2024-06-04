function h = colorbar_noresize

axh = findall(gcf,typ="axes");
nax = length(axh);
origpos = {};
for i = 1:nax
    origpos{i} = axh(i).Position;
end
h = colorbar;
for i = 1:nax
    axh(i).Position = origpos{i};
end
