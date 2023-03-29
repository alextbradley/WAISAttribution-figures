fname = 'outfile.nc';

t = ncread(fname, 'TIME');
x = ncread(fname, 'x');
h = ncread(fname, 'h');
b = ncread(fname, 'b');
gr = ncread(fname, 'grfrac');

for i = 1:3:length(t)
    clf; hold on; title(['t = ' num2str(t(i))])

    bb = squeeze(b(:,25,i));
    hh = squeeze(h(:,25,i));

    plot(x,bb, 'k', 'linewidth', 2)
    
    bottom = -918/1028*hh;
    bottom(bottom < bb) = bb(bottom < bb);
    plot(x,bottom, 'r', 'linewidth', 2)

    drawnow 
    pause
end