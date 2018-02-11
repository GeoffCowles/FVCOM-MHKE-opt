close all; clear all;
optparams;
load mhke;
rectangle('Position',[-0.5*chan_length,-0.5*chan_width,chan_length,chan_width]); drawnow; hold on;
marks = {'r+','b+','g+','k+','m+','rd','bd','gd','kd','md'};

P = [mhke.tpower];
msize = 1 + 10*(P-min(P))/(max(P)-min(P));

nmhke = numel(P);
for i=1:nmhke;
  plot(mhke(i).x,mhke(i).y,char(marks(i)),'MarkerSize',ceil(msize(i)))
end;