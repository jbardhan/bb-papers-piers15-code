qburied=withstatic_asym_buried-withstatic_sym_buried;
qsurface=withstatic_asym_surface-withstatic_sym_surface;
figure;
set(gca,'fontsize',16);
plot(radius, qsurface(:,1), 'b-o', 'linewidth',2,'markersize',10);
hold on;
plot(radius, qburied(:,1),'b-s','linewidth',2,'markersize',10);
plot(radius, qsurface(:,2),'r-o','linewidth',2,'markersize',10);
plot(radius, qburied(:,2),'r-s','linewidth',2,'markersize',10);
xlabel('Sphere radius (Angstroms)');
ylabel('Energy difference (kcal/mol)')
legend('Surface, q = -1', 'Buried, q = -1',...
       'Surface, q = +1','Buried, q = +1','location','southwest');
