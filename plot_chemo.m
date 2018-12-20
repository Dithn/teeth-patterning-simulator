function [] = plot_chemo(x,L,I,A,C,fig_loc)

  axes_loc = axes('Parent',fig_loc);
  plot_loc = plot(x,I,'--',x,A,x,C,':');
  
  set(plot_loc(1),'DisplayName','Inhibitor');
  set(plot_loc(2),'DisplayName','Activator');
  set(plot_loc(3),'DisplayName','Cells density');

  set(plot_loc(1),'LineWidth',3,'Color','blue');
  set(plot_loc(2),'LineWidth',3,'Color','red');
  set(plot_loc(3),'LineWidth',3,'Color','green');

  ylim([-1.5 3])
  yticks([])
  yticklabels({})
  xlim([0 L])
  xticks([0 L])
  
  set(axes_loc, 'YLimMode', 'manual','XLimMode', 'manual' )
  set(axes_loc,'FontSize',40);
  box(axes_loc,'on');
  legend(axes_loc,'show','location','south','orientation','horizontal')
  hold off
  drawnow
