function [] = plot_grow(x,R,I,A,M,S,fig_loc)

  axes_loc = axes('Parent',fig_loc);
  plot_loc = plot(x*R,I,'--',x*R,A,x*R,M,':',x*R,S,':');
  
  set(plot_loc(1),'DisplayName','Inhibitor');
  set(plot_loc(2),'DisplayName','Activator');
  set(plot_loc(3),'DisplayName','Maturation');
  set(plot_loc(4),'DisplayName','Mesenchyme');

  set(plot_loc(1),'LineWidth',3,'Color','blue');
  set(plot_loc(2),'LineWidth',3,'Color','red');
  set(plot_loc(3),'LineWidth',3,'Color','green');
  set(plot_loc(4),'LineWidth',3,'Color','yellow');

  ylim([-2 2])
  yticks([])
  yticklabels({})
  xlim([0 R])
  xticks([0 R])
  xticklabels({'0','R(t)'})
  
  set(axes_loc, 'YLimMode', 'manual','XLimMode', 'manual' )
  set(axes_loc,'FontSize',40);
  box(axes_loc,'on');
  legend(axes_loc,'show','location','south','orientation','horizontal')
  hold off
  drawnow
