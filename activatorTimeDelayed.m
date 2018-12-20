function [vd] = activatorTimeDelayed(x,R,nt,ntau,v)
  nx = size(v,1);
  vd = zeros(nx,1);
  if nt > ntau
    x_delay       = R(nt-ntau)/R(nt);
    index_x_delay = find(x<=x_delay,1,'last');
    for ii=1:index_x_delay
      new_index   = find(x<=x(ii)/x_delay,1,'last');
      vd(ii) = v(new_index,nt-ntau);
    end
  end
end 
