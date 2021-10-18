


function pn = get_pn()
  syms x;
  
  
  xv = [-2 0 1];
  yv = [0 1 -1];
  n = length(xv)-1;
  pn = 0;
  for k = 0:n
    pn = pn + yv(k+1) * get_Lk(xv,k);
  endfor
  pn = expand(pn);
end




function Lk = get_Lk(xv, k)
  syms x;
  n = length(xv) - 1;
  Lk = 1;
  for j = 0:n
    if j ~= k
      Lk = Lk * (x- xv(j+1)) / (xv(k+1) - xv(j+1));
    end
  endfor
  Lk = expand(Lk);
end  