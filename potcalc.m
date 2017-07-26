directory = './out/bin_ad_test/';
mesh = createMesh(directory);
field = createField(sprintf('%sgasdens0.dat', directory), mesh);
ab = 0.2288;
eb = 0.52;
per(1) = 0.0;
per(2) = pi;
aspect_ratio = 0.05;
flaring_index = 0.143;
G = 1;

m(2) = 1E-32;
m(1) = 1-m(2);
q(1) = m(1)/m(2);
q(2) = m(2)/m(1);
x(1) = (1.0-m(1))*ab*cos(per(1))*(1.0-eb);
x(2) = (1.0-m(2))*ab*cos(per(2))*(1.0-eb);
y(1) = 0.0;
y(2) = 0.0;
bsmooth(1) = 0.4;
bsmooth(2) = 0.4;

drx = x(1) - x(2);
dry = y(1) - y(2);
A = sqrt(drx*drx + dry*dry);
for s = 1:2
  rs = sqrt(x(s)*x(s) + y(s)*y(s));
  
  Rstar(s) = (0.49*power(q(s), 2.0/3.0))/(0.6*power(q(s), 2.0/3.0) + log(1.0 + power(q(s), 1.0/3.0)))*A;
  smoothing = bsmooth(s) * aspect_ratio * power(rs, 1.0+flaring_index);
  ssmooth[s] = smoothing*smoothing;
end

Pot_rr = zeros(mesh.nx, mesh.ny);
Pot_ts = zeros(mesh.nx, mesh.ny);
Pot_s = zeros(mesh.nx, mesh.ny);

for s = 1:2
  rs = sqrt(x(s)*x(s) + y(s)*y(s));
  for i = 1:mesh.ny
    InvD = 1.0/mesh.ymed(i);
    for j = 1:mesh.nx;
      ang = (j/mesh.nx)*2.0*PI;
      x = mesh.ymed(i)*cos(ang);
      y = mesh.ymed(i)*sin(ang);
      d = (x-x(i))*(x-x(i)) + (y-y(i))*(y-y(i));
      if sqrt(d) < Rstar(s)
        pot_rr = (-1.0)*G*m(s)*(3.0*power(Rstar(s), 2.0) - d)/(2.0*power(Rstar(s), 3.0));
      else
        pot_rr = (-1.0)*G*m(s)/sqrt(d);
      end
      dsmooth = sqrt(d + ssmooth(s));
      pot_ts = (-1.0)*G*m(s)/dsmooth;

      Pot_rr(j,i) = Pot_rr(j,i) + pot_rr;
      Pot_ts(j,i) = Pot_ts(j,i) + pot_ts;

      pot = (-1.0)*G*1.0*InvD;
     	Pot_s(j,i) = pot;
    end
  end
end

Pot_spheres = field; Pot_spheres.data = Pot_rr;
Pot_smooth = field; Pot_smooth.data = Pot_ts;
Pot_single = field; Pot_single.data = Pot_s;


figure(1); clf;
plotField(Pot_spheres);

figure(2); clf;
plotField(Pot_smooth);

figure(3); clf;
plotField(Pot_single);
