syms n1 r w theta x y delta n2 B

v = 5*(1-n1)*r*w;
f1 = v*cos(theta)
f2 = v*sin(theta)
f3 = v/B*tan(delta+n2)
diff(f3,n1)
diff(f3,n2)
