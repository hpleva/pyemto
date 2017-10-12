import numpy as np
x=np.linspace(0.0,0.05,6)
y=ens = np.array([
 -0.23396358,-0.23392576,-0.23372901,-0.23334956,-0.23275169,-0.23200480
])
y1 = y - y[0]
y2 = y1 / 2.0
import subprocess
gnuplot = subprocess.Popen(["/usr/bin/gnuplot"], 
                           stdin=subprocess.PIPE)
gnuplot.stdin.write("set term dumb 70 20\n")
gnuplot.stdin.write("plot '-' using 1:2 title 'cprime' with lines, '-' using 1:3 title 'c44' with points\n")
#gnuplot.stdin.write("'-' using 1:3 title 'c44' with points\n")
for i,j,k in zip(x,y1,y2):
   gnuplot.stdin.write("%f %f %f\n" % (i,j,k))
gnuplot.stdin.write("e\n")

for i,j,k in zip(x,y1,y2):
   gnuplot.stdin.write("%f %f %f\n" % (i,j,k))
gnuplot.stdin.write("e\n")

#gnuplot.stdin.write("plot sin(x) title 'Sine Function', 2*sin(x) title 'Tangent'\n")



gnuplot.stdin.flush()
