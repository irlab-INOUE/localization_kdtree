set pm3d 
set pm3d map
set xrange [-1:1]
set yrange [-1:1]
#set ticslevel 0
#set cbrange[0:1]
set xlabel 'X[m]'
set ylabel 'Y[m]'
#set palette defined (0 "blue", 0.5 "white", 1 "red")
splot "log" with pm3d

