set term pdfcairo #noenhanced
set output "cmpAccuracy_sdc_ros2.pdf"

k_iter=2
set key outside below maxrow 5
set logscale y 2
set logscale x 2
set format y '10^{%T}'
set format x '10^{%T}'

set ylabel 'max err'
sdcCol="- black blue - - red - green"

set xlabel 'n step'
plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf('sdc/kiter_%d/accuracy_M_%s_P_%s.dat',k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('sdc-M%s-P%s-kiter%d', M,P,2), "ros2/accuracy.dat" u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2"

set logscale x 10
set xlabel 'runtime'
plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf('sdc/kiter_%d/accuracy_M_%s_P_%s.dat',k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('sdc-M%s-P%s-kiter%d', M,P,2), "ros2/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2"
