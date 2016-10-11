set term pdfcairo #noenhanced
set output "cmpAccuracy_sdc_ros2.pdf"

do for [k_iterC in "1 2"] {
k_iter=k_iterC+0 #convert to number
print sprintf("kiter: %d\n", k_iter)
set key outside below maxrow 5
set logscale y 10
set logscale x 2
set format y '10^{%L}'
set format x '2^{%L}'

set ylabel 'max err'
sdcCol="magenta black blue - - red - green"

set xlabel 'n step'
plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf('mrsdc/kiter_%d/accuracy_M_%s_P_%s.dat',k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P,k_iter), "ros2/accuracy.dat" u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf('sdc/kiter_%d/accuracy_M_%s.dat',k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter)

set logscale x 10
set format x '10^{%L}'
set xlabel 'runtime'
plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf('mrsdc/kiter_%d/accuracy_M_%s_P_%s.dat',k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P, k_iter), "ros2/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf('sdc/kiter_%d/accuracy_M_%s.dat',k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter),
}
