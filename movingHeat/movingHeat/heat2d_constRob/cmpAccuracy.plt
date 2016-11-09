set term pdfcairo #noenhanced
set output "cmpAccuracy.pdf"

set title "accuracy"
set key outside below maxrow 5
set logscale y 10
set format y '10^{%L}'

set ylabel 'rel max err'
sdcCol="magenta black blue - - red - green"
ros2File="ros2/accuracy.dat"


set logscale x 10
set format x '10^{%L}'
set xrange [4:200]
set yrange [*:1]
set xlabel 'n step'
do for [k_iterC in "1 2 3 4"] {
k_iter=k_iterC+0 #convert to number
print sprintf("kiter: %d\n", k_iter)

mrFoSrcLapl0='mrsdc/MrY_asSlow_withSlowSrc_lapl0/kiter_%d/accuracy_M_%s_P_%s.dat'
mrFoSrc='mrsdc/MrY_asSlow_withSlowSrc/kiter_%d/accuracy_M_%s_P_%s.dat'
mrFo='mrsdc/MrY_asSlow_withoutSlowSrc/kiter_%d/accuracy_M_%s_P_%s.dat'
sdcFo='sdc/kiter_%d/accuracy_M_%s.dat'

plot for [M in "2 3"] for [P in "2 8"] sprintf(mrFoSrc,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcSrc-M%s-P%s-kiter%d', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter)

plot for [M in "2 3"] for [P in "2 8"] sprintf(mrFoSrcLapl0,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcSrcLapl0-M%s-P%s-kiter%d', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter)

plot for [M in "2 3"] for [P in "2 8"] sprintf(mrFo,k_iter,M,P) u (20.0/$1):3 w lp dt 3 lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P,k_iter), for [M in "2 3"] for [P in "2 8"] sprintf(mrFoSrc,k_iter,M,P) u (20.0/$1):3 w lp dt 2 lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcSrc-M%s-P%s-kiter%d', M,P,k_iter)

}

#set logscale x 10
#set format x '10^{%L}'
#set xlabel 'runtime'
#set xrange [*:*]
#do for [k_iterC in "1 2 3 4"] {
#k_iter=k_iterC+0 #convert to number
#print sprintf("kiter: %d\n", k_iter)
#
#mrFo='mrsdc_alpha_1e-3/kiter_%d/accuracy_M_%s_P_%s.dat'
#mrLapl0Fo='mrsdc_alpha_1e-3_lapl0/kiter_%d/accuracy_M_%s_P_%s.dat'
#sdcFo='sdc_alpha_1e-3/kiter_%d/accuracy_M_%s.dat'
#sdcLapl0Fo='sdc_alpha_1e-3_lapl0/kiter_%d/accuracy_M_%s.dat'
#
#plot for [M in "2 3"] for [P in "2 3"] sprintf('mrsdc_alpha_1e-3/kiter_%d/accuracy_M_%s_P_%s.dat',k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P, k_iter), "ros2_alpha_1e-3/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter),
#
#plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLapl0Fo,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl_{0}', M,P, k_iter), "ros2_alpha_1e-3/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcLapl0Fo,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d-Lapl_{0}', M, k_iter),
#
#}

#set title "runtime per step"
#
#set ylabel "runtime"
#set xlabel "n step"
#set logscale y 2
#set format y '2^{%L}'
#do for [k_iterC in "1 2 3 4"] {
#
#k_iter=k_iterC+0 #convert to number
#plot for [M in "2 3"] for [P in "2 6"] sprintf('mrsdc_alpha_1e-3/kiter_%d/accuracy_M_%s_P_%s.dat',k_iter,M,P) u 4:5 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P, k_iter), for [M in "2 3"] sprintf('sdc_alpha_1e-3/kiter_%d/accuracy_M_%s.dat',k_iter,M) u 4:5 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter),
#
#}
