set term pdfcairo #noenhanced
set output "cmpAccuracy.pdf"

set title "accuracy"
set key outside below maxrow 5
set logscale y 10
set format y '10^{%L}'

set ylabel 'rel max err'
sdcCol="magenta black blue - - red - green"
ros2File="ros2/accuracy.dat"

rkcFile="rkc/stages_8/accuracy.dat"
rkcName="rkc_8"
rkcColor="orange"
rkcDash=2


set logscale x 10
set yrange [*:1]
set format x '10^{%L}'
set xrange [*:200]
set xlabel 'n step'

mrFo='mrsdc/kiter_%d/accuracy_M_%s_P_%s.dat'
mrLapl0Fo='mrsdc_lapl0/kiter_%d/accuracy_M_%s_P_%s.dat'
mrLaplTildeHFo='mrsdc_laplTilde_h/kiter_%d/accuracy_M_%s_P_%s.dat'
mrLaplTildeSqrtHFo='mrsdc_laplTilde_sqrtH/kiter_%d/accuracy_M_%s_P_%s.dat'
mrLaplTildeH2Fo='mrsdc_laplTilde_h2/kiter_%d/accuracy_M_%s_P_%s.dat'

#sdcSol='sdc/kiter_%d/accuracy_M_%s.dat'
#sdcName='sdc-M%s-kiter%d'

#sdcSol='sdc_lapl0/kiter_%d/accuracy_M_%s.dat'
#sdcName='sdc-Lapl_{0}-M%s-kiter%d'

sdcSol='sdc_sqrtH/kiter_%d/accuracy_M_%s.dat'
sdcName='sdc-Lapl-h^{1/2}'

#sdcSol='sdc_1/kiter_%d/accuracy_M_%s.dat'
#sdcName='sdc-Lapl-1'

do for [k_iterC in "1 2 3 4"] {
k_iter=k_iterC+0 #convert to number
print sprintf("kiter: %d\n", k_iter)


plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrFo,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf(sdcName, M,k_iter), rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLapl0Fo,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl_{0}', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf(sdcName, M,k_iter), rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLaplTildeHFo,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl\~h', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf(sdcName, M,k_iter), rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLaplTildeSqrtHFo,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl\~h^{1/2}', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf(sdcName, M,k_iter), rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLaplTildeH2Fo,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl\~h^2', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf(sdcName, M,k_iter), rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName
}

set logscale x 10
set format x '10^{%L}'
set xlabel 'runtime'
set xrange [5:200]
do for [k_iterC in "1 2 3 4"] {
k_iter=k_iterC+0 #convert to number
print sprintf("kiter: %d\n", k_iter)

#mrFo='mrsdc/kiter_%d/accuracy_M_%s_P_%s.dat'
#mrLapl0Fo='mrsdc_lapl0/kiter_%d/accuracy_M_%s_P_%s.dat'
#mrLaplTildeSqrtHFo='mrsdc_laplTilde_sqrtH/kiter_%d/accuracy_M_%s_P_%s.dat'
#mrLaplTildeHFo='mrsdc_laplTilde_h/kiter_%d/accuracy_M_%s_P_%s.dat'
#mrLaplTildeH2Fo='mrsdc_laplTilde_h2/kiter_%d/accuracy_M_%s_P_%s.dat'
#sdcFo='sdc/kiter_%d/accuracy_M_%s.dat'
#sdcLapl0Fo='sdc_lapl0/kiter_%d/accuracy_M_%s.dat'

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrFo,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P, k_iter), "ros2/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf(sdcName, M, k_iter), rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLapl0Fo,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl_{0}', M,P, k_iter), "ros2/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf(sdcName, M, k_iter), rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLaplTildeHFo,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl\~h', M,P, k_iter), "ros2/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf(sdcName, M, k_iter), rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLaplTildeSqrtHFo,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl\~h^{1/2}', M,P, k_iter), "ros2/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf(sdcName, M, k_iter), rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLaplTildeH2Fo,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl\~h^2', M,P, k_iter), "ros2/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcSol,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf(sdcName, M, k_iter), rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

}

#set title "runtime per step"
#
#set yrange [*:*]
#set ylabel "runtime"
#set xlabel "n step"
#set logscale y 2
#set format y '2^{%L}'
#do for [k_iterC in "1 2 3 4"] {
#
#k_iter=k_iterC+0 #convert to number
#plot for [M in "2 3"] for [P in "2 6"] sprintf(mrFo,k_iter,M,P) u 4:5 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P, k_iter), for [M in "2 3"] sprintf(sdcFo,k_iter,M) u 4:5 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter),
#
#}
