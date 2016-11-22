set term pdfcairo #noenhanced
set output "cmpAccuracy.pdf"

set title "accuracy"
set key outside below maxrow 5
set logscale y 10
set format y '10^{%L}'

set ylabel 'rel max err'
sdcCol="magenta black - blue - red - green"
ros2File="ros2/accuracy.dat"

rkcFile="rkc/stages_4/accuracy.dat"
rkcName="rkc_4"
rkcColor="orange"
rkcDash=2

mrFoSrcLapl0='mrsdc_lapl0/kiter_%d/accuracy_M_%s_P_%s.dat'
mrFoLaplh='mrsdc_laplTilde_h/kiter_%d/accuracy_M_%s_P_%s.dat'

mrFoLaplSqrth='mrsdc_laplTilde_sqrtH/kiter_%d/accuracy_M_%s_P_%s.dat'
mrLaplSqrthName='mrsdcLapl\~h^{1/2}-M%s-P%s-kiter%d'

mrFoLaplSqrthTh0='mrsdc_laplTilde_sqrtH_theta_0/kiter_%d/accuracy_M_%s_P_%s.dat'
mrLaplSqrthTh0Name='mrsdcLapl\~h^{1/2}{/Symbol q}-M%s-P%s-kiter%d'

mrFoLaplh2='mrsdc_laplTilde_h2/kiter_%d/accuracy_M_%s_P_%s.dat'
mrFoSrc='mrsdc/kiter_%d/accuracy_M_%s_P_%s.dat'
#mrFo='mrsdc/MrY_asSlow_withoutSlowSrc/kiter_%d/accuracy_M_%s_P_%s.dat'
sdcFo='sdc/kiter_%d/accuracy_M_%s.dat'

set logscale x 10
set format x '10^{%L}'
set xrange [4:500]
set yrange [*:1]
set xlabel 'n step'
do for [k_iterC in "1 2 3 4"] {
k_iter=k_iterC+0 #convert to number
print sprintf("kiter: %d\n", k_iter)


plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoSrc,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

#plot for [M in "2 3"] for [P in "2 8"] sprintf(mrFoSrcLapl0,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcLapl0-M%s-P%s-kiter%d', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoLaplh,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcLapl\~h-M%s-P%s-kiter%d', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 4 8"] sprintf(mrFoLaplSqrth,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf(mrLaplSqrthName, M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 4 8"] sprintf(mrFoLaplSqrthTh0,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf(mrLaplSqrthTh0Name, M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoLaplh2,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcLapl\~h^2-M%s-P%s-kiter%d', M,P,k_iter), ros2File u (20.0/$1):3 w l lc rgbcolor 'orange' t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,1) pt (M+0) title sprintf('sdc-M%s-kiter%d', M,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoLaplh2,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcSrcLapl\~h^2-M%s-P%s-kiter%d', M,P,k_iter),for [M in "2 3"] for [P in "2 4"] sprintf(mrFoSrc,k_iter,M,P) u (20.0/$1):3 w l dt 2 lc rgbcolor word(sdcCol,P+0) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoLaplh,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcSrcLapl\~h-M%s-P%s-kiter%d', M,P,k_iter),for [M in "2 3"] for [P in "2 4"] sprintf(mrFoSrc,k_iter,M,P) u (20.0/$1):3 w l dt 2 lc rgbcolor word(sdcCol,P+0) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

plot for [M in "2 3"] for [P in "2 4 8"] sprintf(mrFoLaplSqrth,k_iter,M,P) u (20.0/$1):3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf(mrLaplSqrthName, M,P,k_iter),for [M in "2 3"] for [P in "2 4"] sprintf(mrFoSrc,k_iter,M,P) u (20.0/$1):3 w l dt 2 lc rgbcolor word(sdcCol,P+0) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

#plot for [M in "2 3"] for [P in "2 8"] sprintf(mrFo,k_iter,M,P) u (20.0/$1):3 w lp dt 3 lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P,k_iter), for [M in "2 3"] for [P in "2 8"] sprintf(mrFoSrc,k_iter,M,P) u (20.0/$1):3 w lp dt 2 lc rgbcolor word(sdcCol,P+0) pt (M+0) title sprintf('mrsdcSrc-M%s-P%s-kiter%d', M,P,k_iter) #, rkcFile u 4:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName

}

set logscale x 10
set format x '10^{%L}'
set xlabel 'runtime'
set xrange [*:*]
set xrange [100:*]
do for [k_iterC in "1 2 3 4"] {
k_iter=k_iterC+0 #convert to number
print sprintf("kiter: %d\n", k_iter)

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoSrc,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d', M,P, k_iter), ros2File u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter) #, rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoLaplh,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdcLapl\~h-M%s-P%s-kiter%d', M,P, k_iter), ros2File u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter) #, rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoLaplSqrth,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf(mrLaplSqrthName, M,P, k_iter), ros2File u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter) #, rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoLaplSqrthTh0, k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf(mrLaplSqrthTh0Name, M,P, k_iter), ros2File u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter) #, rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

plot for [M in "2 3"] for [P in "2 4"] sprintf(mrFoLaplh2,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdcLapl\~h^2-M%s-P%s-kiter%d', M,P, k_iter), ros2File u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcFo,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d', M, k_iter) #, rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

#plot for [M in "2 3"] for [P in "2 3 6 8"] sprintf(mrLapl0Fo,k_iter,M,P) u 5:3 w lp lc rgbcolor word(sdcCol,P+0) pt (M+1) title sprintf('mrsdc-M%s-P%s-kiter%d-Lapl_{0}', M,P, k_iter), "ros2_alpha_1e-3/accuracy.dat" u 5:3 w l lc rgbcolor 'orange' lw 4 t "ros2", for [M in "2 3"] sprintf(sdcLapl0Fo,k_iter,M) u 5:3 w lp lc rgbcolor word(sdcCol,1) pt (M+1) title sprintf('sdc-M%s-kiter%d-Lapl_{0}', M, k_iter) #, rkcFile u 5:3 w l dt rkcDash lw 2 lc rgbcolor rkcColor t rkcName,

}

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
