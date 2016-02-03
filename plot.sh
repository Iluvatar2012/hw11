#!/usr/bin/gnuplot

# set terminal and output, use input filename as base for new filename
set terminal pngcairo size 1920,1080 enhanced font 'Verdana,30'

# set a grid for the graph
set grid

set label "dt = 0.1 * dx" center at -25,0.2

set output "0_to_3.png"
plot for [i=0:3] 'psi_'.i u 1:2 w l lw 2 t 'psi_{'.i.'}'

set output "4_to_7.png"
plot for [i=4:7] 'psi_'.i u 1:2 w l lw 2 t 'psi_{'.i.'}'

set output "8_to_10.png"
plot for [i=8:10] 'psi_'.i u 1:2 w l lw 2 t 'psi_{'.i.'}'

set output "0_to_10.png"
plot for [i=0:10] 'psi_'.i u 1:2 w l lw 2 t 'psi_{'.i.'}'
