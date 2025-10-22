
IM_LIST=(
    "kodim01.pgm"
    "kodim02.pgm"
    "kodim03.pgm"
    "kodim04.pgm"
    "kodim05.pgm"
    "kodim06.pgm"
    "kodim07.pgm"
    "kodim08.pgm"
    "kodim09.pgm"
    "kodim10.pgm"
    "kodim11.pgm"
    "kodim12.pgm"
    "kodim13.pgm"
    "kodim14.pgm"
    "kodim15.pgm"
    "kodim16.pgm"
    "kodim17.pgm"
    "kodim18.pgm"
    "kodim19.pgm"
    "kodim20.pgm"
    "kodim21.pgm"
    "kodim22.pgm"
    "kodim23.pgm"
    "kodim24.pgm"
)

mkdir -p dat
mkdir -p pdf

 totORG=0
totGZIP=0
totLZIP=0
totBZIP=0
 totPNG=0
totVER1=0

PDF_LIST=()
for ORG in ${IM_LIST[@]}
do
     PNG="png/${ORG%.pgm}.png"
    VER1="binV1/${ORG%.pgm}.bin"
    GZIP="gzip/${ORG}.gz"
    BZIP="bzip2/${ORG}.bz2"
    LZIP="lzip/${ORG}.lz"

    DATA="dat/${ORG%.pgm}.dat"

     fsORG=$(ls -l $ORG  | awk '{print $5}')
    fsGZIP=$(ls -l $GZIP | awk '{print $5}')
    fsLZIP=$(ls -l $LZIP | awk '{print $5}')
    fsBZIP=$(ls -l $BZIP | awk '{print $5}')
     fsPNG=$(ls -l $PNG  | awk '{print $5}')
    fsVER1=$(ls -l $VER1 | awk '{print $5}')
    echo   "${ORG%.pgm}" > $DATA
    printf "pgm\t%d\t%.2f\n"   $fsORG  $(echo "100.0*$fsORG /$fsORG" | bc -l ) >> $DATA
    printf "gzip\t%d\t%.2f\n"  $fsGZIP $(echo "100.0*$fsGZIP/$fsORG" | bc -l ) >> $DATA
    printf "lzip\t%d\t%.2f\n"  $fsLZIP $(echo "100.0*$fsLZIP/$fsORG" | bc -l ) >> $DATA
    printf "bzip2\t%d\t%.2f\n" $fsBZIP $(echo "100.0*$fsBZIP/$fsORG" | bc -l ) >> $DATA
    printf "png\t%d\t%.2f\n"   $fsPNG  $(echo "100.0*$fsPNG /$fsORG" | bc -l ) >> $DATA
    printf "VER1\t%d\t%.2f\n"  $fsVER1 $(echo "100.0*$fsVER1/$fsORG" | bc -l ) >> $DATA

     totORG=$(($totORG+$fsORG))
    totGZIP=$(($totGZIP+$fsGZIP))
    totLZIP=$(($totLZIP+$fsLZIP))
    totBZIP=$(($totBZIP+$fsBZIP))
     totPNG=$(($totPNG+$fsPNG))
    totVER1=$(($totVER1+$fsVER1))

    PDF="pdf/${ORG%.pgm}.pdf"
    GPL="pdf/${ORG%.pgm}.gp"
    echo -e "#!/usr/bin/gnuplot -persist" > $GPL
    echo -e "set terminal pdfcairo transparent enhanced font \"Arial,16\" size 4.50in, 3.20in" >> $GPL
    echo -e "set output '$PDF'" >> $GPL
    echo -e "set bar 1.000000 front" >> $GPL
    echo -e "set border 1 front lt black linewidth 1.000 dashtype solid" >> $GPL
    echo -e "set boxwidth 0.5 absolute" >> $GPL
    echo -e "set style fill   solid 1.00 border" >> $GPL
    echo -e "unset grid" >> $GPL
    echo -e "set xtics border in scale 1,0.5 nomirror norotate  autojustify" >> $GPL
    echo -e "set xtics  norangelimit " >> $GPL
    echo -e "set xtics   ()" >> $GPL
    echo -e "unset ytics" >> $GPL
    echo -e "set title \"${ORG%.pgm}\" " >> $GPL
    echo -e "set yrange [ 0.00000 : 105.000 ] noreverse nowriteback" >> $GPL
    echo -e "mycolor(x)=(x<=1?0:(x<=4?255:(x<=5?32768:16711680)))" >> $GPL
    echo -e "plot '$DATA' using 0:3:(mycolor(\$0)):xtic(1) with boxes lc rgb variable notitle, '' using 0:(\$3+5):(\$3) with labels notitle" >> $GPL

    gnuplot $GPL
    rm -f   $GPL

    PDF_LIST=( "${PDF_LIST[@]}" "$PDF" )
done


DATA="dat/total.dat"
echo "total" > $DATA
printf "pgm\t%d\t%.2f\n"   $totORG  $(echo "100.0*$totORG /$totORG" | bc -l ) >> $DATA
printf "gzip\t%d\t%.2f\n"  $totGZIP $(echo "100.0*$totGZIP/$totORG" | bc -l ) >> $DATA
printf "lzip\t%d\t%.2f\n"  $totLZIP $(echo "100.0*$totLZIP/$totORG" | bc -l ) >> $DATA
printf "bzip2\t%d\t%.2f\n" $totBZIP $(echo "100.0*$totBZIP/$totORG" | bc -l ) >> $DATA
printf "png\t%d\t%.2f\n"   $totPNG  $(echo "100.0*$totPNG /$totORG" | bc -l ) >> $DATA
printf "VER1\t%d\t%.2f\n"  $totVER1 $(echo "100.0*$totVER1/$totORG" | bc -l ) >> $DATA

PDF="pdf/total.pdf"
GPL="pdf/total.gp"
echo -e "#!/usr/bin/gnuplot -persist" > $GPL
echo -e "set terminal pdfcairo transparent enhanced font \"Arial,16\" size 4.50in, 3.20in" >> $GPL
echo -e "set output '$PDF'" >> $GPL
echo -e "set bar 1.000000 front" >> $GPL
echo -e "set border 1 front lt black linewidth 1.000 dashtype solid" >> $GPL
echo -e "set boxwidth 0.5 absolute" >> $GPL
echo -e "set style fill   solid 1.00 border" >> $GPL
echo -e "unset grid" >> $GPL
echo -e "set xtics border in scale 1,0.5 nomirror norotate  autojustify" >> $GPL
echo -e "set xtics  norangelimit " >> $GPL
echo -e "set xtics   ()" >> $GPL
echo -e "unset ytics" >> $GPL
echo -e "set title \"all 24 Kodak images\" " >> $GPL
echo -e "set yrange [ 0.00000 : 105.000 ] noreverse nowriteback" >> $GPL
echo -e "mycolor(x)=(x<=1?0:(x<=4?255:(x<=5?32768:16711680)))" >> $GPL
echo -e "plot '$DATA' using 0:3:(mycolor(\$0)):xtic(1) with boxes lc rgb variable notitle, '' using 0:(\$3+5):(\$3) with labels notitle" >> $GPL

gnuplot $GPL
rm -f   $GPL

PDF_LIST=( "$PDF" "${PDF_LIST[@]}" )


pdfjam --fitpaper true --rotateoversize true -o results.pdf --landscape "${PDF_LIST[@]}"






