
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

BIN_V1="../bin/GNU-11.3.0/release/losslessImageVersion1"

mkdir -p binV1
mkdir -p bzip2
mkdir -p gzip
mkdir -p lzip
mkdir -p png
mkdir -p recV1


for ORG in ${IM_LIST[@]}
do
    STR="binV1/${ORG%.pgm}.bin"
    REC="recV1/$ORG"

    $BIN_V1 -e $ORG $STR
    $BIN_V1 -d      $STR $REC

    printf "diff: "
    diff $ORG $REC
    printf "\n\n"

    GZIP="gzip/${ORG}.gz"
    BZIP="bzip2/${ORG}.bz2"
    LZIP="lzip/${ORG}.lz"
    PNG="png/${ORG%.pgm}.png"
    gzip  -k -c $ORG > $GZIP
    lzip  -k -c $ORG > $LZIP
    bzip2 -k -c $ORG > $BZIP
    convert     $ORG $PNG
done


./genPlot.sh






