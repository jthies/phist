SMATS=`ls S*.mm`
DMATS=`ls D*.mm`
CMATS=`ls C*.mm`
ZMATS=`ls Z*.mm`

for mat in ${SMATS}
        do
                echo $mat
                ./mm2binCRS.py ${mat} 5
                mv out.crs $(basename "${mat}" .mm).bin
        done

for mat in ${DMATS}
        do
                echo $mat
                ./mm2binCRS.py ${mat} 6
                mv out.crs $(basename "${mat}" .mm).bin
        done

for mat in ${CMATS}
        do
                echo $mat
                ./mm2binCRS.py ${mat} 9
                mv out.crs $(basename "${mat}" .mm).bin
        done

for mat in ${ZMATS}
        do
                echo $mat
                ./mm2binCRS.py ${mat} 10
                mv out.crs $(basename "${mat}" .mm).bin
        done

./mm2binCRS.py mhd4800b.mm 6
mv out.crs mhd4800b.bin

./mm2binCRS.py mhd1280b.mm 10
mv out.crs mhd1280b.bin
