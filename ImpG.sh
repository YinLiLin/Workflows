#! /bin/bash

G=$1; T=$2; W=$3

if [ ! $G ]
then
echo "Error: Please provide binary files!"
exit
fi

if [ ! $T ]
then
echo "Error: Please provide Z-score file!"
exit
fi

echo "Running ImpG..."
rm -rf ImpG_temp
mkdir ImpG_temp
echo "Transform file format by plink"
plink --bfile $G --allow-no-sex --recode beagle --out ImpG_temp/bgl > /dev/null
head -n 2 ImpG_temp/*.dat | tail -n 1 > ImpG_temp/ImpG.bgl
tail -n +4 ImpG_temp/*.dat >> ImpG_temp/ImpG.bgl
rm ImpG_temp/*.dat
printf 'SNP\tPos\ta0\ta1\n' > ImpG_temp/Map
paste -d'\t' <(cut -f 1-2 ImpG_temp/*.map) <(cut -f 3-4 ImpG_temp/*.map | tr "12" "AC") >> ImpG_temp/Map
mv ImpG_temp/Map ImpG_temp/*.map

if [ ! $W ]
then
W=1000000
fi

start=$(date +%s)

rm ImpG.log
mkdir ImpG_temp/maps
echo "Divide genome into windows"
GenMaps -m ImpG_temp/*.map -p ImpG_temp/maps/Wind -w $W -s ImpG_temp/chr_wind > ImpG.log

if [ $? -eq 0 ]
then
mkdir ImpG_temp/haps/
head -1 ImpG_temp/ImpG.bgl | tr ' ' '\n' | tail -n +3 | uniq | tr -s '\n' > ImpG_temp/samples
echo "Derive haps for windows"
GenHaps -d ImpG_temp/maps/ -s ImpG_temp/chr_wind -a ImpG_temp/ImpG.bgl -p ImpG_temp/samples -o ImpG_temp/haps/ >> ImpG.log
fi

if [ $? -eq 0 ]
then
mkdir ImpG_temp/typed/
awk '$3=="1"{$3="A"}{print}' $T | awk '$3=="2"{$3="C"}{print}' | awk '$4=="1"{$4="A"}{print}' | awk '$4=="2"{$4="C"}{print}' > ImpG_temp/y
echo "Divide typed SNPs within windows"
GenMaps-Typed -m ImpG_temp/maps/ -s ImpG_temp/chr_wind -y ImpG_temp/y -d ImpG_temp/typed/ -o ImpG_temp/chr_wind2 >> ImpG.log
fi

if [ $? -eq 0 ]
then
mkdir ImpG_temp/betas/
echo "Calculate Beta for windows"
ImpG-Summary-GenBeta-Chr -b ImpG-Summary-GenBeta -s ImpG_temp/chr_wind -d ImpG_temp/ >> ImpG.log
fi

if [ $? -eq 0 ]
then
mkdir ImpG_temp/imp/
echo "Impute Z values for missing SNPs"
ImpG-Summary-Chr -b ImpG-Summary -d ImpG_temp/ -s ImpG_temp/chr_wind >> ImpG.log
fi

if [ $? -eq 0 ]
then
echo "Merge Z values"
MergeZsc -d ImpG_temp/imp/ -s ImpG_temp/chr_wind -o all.impz >> ImpG.log
fi

end=$(date +%s)
echo "Total run time: "$(( $end - $start ))"s"
