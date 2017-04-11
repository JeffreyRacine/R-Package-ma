#!/bin/sh

maver=`cat maver`
miver=${maver##*-}
maver=${maver%-*}

if [ $1 ] 
    then 
    let 'miver+=1' 
    echo $maver-$miver > maver
fi

echo updating DESCRIPTION
sed -i.bak -e 's/Version[:].*/'"Version: $maver-$miver/" -e 's/Date[:].*/'"Date: `date +%Y-%m-%d`/" DESCRIPTION

rm DESCRIPTION.bak


