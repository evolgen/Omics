var=`cat *.ffn`
for i in $var
do
if echo "$i" | grep -q ">"
then
str3=${i##*">"}
echo $i >$str3

else
echo $i >>$str3
fi

done
