while read -r line; do echo ${#line}; done < eachlinecount.txt
# awk '{ print length($0); }' abc.txt
