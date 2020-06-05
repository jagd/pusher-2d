sed -i -e 's/\<double\>/long double/g' `find benchmark tests pusher -iname '*.cc' -o -iname '*.h'`
