line=`ps -aef | grep virtuoso-t`
pid=`echo $line | cut -d " " -f 2`
kill -9 $pid
cd /usr/local/var/lib/virtuoso/db/
/usr/local/bin/virtuoso-t

