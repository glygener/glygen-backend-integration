C
C       To map the amino acid position PDB against Uniprot
C	pdb2uniprot 
C	Author: Nagarajan Pattabiraman
C	Version: 1.0
C	Date 01-Aug-2018
C

	dimension pid(100),amac(100),chid(100),ids(100)
	dimension dchid(100),ndbres(100),dbref(100),iuni(100)
	dimension lchid(100),lseq(100)
	dimension mdb(100),las(100)
	dimension imdif(100)
	dimension incmod(100),inclnk(100)
c
        character*80 line
	character*8 type,dbref,mdb,mdbp
	character*4 pid,pdbid1
	character*3 amac,las
	character*1 chid,dchid,lchid
C	
        idb = 0
        imd = 0
        ili = 0
	distmin = 1.35
	distmax = 1.50
C
100     continue

        read(5,10,end=101) line
C 	write(6,10)line
10      format(a)
c
C	Select only Uniprot Database in DBREF record from PDB
C
        if(line(1:5).eq.'DBREF'.and.line(27:29).eq.'UNP') then
            idb = idb + 1
            read(line,301,end=101) pdbid1,dchid(idb),ndbres(idb)
     1      ,dbref(idb),iuni(idb)
C           write(6,301) pdbid1,dchid(idb),ndbres(idb),dbref(idb)
C    1      ,iuni(idb)
301        format(7X,A4,1x,A1,1X,I4,15X,A8,14X,I5)
            go to 100
C
C	Read MODRES record form PDB
C
        else if(line(1:6).eq.'MODRES') then
            if(line(30:47).ne.'GLYCOSYLATION SITE') go to 100
            imd = imd + 1
            read(line,201,end=101) pid(imd),amac(imd),chid(imd),ids(imd)
201     format(7X,A4,1X,A3,1X,A1,1X,I4)
C            write(6,201) pid(imd),amac(imd),chid(imd),ids(imd)
             go to 100
C
C	Read LINK record from PDB
C
        else if(line(1:4).eq.'LINK') then
C
	     read(line,209) dist
209          format(73x,f5.2)
C	
C	Eleminate the LINK record with the last column lt distmin and 
C       gt than distmax
C
	     if(dist .lt. distmin .or. dist .gt. distmax) go to 100
C
C	Check the column containing ASN in LINK record
C
	     if(line(18:20).eq.'ASN') then
                ili = ili + 1
                read(line,401)	las(ili),lchid(ili),lseq(ili) 
C               write(6,404) lchid(ili),lseq(ili) 
401	format(17x,A3, 1X, A1,I4)
404	format('1 ','ASN',21x,A1,I4)
	        go to 100
	     else if(line(48:50).eq.'ASN') then
                ili = ili + 1
                read(line,402)	lchid(ili),lseq(ili) 
C               write(6,403) lchid(ili),lseq(ili) 
402	format(51x,A1,I4)
403	format('2 ASN',51x,A1,I4)
                go to 100
            endif
C
C	Check the column containing SER in LINK record
C
	     if(line(18:20).eq.'SER') then
                ili = ili + 1
                read(line,401)	las(ili),lchid(ili),lseq(ili) 
C               write(6,404) lchid(ili),lseq(ili) 
	        go to 100
	     else if(line(48:50).eq.'SER') then
                ili = ili + 1
                read(line,402)	lchid(ili),lseq(ili) 
C               write(6,403) lchid(ili),lseq(ili) 
                go to 100
            endif
C
C	Check the column containing THR in LINK record
C
	     if(line(18:20).eq.'THR') then
                ili = ili + 1
                read(line,401)	las(ili),lchid(ili),lseq(ili) 
C               write(6,404) lchid(ili),lseq(ili) 
	        go to 100
	     else if(line(48:50).eq.'THR') then
                ili = ili + 1
                read(line,402)	lchid(ili),lseq(ili) 
C               write(6,403) lchid(ili),lseq(ili) 
                go to 100
               endif
	    go to 100
        endif
c
	go to 100
101     continue
C
	if(idb.eq.0 ) then
           print *, 'NO DBREF IDB = 0', pdbid1
	   go to 1001
	endif
C
C	Keep only one chain from the MODRES record from PDB
C
	if(imd .eq. 0) go to 1002
	do 700 i = 1,imd
	   isel = 0 
	   do 701 j = 1,idb
	      if(dchid(j).eq.chid(i)) then
	         mdb(i) = dbref(j)
                 imdif(i) = iuni(j) - ndbres(j)
	         isel = 1
	         go to 867
	      endif
701	   continue
867	continue
	incmod(i) = isel
700	continue
C
	do 801 j = 1,imd
	   if(incmod(j).eq.0) go to 801
	   ieq = 0
	   do 800 i = 1,j
	      if(i.eq.j) go to 800
              if((ids(j).eq.ids(i)).and. (mdb(j).eq.mdb(i)))then
	         ieq = 1
                 go to 802
	      endif
800	  continue
802	  continue
          if(ieq.eq.0) then
	     if(amac(j).eq.'ASN') type = 'N-Linked'
	     if(amac(j).eq.'SER'.or.amac(j).eq.'THR') 
     1            type = 'O-Linked'
    	      print *,pdbid1,' ',amac(j),' ',chid(j),ids(j),mdb(j),
     1        ids(j)+imdif(j),imdif(j),type,' RESULTS MODRES '
	   endif
801	continue
	go to 1001
C
C	Keep only one chain from the LINK record from PDB
C
1002	continue
	do 1700 i = 1,ili
	   isel = 0
	   do 1701 j = 1,idb
	      if(dchid(j).eq.lchid(i)) then
      	         mdb(i) = dbref(j)
                 imdif(i) = iuni(j) - ndbres(j)
	         isel = 1
                 go to 1901
	      endif
1701	   continue
1901	continue
	inclnk(i) = isel
1700	continue
C
C	Keep only one chain from PDB
	do 901 j = 1,ili
	   if(inclnk(j).eq.0)  go to 901
	   ieq = 0
	   do 900 i = 1,j
	      if(i.eq.j) go to 900
              if((lseq(j).eq.lseq(i)).and. (mdb(j).eq.mdb(i)))then
	         ieq = 1
                 go to 902
	      endif
900	    continue
902	    continue
            if(ieq.eq.0) then
	       if(las(j).eq.'ASN') type = 'N-Linked'
	       if(las(j).eq.'SER'.or.las(j).eq.'THR') 
     1             type = 'O-Linked'
	       if(inclnk(j).eq.1) then
 	          print *,pdbid1,' ',las(j),' ',lchid(j),lseq(j),mdb(j),
     1                   lseq(j)+imdif(j),imdif(j),type,' RESULTS LINK '
                endif
	    endif
901	continue
C
1001	continue
        stop
        end
