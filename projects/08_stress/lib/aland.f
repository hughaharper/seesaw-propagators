c******************************************
c
      subroutine aland (rlat,rlon,ysw,iask)
c
c******************************************
      integer*2 mask,mas,n2
      character*80 namr
      logical ysw,y1st
      common/info/lcrt(2),lin(10),nin,lout(10),nout
      dimension mask(3600)
      save y1st,mask
      data namr /'\/user2\/sandwell\/alt\/data\/land1.bin'/
c  rlat is latitude in degrees
c  rlon is east longitude in degrees
c  ysw is returned as true if over land and false over the sea as
c    determined by 1 degree boxes stored as bits in file land#.bin which
c    is automatically read in to memory on the first call, requiring
c    3600 words.
c    anything above 81 or below -79 degrees latitude is defined as land.
c    iask if >0 then asks for land-sea filename
    1 format (' enter land/sea input file name')
    2 format (a80)
      data y1st/.true./
      if (y1st) then
         y1st= .false.
         if (iask.gt.0) then
            write (lcrt(2),1)
            read (lcrt(1),2) namr
            lup= lcrt(2)
         else
            lup= lcrt(2)
            endif
         open(38,file=namr,iostat=ios,err=9000,form='unformatted',
     1   status='old')
         read (38,iostat=ios,err=9010) (mask(j),j=1,3600)
         close (38,iostat=ios,err=9020)
         endif
 1100 continue
      ysw= .true.
      i4= rlat+79.
      if (i4.lt.0) return
      if (i4.ge.160) return
      blon= rlon + 360.
      n= amod(blon,360.) + 1.
c     the use of 2 word integers is necessary to use irbit as is:
      i4= i4*360 + n
      i4n= (i4-1)/16
      n2= i4 - i4n*16
      mas= mask(i4n+1)
      if (irbit(mas,n2).eq.0) ysw= .false.
      return
 9000 write (lup,9001) ios,namr
 9001 format (/' error no.',i6,' while opening file ',a6)
      stop
 9010 write (lup,9011) ios,namr,i,k
 9011 format (/' error no.',i6,' while reading file ',a6,5x,2i5)
      stop
 9020 write (lup,9021) ios,namr
 9021 format (/' error no.',i6,' while closing file ',a6)
      go to 1100
      end
      function irbit(iaddr,n)
c  returns as a 0 or 1 depending on the nth bit of array iaddr
      implicit integer*2 (i-n)
      integer*4 lcrt,lin,lout,nin,nout
      integer*2 bit(16),iaddr(1)
      common/info/lcrt(2),lin(10),nin,lout(10),nout
      data izero/0/
      data bit/o'100000',o'040000',o'020000',o'010000',o'004000',
     +  o'002000',o'001000',o'000400',o'000200',o'000100',
     +  o'000040',o'000020',o'000010',o'000004',o'000002',
     +  o'000001'/
      if(n.le.0) go to 100
      mword=(n-1)/16+1
      kbit=n-(mword-1)*16
      iwork=and(iaddr(mword),bit(kbit))
      if(iwork.ne.izero)iwork=1
      irbit=iwork
      return
  100 continue
      write(lcrt(2),101)n
  101 format(' tried to obtain value of bit ',i10,' in function',
     .  ' irbit -- stop')
      stop
      end
