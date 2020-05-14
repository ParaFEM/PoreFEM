      subroutine propst(p)
      real p(*)
      data zero/0.0/, half/0.5/, one/1.0/, two/2.0/, p15/1.5/,p25/2.5/,p35/3.5/
	if(p(3)>p15.and.p(3)<=p25)then
         p(5) = alog( one + p(2)*p(2)/(p(1)*p(1)) )	! var  of log-property
         p(4) = alog(p(1)) - half*p(5)			! mean of log-property
         p(5) = sqrt(p(5))				! sd   of log-property
      elseif( p(3)>p25.and.p(3)<=p35.and.p(6).eq.zero)then
            p(1) = half*(p(4)+p(5))
      endif
      end
