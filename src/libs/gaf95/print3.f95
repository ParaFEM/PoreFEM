      subroutine print3( k, fmt, val1, val2, val3 )
      parameter (MP = 3)
      character*(*) fmt
      real val1, val2, val3, tmp(MP)

      tmp(1) = val1
      tmp(2) = val2
      tmp(3) = val3

      call printv( k, fmt, tmp, MP )

      return
      end
