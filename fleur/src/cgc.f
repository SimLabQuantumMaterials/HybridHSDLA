      FUNCTION cgc(l1,l2,l3,m1,m2,m3)

      IMPLICIT NONE  
      INTEGER :: l1, l2, l3, m1, m2, m3
      REAL  :: two_l1p1, two_l1p2, l1pm3, l1pm3p1, l1mm3p1, l1mm3, cgc 

      IF (m3 /= m1 + m2) THEN
       cgc = 0.0
       RETURN
      END IF 
!     gb  m3 = m1 + m2
      two_l1p1 = 2 * l1 + 1
      two_l1p2 = 2 * l1 + 2
      l1pm3 = l1 + m3
      l1pm3p1 = l1 + m3 + 1
      l1mm3p1 = l1 - m3 + 1
      l1mm3 = l1 - m3 
      cgc = 0.0 
      IF (l3 == l1 + 1) THEN
          IF (m2 == 1) then
           cgc = sqrt( (l1pm3 * l1pm3p1) / (two_l1p1 * two_l1p2))   
          ELSEIF (m2 == 0) THEN
           cgc = sqrt( (l1mm3p1 * l1pm3p1) / (two_l1p1 * (l1 + 1)))
          ELSEIF (m2 == -1) THEN
           cgc = sqrt( (l1mm3 * l1mm3p1) / (two_l1p1 * two_l1p2))   
          END IF
      ELSE IF(l3 == l1 -1) THEN
          IF (m2 == 1) then
           cgc = sqrt( (l1mm3 * l1mm3p1) / (2.d0 * l1 * two_l1p1))   
          ELSEIF (m2 == 0) THEN
           cgc = -sqrt( (l1mm3 * l1pm3) / (l1 * (two_l1p1)))   
          ELSEIF (m2 == -1) THEN
           cgc = sqrt( (l1pm3p1 * l1pm3) / (2.0d0 * l1 * two_l1p1))   
          END IF
      END IF
      END FUNCTION cgc
