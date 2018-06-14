      SUBROUTINE html_out(
     >                    where,it,iter,ivers)
c
c writes html elements into the inf-file to make it readable by a html-viewer
c called from the main program a beginning & end and before every major call.
c
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: it,iter
      CHARACTER(len=3),INTENT (IN) :: where
      CHARACTER(len=9),INTENT (IN) :: ivers

      INTEGER ind,itnext,i,j,itab
      CHARACTER(len=3)  mark,next(6)
      CHARACTER(len=14) name(6)
c
c some definitions first
c
      name(1) = 'Potential     ' ; next(1) = 'eig'
      name(2) = 'Eigenvalues   ' ; next(2) = 'fer'
      name(3) = 'Fermi Energy  ' ; next(3) = 'rho'
      name(4) = 'Charge Density' ; next(4) = 'tot'
      name(5) = 'Total Energy  ' ; next(5) = 'mix'
      name(6) = 'Mixing        ' ; next(6) = 'pot'
c
c from where did we get the call?
c
      mark='mid'
      IF (where.EQ.'top') THEN
        mark='top'
      ELSEIF (where.EQ.'pot') THEN
        ind = 1
      ELSEIF (where.EQ.'eig') THEN
        ind = 2
      ELSEIF (where.EQ.'fer') THEN
        ind = 3
      ELSEIF (where.EQ.'rho') THEN
        ind = 4
      ELSEIF (where.EQ.'tot') THEN
        ind = 5
      ELSEIF (where.EQ.'mix') THEN
        ind = 6
      ELSEIF (where.EQ.'bot') THEN
        mark='bot'
      ELSE
        WRITE (*,*) 'Do not call html_out with ',where
      ENDIF
      itnext = it+1
      IF (itnext.GT.iter) itnext = 1

      IF (mark.EQ.'top') THEN
c
c        write title section
c
        WRITE (16,*) '<HEAD> <TITLE>  f l a p w  version ',
     +                            ivers,' </TITLE> </HEAD>'
        WRITE (16,*) '<BODY>'
        WRITE (16,*) '<H2><FONT COLOR="#800000"> f l a p w  version ',
     +               ivers,' </FONT>  </H2> <P>'
        WRITE (16,*) '<H3><A HREF="#dum0"> Input Summary </A></H3> <P>'
        WRITE (16,*) '<table border>'
        DO itab = 1,6
          IF (mod(itab,2).EQ.1) WRITE (16,*) '<tr>'
          WRITE (16,*) '<td align=right><H3>',name(itab) 
          j = itab-1 ; IF (j.EQ.0) j=6
          DO i = 1,iter
            WRITE (16,7000) next(j),i,i
          ENDDO
          WRITE (16,*) '</H3>'
        ENDDO
        WRITE (16,*) '</table>'
 7000   FORMAT('<A HREF="#',a3,i1,'">',i1,'</A>')

        WRITE (16,8000) 'dum',0,'Input Summary ',0
        WRITE (16,8010) 'pot',1,'dum',0,0
        WRITE (16,*) '<PRE>'
      ELSEIF (mark.EQ.'mid') THEN
c
c       write a marker inside text
c
        WRITE (16,*) '</PRE>'
        WRITE (16,8000) where,it,name(ind),it
        IF (next(ind).EQ.'pot') THEN
          WRITE (16,8010) next(ind),itnext,where,itnext,itnext
        ELSE
          WRITE (16,8010) next(ind),it,where,itnext,itnext
        ENDIF
        WRITE (16,*) '<PRE>'
      ELSEIF (mark.EQ.'bot') THEN
c
c       at the end of file:
c
        WRITE (16,*) '</PRE>'
        WRITE (16,*) '</BODY>'
        WRITE (16,*) '</HTML>'
      ENDIF
 8000 FORMAT('<P><A NAME="',a3,i1,
     +       '"</A><B><FONT COLOR="#208000">  -------- ',
     +       a14,' # ',i1,' ----------</FONT></B>')
 8010 FORMAT('<DIV ALIGN=right> <A HREF="#top">top</A> <A HREF="#',a3,
     +       i1,'">next</A> <A HREF="#',a3,i1,'">it#',i1,'</A> </DIV>')

      RETURN
      END
