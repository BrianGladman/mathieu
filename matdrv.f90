program matdrv
    use param
    use mathieu

    implicit none
    integer   ndec, i, j, k, l, nex, kindd, kindq, icq, isq, minacc, narg, &
              lnum, izxi, ioprad, iopang, max, maxd, maxj, maxk, maxn, maxp, &
              maxterm, maxbp, maxkbp, maxlp, ismax
    real(knd) cm, q, arg1, darg, r, x1, z, xi

    real(knd), dimension(:), allocatable ::     arg
    real(knd), dimension(:), allocatable ::     mc1c, mc1dc, mc23c, mc23dc, &
                                                ms1c, ms1dc, ms23c, ms23dc
    integer, dimension(:), allocatable ::       mc1e, mc1de, mc23e, mc23de, &
                                                ms1e, ms1de, ms23e, ms23de, &
                                                acc_rc, acc_rs
    integer, dimension(:,:), allocatable ::     acc_acs
    real(knd), dimension (:,:), allocatable ::  ce, ced, se, sed
    character (len = 9) kstr, fstr, astr
    character (len = 256) tstr

    kindd = 8
    kindq = 16

    if(knd == kindd) then
        kstr = "e23.14"
        fstr = "f17.14"
    else if (knd == kindq) then
        kstr = "e39.30"
        fstr = "f17.14"   ! normally "f33.30"
    end if
    astr = "e23.14"

    ndec = precision(cm)
    nex = range(cm) - 1

!   open input and output files
    open(1, file='matfcn.dat')
    open(20 ,file='fort.20')
    open(30, file='fort.30')
    open(40 ,file='fort.40')
    open(50, file='fort.50')
    open(60, file='fort.60')

    read(1, *) lnum, ioprad, iopang, izxi, icq, isq

    if (icq == 1) read(1,*) q
    if (icq == 2) read(1,*) cm
    read(1, *) r
    if (iopang /= 0) read(1, *) arg1, darg, narg

    select case (icq)
    case (1)
        cm = sqrt(2.0e0_knd * abs(q))
        if (q < 0.0e0_knd) isq = -1
    case (2)
          q = cm * cm / 2.0e0_knd
          if (isq == -1) q = -q
    end select

    select case (izxi)
    case (1)
          x1 = 2.0e0_knd * (sinh(0.5e0_knd * z)) ** 2
          z = r
    case (2)
          x1 = r
          z = 0.0e0_knd
    end select
    if(iopang  ==  0) narg = 1

    maxj = 1
    if(ioprad /= 0) maxj = 1.5 * lnum + 4 * ndec + int(cm) + 105


    if(isq == 1) then
        ismax = int(0.85e0_knd * cm)
        if(x1 < 1.0e0_knd) ismax = 1.5 * ismax
        if(x1 < 0.1e0_knd .and. x1 >= 0.00001e0_knd) &
            ismax = ismax * int(3.0e0_knd ** (-1.0e0_knd - log10(x1)))
        if(x1 < 0.00001e0_knd) ismax = 81 * ismax

        maxj = max(maxj, ismax + int(cm / 3.14159 )+ 4 * ndec + int(cm) + 105)
    end if

    maxp = 2 * lnum + 2 * cm + 4 * ndec + 105
    maxn = maxj
    maxk = 1
    maxkbp = 1

    if (ioprad /= 0 .and. isq /= 1 .and. x1 /= 0.0e0_knd) then
        maxkbp = lnum + 4 * ndec + int(cm) + 103
        if(minacc > ndec - 2 .or. (minacc > ndec - 4 .and. cm >  5.0e0_knd) &
                             .or. (minacc > ndec - 6 .and. cm > 10.0e0_knd) &
                             .or. (minacc > ndec - 8 .and. cm > 15.0e0_knd) &
                             .or. cm > 20e0_knd) then
            maxk = 1.5 * cm + 100
            if(x1 < 1.0e0_knd .and. x1 >= 0.01e0_knd) &
                maxk = 40 / x1 + cm / sqrt(x1) + 1.5 * cm + 200
            if(x1 < 0.01e0_knd.and.x1 > 0.0e0_knd) &
                maxk = 34 / x1 + 1.5 * cm + 1.4 * cm / sqrt(x1) + 200
            maxterm = 10000000
            maxk = min(maxk, maxterm)
            maxk = max(maxk, maxj)
        end if
    end if

    maxd = max(maxj, maxn, maxk, maxkbp, maxp) / 2 + 1

    maxlp = lnum + 3
    if(isq == -1) maxlp = max(maxlp, ismax + int(cm /3.14159) + 3)

    xi = x1 + 1.0e0_knd

    if(ioprad /= 0 .and. izxi == 1) write(20, "(1x,' z = '," // kstr // ")") z

    if(ioprad /= 0 .and. izxi == 2) write(20, "(1x,'xi = '," // kstr // ")") xi

    if(icq == 1) then
        if (ioprad /= 0) write(20, "(1x,' q = '," // kstr // ")") q
        if (iopang /= 0) write(30, "(1x,' q = '," // kstr // ")") q
    end if

    if(icq == 2) then
        if (ioprad /= 0) then
            if (isq == 1)  write(20, "(1x,' c = '," // kstr // ")") cm
            if (isq == -1) write(20, "(1x,' c = i times'," // kstr // ")") cm
        end if
        if (iopang /= 0) then
            if (isq == 1)  write(30, "(1x,' c = '," // kstr // ")") cm
            if (isq == -1) write(30, "(1x,' c = i times'," // kstr // ")") cm
        end if
    end if

    allocate (arg(narg), acc_rc(lnum), acc_rs(lnum), acc_acs(lnum, narg))
    allocate (mc1c(lnum),mc1dc(lnum),mc23c(lnum),mc23dc(lnum), &
              ms1c(lnum),ms1dc(lnum),ms23c(lnum),ms23dc(lnum))
    allocate (mc1e(lnum),mc1de(lnum),mc23e(lnum),mc23de(lnum), &
              ms1e(lnum),ms1de(lnum),ms23e(lnum),ms23de(lnum))
    allocate (ce(lnum,narg),ced(lnum,narg), &
              se(lnum,narg),sed(lnum,narg))

    if(iopang /= 0) then
        do j = 1, narg
            arg(j) = arg1 + (j - 1) * darg
        end do
    end if

    call matfcn(lnum, ioprad, izxi, icq, isq, cm, r, iopang, narg, arg, &
                mc1c, mc1e, mc1dc, mc1de, mc23c, mc23e, mc23dc, mc23de, acc_rc, &
                ms1c, ms1e, ms1dc, ms1de, ms23c, ms23e, ms23dc, ms23de, acc_rs, &
                ce, ced, se, sed, acc_acs)

    do i = 1, lnum
        l = i - 1

        if(iopang  /=  0) write(30, 70) l
70      format(1x,'l = ',i6)

        if(isq /= -1) then
            tstr = "(1x,i5,2x,4(" // fstr // ",i6,2x), i2)"
            write(20, tstr) l,mc1c(i),mc1e(i),mc1dc(i),mc1de(i),mc23c(i),mc23e(i),mc23dc(i),mc23de(i),acc_rc(i)
            tstr = "(8x,4(" // fstr // ",i6,2x), i2)"
            if(l > 0) write(20, tstr) ms1c(i),ms1e(i),ms1dc(i),ms1de(i),ms23c(i),ms23e(i),ms23dc(i),ms23de(i),acc_rs(i)
        else
            tstr = "(1x,i5,2x,4(" // fstr // ",i7,2x), i2)"
            write(20, tstr) l,mc1c(i),mc1e(i),mc1dc(i),mc1de(i),mc23c(i),mc23e(i),mc23dc(i),mc23de(i),acc_rc(i)
            tstr = "(8x,4(" // fstr // ",i7,2x), i2)"
            if(l > 0) write(20, tstr) ms1c(i),ms1e(i),ms1dc(i),ms1de(i),ms23c(i),ms23e(i),ms23dc(i),ms23de(i),acc_rs(i)
        end if

       do k = 1,narg
            if(iopang == 1) then
                tstr = "(4x," // fstr // ",5x," // astr // ",2x," // astr // ",2x,i2)"
                write(30, tstr) arg(k),ce(i,k),se(i,k), acc_acs(i,k)
            else if(iopang == 2) then
                tstr = "(4x," // fstr // ",5x," // astr // ",2x," // astr // ",2x,/,26x," // astr // ",2x," // astr // ",2x,i2)"
                write(30, tstr) arg(k),ce(i,k),ced(i,k),se(i,k),sed(i,k), acc_acs(i,k)
            endif
        end do
    end do

    deallocate (mc1c,mc1dc,mc23c,mc23dc,ms1c,ms1dc,ms23c,ms23dc)
    deallocate (mc1e,mc1de,mc23e,mc23de,ms1e,ms1de,ms23e,ms23de)
    deallocate (ce,ced,se,sed)
    deallocate (acc_rc, acc_rs, acc_acs)
    close(60)
    close(50)
    close(40)
    close(30)
    close(20)
    close(1)

end program matdrv
