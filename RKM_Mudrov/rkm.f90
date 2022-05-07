program rkm_method
    external rp
    common p
    real y(8)
1   type*, ' x, x9, h, eh, eps, h0, y(1), y(2), p?'
    accept*, x, x9, h, eh, eps, h0, y(1), y(2), p
    
    open(1, file='test.dat')
    
2   x1 = x
3   call rkm(2, x, h0, eps, y, rp)
    h0 = sign(h0, h)
    d = x1 + h - x
    if (abs(d) .le. eh*abs(x)) goto 4
    if ((d .lt. 0) .eq. (h .gt. 0)) h0 = d
    goto 3
4   write(1,'(3f18.10)'), x, y(1), y(2)
    if ((x .lt. x9) .eq. (h .gt. 0)) goto 2
    close(1)
    goto 1
end program rkm_method

subroutine rkm(n, x, h, e, y, rp)
    real y(8), z(8), f0(8), f(8), k1(8), k3(8)

    call rp(x, y, f0)

    do 11 i = 1, n
        z(i) = y(i)
11  end do

12  h3 = h/3.
    h4 = 4*h3

    do 13 i = 1, n
        k1(i) = h3*f0(i)
        y(i) = z(i) + k1(i)
13  end do

    call rp(x + h3, y, f)

    do 14 i = 1, n
        y(i) = z(i) + (k1(i) + h3*f(i))/2
14  end do

    call rp(x + h3, y, f)

    do 15 i = 1, n
        k3(i) = h*f(i)
        y(i) = z(i) + 0.375*(k1(i) + k3(i))
15  end do

    call rp(x + h/2, y, f)

    do 16 i = 1, n
        k1(i) = k1(i) + h4*f(i)
        y(i) = z(i) + 1.5*(k1(i) - k3(i))
16  end do

    call rp(x + h, y, f)

    r = 0.
    do 17 i = 1, n
        a = h3*f(i)
        y(i) = z(i) + (k1(i) + a)/2
        a = 2*k1(i) - 3.*k3(i) - a
        if (y(i) .ne. 0) a = a/y(i)
        if (abs(a) .gt. r) r = abs(a)
17  end do

    if (r .le. e) goto 18
    h = h/2
    goto 12
18  x = x + h
    if (32*r .lt. e) h = 2*h
    return
end subroutine rkm

subroutine rp(x, y, f)
    common p
    real y(2), f(2)

    f(1) = y(2)
    f(2) = p*(1.-y(1)**2)*y(2) - y(1)
    return
end subroutine rp
