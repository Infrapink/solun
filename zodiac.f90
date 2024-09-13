module sidereal
  implicit none

  ! Calculate right ascension and declination of a given star,
  ! taking into account its proper motion and the precession of the equinoxes
  
  ! Based on cpater 21 of *Astronomical Algorithms* by Jean Meeus
  ! 2nd Edition, Willman-Bell, 1991

contains
  
  subroutine propmot(jday, ra2000, dec2000, distance, rv, deltara, deltadec, answer)
    ! Calculate the effect of proper motion on a star's right ascension and declination
    
    real(8), intent(in) :: jday ! Julian Day we're interetested in
    real(8), intent(in) :: ra2000 ! right ascension at J2000.0
    real(8), intent(in) :: dec2000 ! declination at J2000.0
    real(8), intent(in) :: distance ! distance from the sun, in parsecs
    real(8), intent(in) :: rv ! radial velocity, in parsecs per year
    real(8), intent(in) :: deltara ! right ascension component of proper motion, in seconds of arc. This has to be looked up
    real(8), intent(in) :: deltadec ! declination component of proper motion, in seconds of time. This has to be looked up.
    real(8), dimension(2), intent(out) :: answer
    
    real(8) :: ra ! right ascension at the time of interest
    real(8) :: dec ! declination at the time of interest

    real(8) :: t
    real(8) :: u
    
    real(8) :: x
    real(8) :: y
    real(8) :: z

    real(8) :: xdelta
    real(8) :: ydelta
    real(8) :: zdelta

    real(8) :: xprime
    real(8) :: yprime
    real(8) :: zprime

    real(8) :: pi
    real(8) :: d2r ! convert degrees to radians
    real(8) :: r2d ! convert radians to degrees

    pi = 4.0 * atan(1.0)
    d2r = pi / 180.0 ! convert degrees to radians
    r2d = 180.0 / pi ! convert radians to degrees

    x = distance * cos(dec2000 * d2r) * cos(ra2000 * d2r)
    y = distance * cos(dec2000 * d2r) * sin(ra2000 * d2r)
    z = distance * sin(dec2000 * d2r)
    !print *, x, y, z

    xdelta = ((x / distance) * rv) - (z * deltadec * cos(ra2000 * d2r)) - (y * deltara)
    ydelta = ((y / distance) * rv) - (z * deltadec * sin(ra2000 * d2r)) + (x * deltara)
    zdelta = ((z / distance) * rv) - (distance * deltadec * cos(dec2000 * d2r))

    t = (jday - 2451545.0) / 365.25 ! julian years since the year 2000 began at Greenwich Observatory

    xprime = x + (t * xdelta)
    yprime = y + (t * ydelta)
    zprime = z + (t * zdelta)

    ra = atan(yprime / xprime) * r2d

    u = sqrt((xprime * xprime) + (yprime * yprime))
    dec = atan(zprime / u) * r2d

    !print *, ra
    !print *, dec
    answer(1) = ra
    answer(2) = dec
  end subroutine propmot

  subroutine precession(jday, ra2000, dec2000, distance, rv, deltara, deltadec, answer)
    ! apply the effect of precession to obtain
    ! the actual right ascension and declination

    real(8), intent(in) :: jday ! Julian Day in question
    real(8), intent(in) :: ra2000 ! right ascension at J2000.0, in DEGREES
    real(8), intent(in) :: dec2000 ! declination at J2000.0, in DEGREES
    real(8), intent(in) :: distance ! distance from the sun, in PARSECS
    real(8), intent(in) :: rv ! radial velocity, in parsecs per year
    real(8), intent(in) :: deltara ! right ascension component of proper motion, in ARCSECONDS. This has to be looked up
    real(8), intent(in) :: deltadec ! declination component of proper motion,m in ARCSECONDS. This has to be looked up.
    real(8), dimension(2), intent(out) :: answer ! RA and Dec for use in calculations

    real(8) :: pi
    real(8) :: d2r ! convert degrees to radians
    real(8) :: r2d ! convert radians to degrees

    real(8) :: t ! Julian centuries since J2000.0
    real(8) :: zeta
    real(8) :: z
    real(8) :: theta

    real(8) :: a ! placeholder
    real(8) :: b ! placeholder
    real(8) :: c ! placeholder

    real(8), dimension(2) :: radec ! RA and Dec after accounting for proper motion but before accounting for precession
    real(8) :: ra ! right ascension after taking proper motion into account
    real(8) :: dec ! declination after taking proper motion into account

    pi = 4.0 * atan(1.0)
    d2r = pi / 180.0
    r2d = 180.0 / pi

    !print *, jday, ra2000, dec2000, distance, rv, deltara, deltadec
    !if (distance /= 0) then
    call propmot(jday, ra2000, dec2000, distance, rv, deltara, deltadec, radec)
    !end if
    ra = radec(1)
    dec = radec(2)
    !print *, radec
    print *, ra
    print *, dec

    t = (jday - 2451545.0) / 36525.0

    zeta = (2306.2181 * t) + (0.30188 * t * t) + (0.017998 * (t ** 3))
    z = (2306.2181 * t) + (1.09468 * t * t) + (0.018203 * (t ** 3))
    theta = (2004.3109 * t) + (0.42665 * t * t) + (0.041833 * (t ** 3))
    !print *, zeta
    !print *, z
    !print *, theta

    ! zeta, z, and theta are in ARCSECONDS, and so need to be converted into radians for the next bit

    a = cos(d2r * dec) * sin(d2r * (ra + (zeta / 3600.0)))    
    b = cos(d2r * (theta / 3600.0))
    b = b * cos(d2r * dec)
    b = b * cos(d2r * (ra + (zeta / 3600.0)))
    b = b - (sin(d2r * (theta / 3600.0)) * sin(d2r * dec))
    c = sin(d2r * (theta / 3600.0))
    c = c * cos(d2r * dec)
    c = c * cos(d2r * (ra + (zeta / 3600.0)))
    c = c + (cos(d2r * (theta / 3600.0)) * sin(d2r * dec))

    !print *, a, b, c
    !print *, a/b
    !print *, atan(a/b)

    answer(1) = z + (r2d * atan(a / b)) ! right ascension, in degrees, taking both proper motion and precession into accont
    answer(2) = r2d * asin(c) ! declination, in degrees, taking both proper motion and precession into account
    !print *, answer(2)    
    !answer(1) = mod(answer(1), 360.0)
    !answer(2) = mod(answer(2), 360.0)
    !print *, answer
    !print *, "precession"
    !print *, "" ! For some reason, the numbers go weird if you don't print a blank line here. I don't know why.
    
  end subroutine precession

  subroutine nutation(jday, nut)
    ! Calculate the nutation of the obliquity to the ecplitic
    real(8), intent(in) :: jday ! Julian Day in question
    real(8), dimension(2), intent(out) :: nut ! Nutation to the ecliptic and of longitude

    real(8) :: T ! Julian Centuries since J2000.0
    real(8) :: D ! mean elonation of the moon from the sun
    real(8) :: M ! mean anomaly of the sun
    real(8) :: Mprime ! mean anomaly of the moon
    real(8) :: F ! moon's argument of latitude
    real(8) :: omega ! longitude of the ascending node of the moon's mean orbit on the ecliptic

    real(8), dimension(64,5) :: args
    real(8), dimension(64,2) :: psi_coeffs
    real(8), dimension(64,2) :: eps_coeffs
    real(8) :: epsilon0 !mean obliquity of the ecliptic, in degrees
    real(8) :: delta_epsilon ! variation in the obliquity of the ecliptic, in arcseconds
    real(8) :: delta_psi ! nutation in longitude, in arcseconds
    real(8) :: U ! Julian Decamillennia since J2000.0

    integer :: i
    real(8) :: a ! placeholder

    real(8) :: pi
    real(8) :: d2r ! convert degrees to radians
    real(8) :: r2d ! convert radians to degrees

    pi = 4.0 * atan(1.0)
    d2r = 180.0 / pi
    r2d = pi / 180.0

    T = (jday - 2451545.0) / 36525.0
    D = 297.85036 + (445267.111480 * T) - (0.00191427 * T * T) + ((T ** 3) / 189242.0)
    M = 357.52772 + (35999.0503407 * T) - (0.0001603 * T * T) - ((T ** 3) / 300000.0)
    Mprime = 134.96298 + (477198.867398 * T) + (0.0086972 * T * T) + ((T ** 3) / 56250.0)
    F = 93.27191 + (483202.017538 * T) - (0.0036825 * T * T) + ((T ** 3) / 327270.0)
    omega = 125.04452 - (1934.136261 * T) + (0.0020708 * T * T) + ((T ** 3) / 450000.0)

    args(1,1) = 0
    args(1,2) = 0
    args(1,3) = 0
    args(1,4) = 0
    args(1,5) = 1
    args(2,1) = -2
    args(2,2) = 0
    args(2,3) = 0
    args(2,4) = 2
    args(2,5) = 2
    args(3,1) = 0
    args(3,2) = 0
    args(3,3) = 0
    args(3,4) = 2
    args(3,5) = 2
    args(4,1) = 0
    args(4,2) = 0
    args(4,3) = 0
    args(4,4) = 0
    args(4,5) = 2
    args(5,1) = 0
    args(5,2) = 1
    args(5,3) = 0
    args(5,4) = 0
    args(5,5) = 0
    args(6,1) = 0
    args(6,2) = 0
    args(6,3) = 1
    args(6,4) = 0
    args(6,5) = 0
    args(7,1) = -2
    args(7,2) = 1
    args(7,3) = 0
    args(7,4) = 2
    args(7,5) = 2
    args(8,1) = 0
    args(8,2) = 0
    args(8,3) = 0
    args(8,4) = 2
    args(8,5) = 1
    args(9,1) = 0
    args(9,2) = 0
    args(9,3) = 1
    args(9,4) = 2
    args(9,5) = 2
    args(10,1) = -2
    args(10,2) = -1
    args(10,3) = 0
    args(10,4) = 2
    args(10,5) = 2
    args(11,1) = -2
    args(11,2) = 0
    args(11,3) = 1
    args(11,4) = 0
    args(11,5) = 0
    args(12,1) = -2
    args(12,2) = 0
    args(12,3) = 0
    args(12,4) = 2
    args(12,5) = 1
    args(13,1) = 0
    args(13,2) = 0
    args(13,3) = -1
    args(13,4) = 2
    args(13,5) = 2
    args(14,1) = 2
    args(14,2) = 0
    args(14,3) = 0
    args(14,4) = 0
    args(14,5) = 0
    args(15,1) = 0
    args(15,2) = 0
    args(15,3) = 1
    args(15,4) = 0
    args(15,5) = 1
    args(16,1) = 2
    args(16,2) = 0
    args(16,3) = -1
    args(16,4) = 2
    args(16,5) = 2
    args(17,1) = 0
    args(17,2) = 0
    args(17,3) = -1
    args(17,4) = 0
    args(17,5) = 1
    args(18,1) = 0
    args(18,2) = 0
    args(18,3) = 1
    args(18,4) = 2
    args(18,5) = 1
    args(19,1) = -2
    args(19,2) = 0
    args(19,3) = 2
    args(19,4) = 0
    args(19,5) = 0
    args(20,1) = 0
    args(20,2) = 0
    args(20,3) = -2
    args(20,4) = 2
    args(20,5) = 1
    args(21,1) = 2
    args(21,2) = 0
    args(21,3) = 0
    args(21,4) = 2
    args(21,5) = 2
    args(22,1) = 0
    args(22,2) = 0
    args(22,3) = 2
    args(22,4) = 2
    args(22,5) = 2
    args(23,1) = 0
    args(23,2) = 0
    args(23,3) = 2
    args(23,4) = 0
    args(23,5) = 0
    args(24,1) = -2
    args(24,2) = 0
    args(24,3) = 1
    args(24,4) = 2
    args(24,5) = 2
    args(25,1) = 0
    args(25,2) = 0
    args(25,3) = 0
    args(25,4) = 2
    args(25,5) = 0
    args(26,1) = -2
    args(26,2) = 0
    args(26,3) = 0
    args(26,4) = 2
    args(26,5) = 0
    args(27,1) = 0
    args(27,2) = 0
    args(27,3) = -1
    args(27,4) = 2
    args(27,5) = 1
    args(28,1) = 0
    args(28,2) = 0
    args(28,3) = -1
    args(28,4) = 2
    args(28,5) = 1
    args(29,1) = 0
    args(29,2) = 2
    args(29,3) = 0
    args(29,4) = 0
    args(29,5) = 0
    args(30,1) = 2
    args(30,2) = 0
    args(30,3) = -1
    args(30,4) = 0
    args(30,5) = 1
    args(31,1) = -2
    args(31,2) = 2
    args(31,3) = 0
    args(31,4) = 2
    args(31,5) = 2
    args(32,1) = 0
    args(32,2) = 1
    args(32,3) = 0
    args(32,4) = 0
    args(32,5) = 1
    args(33,1) = -2
    args(33,2) = 0
    args(33,3) = 1
    args(33,4) = 0
    args(33,5) = 1
    args(34,1) = 0
    args(34,2) = -1
    args(34,3) = 0
    args(34,4) = 0
    args(34,5) = 1
    args(35,1) = 0
    args(35,2) = 0
    args(35,3) = 2
    args(35,4) = -2
    args(35,5) = 0
    args(36,1) = 2
    args(36,2) = 0
    args(36,3) = -1
    args(36,4) = 2
    args(36,5) = 1
    args(37,1) = 2
    args(37,2) = 0
    args(37,3) = 1
    args(37,4) = 2
    args(37,5) = 2
    args(38,1) = 0
    args(38,2) = 1
    args(38,3) = 0
    args(38,4) = 2
    args(38,5) = 2
    args(39,1) = -2
    args(39,2) = 1
    args(39,3) = 1
    args(39,4) = 0
    args(39,5) = 0
    args(40,1) = 0
    args(40,2) = -1
    args(40,3) = 0
    args(40,4) = 2
    args(40,5) = 2
    args(41,1) = 2
    args(41,2) = 0
    args(41,3) = 0
    args(41,4) = 2
    args(41,5) = 1
    args(42,1) = 2
    args(42,2) = 0
    args(42,3) = 1
    args(42,4) = 0
    args(42,5) = 0
    args(43,1) = -2
    args(43,2) = 0
    args(43,3) = 2
    args(43,4) = 2
    args(43,5) = 2
    args(44,1) = -2
    args(44,2) = 0
    args(44,3) = 1
    args(44,4) = 2
    args(44,5) = 1
    args(45,1) = 2
    args(45,2) = 0
    args(45,3) = -2
    args(45,4) = 0
    args(45,5) = 1
    args(46,1) = 2
    args(46,2) = 0
    args(46,3) = 0
    args(46,4) = 0
    args(46,5) = 1
    args(47,1) = 0
    args(47,2) = -1
    args(47,3) = 1
    args(47,4) = 0
    args(47,5) = 0
    args(48,1) = -2
    args(48,2) = -1
    args(48,3) = 0
    args(48,4) = 2
    args(48,5) = 1
    args(49,1) = -2
    args(49,2) = 0
    args(49,3) = 0
    args(49,4) = 0
    args(49,5) = 1
    args(50,1) = 0
    args(50,2) = 0
    args(50,3) = 2
    args(50,4) = 2
    args(50,5) = 1
    args(51,1) = -2
    args(51,2) = 0
    args(51,3) = 2
    args(51,4) = 0
    args(51,5) = 1
    args(52,1) = -2
    args(52,2) = 1
    args(52,3) = 0
    args(52,4) = 2
    args(52,5) = 1
    args(53,1) = 0
    args(53,2) = 0
    args(53,3) = 1
    args(53,4) = -2
    args(53,5) = 0
    args(54,1) = -1
    args(54,2) = 0
    args(54,3) = 1
    args(54,4) = 0
    args(54,5) = 0
    args(55,1) = -2
    args(55,2) = 1
    args(55,3) = 0
    args(55,4) = 0
    args(55,5) = 0
    args(56,1) = 1
    args(56,2) = 0
    args(56,3) = 0
    args(56,4) = 0
    args(56,5) = 0
    args(57,1) = 0
    args(57,2) = 0
    args(57,3) = 1
    args(57,4) = 2
    args(57,5) = 0
    args(58,1) = 0
    args(58,2) = 0
    args(58,3) = -2
    args(58,4) = 2
    args(58,5) = 2
    args(59,1) = -1
    args(59,2) = -1
    args(59,3) = 1
    args(59,4) = 0
    args(59,5) = 0
    args(60,1) = 0
    args(60,2) = 1
    args(60,3) = 1
    args(60,4) = 0
    args(60,5) = 0
    args(61,1) = 0
    args(61,2) = -1
    args(61,3) = 1
    args(61,4) = 2
    args(61,5) = 2
    args(62,1) = 2
    args(62,2) = -1
    args(62,3) = -1
    args(62,4) = 2
    args(62,5) = 2
    args(63,1) = 0
    args(63,2) = 0
    args(63,3) = 3
    args(63,4) = 2
    args(63,5) = 2
    args(64,1) = 2
    args(64,2) = -1
    args(64,3) = 0
    args(64,4) = 2
    args(64,5) = 2

    psi_coeffs(1,1) = 171996
    psi_coeffs(1, 2) = -174.27
    psi_coeffs(2, 1) = -13187
    psi_coeffs(2, 2) = -1.67
    psi_coeffs(3, 1) = -2274
    psi_coeffs(3, 2) = -0.27
    psi_coeffs(4, 1) = 2062
    psi_coeffs(4, 2) = 0.27
    psi_coeffs(5, 1) = 1426
    psi_coeffs(5, 2) = -3.47
    psi_coeffs(6, 1) = 712
    psi_coeffs(6, 2) = 0.17
    psi_coeffs(7, 1) = -517
    psi_coeffs(7, 2) = 1.27
    psi_coeffs(8, 1) = -368
    psi_coeffs(8, 2) = -0.47
    psi_coeffs(9, 1) = -301
    psi_coeffs(9, 2) = 0
    psi_coeffs(10, 1) = 217
    psi_coeffs(10, 2) = -0.57
    psi_coeffs(11, 1) = -158
    psi_coeffs(11, 2) = 0
    psi_coeffs(12, 1) = 129
    psi_coeffs(12, 2) = 0.17
    psi_coeffs(13, 1) = 123
    psi_coeffs(13, 2) = 0
    psi_coeffs(14, 1) = 63
    psi_coeffs(14, 2) = 0
    psi_coeffs(15, 1) = 63
    psi_coeffs(15, 2) = 0.17
    psi_coeffs(16, 1) = -59
    psi_coeffs(16, 2) = 0
    psi_coeffs(17, 1) = -58
    psi_coeffs(17, 2) = -0.17
    psi_coeffs(18, 1) = -51
    psi_coeffs(18, 2) = 0
    psi_coeffs(19, 1) = 48
    psi_coeffs(19, 2) = 0
    psi_coeffs(20, 1) = 46
    psi_coeffs(20, 2) = 0
    psi_coeffs(21, 1) = -38
    psi_coeffs(21, 2) = 0
    psi_coeffs(22, 1) = -31
    psi_coeffs(22, 2) = 0
    psi_coeffs(23, 1) = 29
    psi_coeffs(23, 2) = 0
    psi_coeffs(24, 1) = 29
    psi_coeffs(24, 2) = 0
    psi_coeffs(25, 1) = 26
    psi_coeffs(25, 2) = 0
    psi_coeffs(26, 1) = -22
    psi_coeffs(26, 2) = 0
    psi_coeffs(27, 1) = 21
    psi_coeffs(27, 2) = 0
    psi_coeffs(28, 1) = 17
    psi_coeffs(28, 2) = 0.17
    psi_coeffs(29, 1) = 16
    psi_coeffs(29, 2) = 0
    psi_coeffs(30, 1) = -16
    psi_coeffs(30, 2) = 0.17
    psi_coeffs(31, 1) = -15
    psi_coeffs(31, 2) = 0
    psi_coeffs(32, 1) = -13
    psi_coeffs(32, 2) = 0
    psi_coeffs(33, 1) = -12
    psi_coeffs(33, 2) = 0
    psi_coeffs(34, 1) = 11
    psi_coeffs(34, 2) = 0
    psi_coeffs(35, 1) = -10
    psi_coeffs(35, 2) = 0
    psi_coeffs(36, 1) = -8
    psi_coeffs(36, 2) = 0
    psi_coeffs(37, 1) = 7
    psi_coeffs(37, 2) = 0
    psi_coeffs(38, 1) = -7
    psi_coeffs(38, 2) = 0
    psi_coeffs(39, 1) = -7
    psi_coeffs(39, 2) = 0
    psi_coeffs(40, 1) = -7
    psi_coeffs(40, 2) = 0
    psi_coeffs(41, 1) = 6
    psi_coeffs(41, 2) = 0
    psi_coeffs(42, 1) = 6
    psi_coeffs(42, 2) = 0
    psi_coeffs(43, 1) = 6
    psi_coeffs(43, 2) = 0
    psi_coeffs(44, 1) = -6
    psi_coeffs(44, 2) = 0
    psi_coeffs(45, 1) = -6
    psi_coeffs(45, 2) = 0
    psi_coeffs(46, 1) = 5
    psi_coeffs(46, 2) = 0
    psi_coeffs(47, 1) = -5
    psi_coeffs(47, 2) = 0
    psi_coeffs(48, 1) = -5
    psi_coeffs(48, 2) = 0
    psi_coeffs(49, 1) = -5
    psi_coeffs(49, 2) = 0
    psi_coeffs(50, 1) = 4
    psi_coeffs(50, 2) = 0
    psi_coeffs(51, 1) = 4
    psi_coeffs(51, 2) = 0
    psi_coeffs(52, 1) = 4
    psi_coeffs(52, 2) = 0
    psi_coeffs(53, 1) = -4
    psi_coeffs(53, 2) = 0
    psi_coeffs(54, 1) = -4
    psi_coeffs(54, 2) = 0
    psi_coeffs(55, 1) = -4
    psi_coeffs(55, 2) = 0
    psi_coeffs(56, 1) = 3
    psi_coeffs(56, 2) = 0
    psi_coeffs(57, 1) = -3
    psi_coeffs(57, 2) = 0
    psi_coeffs(58, 1) = -3
    psi_coeffs(58, 2) = 0
    psi_coeffs(59, 1) = -3
    psi_coeffs(59, 2) = 0
    psi_coeffs(60, 1) = -3
    psi_coeffs(60, 2) = 0
    psi_coeffs(61, 1) = -3
    psi_coeffs(61, 2) = 0
    psi_coeffs(62, 1) = -3
    psi_coeffs(62, 2) = 0
    psi_coeffs(63, 1) = -3
    psi_coeffs(63, 2) = 0
    psi_coeffs(64, 1) = -3
    psi_coeffs(64, 2) = 0
    
    eps_coeffs(1,1) = 92025
    eps_coeffs(1,2) = 8.9
    eps_coeffs(2,1) = 5736
    eps_coeffs(2,2) = 3.1
    eps_coeffs(3,1) = 977
    eps_coeffs(3,2) = 0.5
    eps_coeffs(4,1) = -895
    eps_coeffs(4,2) = 0.5
    eps_coeffs(5,1) = 54
    eps_coeffs(5,2) = 0.17
    eps_coeffs(6,1) = -7
    eps_coeffs(6,2) = 0
    eps_coeffs(7,1) = 224
    eps_coeffs(7,2) = -0.6
    eps_coeffs(8,1) = 200
    eps_coeffs(8,2) = 0
    eps_coeffs(9,1) = 129
    eps_coeffs(9,2) = -0.1
    eps_coeffs(10,1) = -95
    eps_coeffs(10,2) = 0.3
    eps_coeffs(11,1) = 0
    eps_coeffs(11,2) = 0
    eps_coeffs(12,1) = -70
    eps_coeffs(12,2) = 0
    eps_coeffs(13,1) = -53
    eps_coeffs(13,2) = 0
    eps_coeffs(14,1) = 0
    eps_coeffs(14,2) = 0
    eps_coeffs(15,1) = -33
    eps_coeffs(15,2) = 0
    eps_coeffs(16,1) = 26
    eps_coeffs(16,2) = 0
    eps_coeffs(17,1) = 32
    eps_coeffs(17,2) = 0
    eps_coeffs(18,1) = 27
    eps_coeffs(18,2) = 0
    eps_coeffs(19,1) = 0
    eps_coeffs(19,2) = 0
    eps_coeffs(20,1) = -24
    eps_coeffs(20,2) = 0
    eps_coeffs(21,1) = 16
    eps_coeffs(21,2) = 0
    eps_coeffs(22,1) = 13
    eps_coeffs(22,2) = 0
    eps_coeffs(23,1) = 0
    eps_coeffs(23,2) = 0
    eps_coeffs(24,1) = -12
    eps_coeffs(24,2) = 0
    eps_coeffs(25,1) = 0
    eps_coeffs(25,2) = 0
    eps_coeffs(26,1) = 0
    eps_coeffs(26,2) = 0
    eps_coeffs(27,1) = -10
    eps_coeffs(27,2) = 0
    eps_coeffs(28,1) = 0
    eps_coeffs(28,2) = 0
    eps_coeffs(29,1) = -8
    eps_coeffs(29,2) = 0
    eps_coeffs(30,1) = 7
    eps_coeffs(30,2) = 0
    eps_coeffs(31,1) = 9
    eps_coeffs(31,2) = 0
    eps_coeffs(32,1) = 7
    eps_coeffs(32,2) = 0
    eps_coeffs(33,1) = 6
    eps_coeffs(33,2) = 0
    eps_coeffs(34,1) = 0
    eps_coeffs(34,2) = 0
    eps_coeffs(35,1) = 5
    eps_coeffs(35,2) = 0
    eps_coeffs(36,1) = 3
    eps_coeffs(36,2) = 0
    eps_coeffs(37,1) = -3
    eps_coeffs(37,2) = 0
    eps_coeffs(38,1) = 0
    eps_coeffs(38,2) = 0
    eps_coeffs(39,1) = 3
    eps_coeffs(39,2) = 0
    eps_coeffs(40,1) = 3
    eps_coeffs(40,2) = 0
    eps_coeffs(41,1) = 0
    eps_coeffs(41,2) = 0
    eps_coeffs(42,1) = -3
    eps_coeffs(42,2) = 0
    eps_coeffs(43,1) = -3
    eps_coeffs(43,2) = 0
    eps_coeffs(44,1) = 3
    eps_coeffs(44,2) = 0
    eps_coeffs(45,1) = 3
    eps_coeffs(45,2) = 0
    eps_coeffs(46,1) = 0
    eps_coeffs(46,2) = 0
    eps_coeffs(47,1) = 3
    eps_coeffs(47,2) = 0
    eps_coeffs(48,1) = 3
    eps_coeffs(48,2) = 0
    eps_coeffs(49,1) = 3
    eps_coeffs(49,2) = 0
    eps_coeffs(50,1) = 0
    eps_coeffs(50,2) = 0
    eps_coeffs(51,1) = 0
    eps_coeffs(51,2) = 0
    eps_coeffs(52,1) = 0
    eps_coeffs(52,2) = 0
    eps_coeffs(53,1) = 0
    eps_coeffs(53,2) = 0
    eps_coeffs(54,1) = 0
    eps_coeffs(54,2) = 0
    eps_coeffs(55,1) = 0
    eps_coeffs(55,2) = 0
    eps_coeffs(56,1) = 0
    eps_coeffs(56,2) = 0
    eps_coeffs(57,1) = 0
    eps_coeffs(57,2) = 0
    eps_coeffs(58,1) = 0
    eps_coeffs(58,2) = 0
    eps_coeffs(59,1) = 0
    eps_coeffs(59,2) = 0
    eps_coeffs(60,1) = 0
    eps_coeffs(60,2) = 0
    eps_coeffs(61,1) = 0
    eps_coeffs(61,2) = 0
    eps_coeffs(62,1) = 0
    eps_coeffs(62,2) = 0
    eps_coeffs(63,1) = 0
    eps_coeffs(63,2) = 0
    eps_coeffs(64,1) = 0
    eps_coeffs(64,2) = 0

    delta_epsilon = 0.0
    delta_psi = 0.0
    do i = 1, 64
       a = (eps_coeffs(i,1) + (eps_coeffs(i,2) * T))
       a = a * cos(d2r * ((args(i,1) * D) + (args(i,2) * M) + (args(i,3) * Mprime) + (args(i,4) * F) + (args(i,5) * omega)))
       delta_epsilon = delta_epsilon + a

       a =  args(i,1) * D
       a = a + (args(i,2) * M)
       a = a + (args(i,3) * Mprime)
       a = a + (args(i,4) * F)
       a = a + (args(i,5) * omega)
       a = sin(d2r * a)
       a = a * (eps_coeffs(i,1) + (eps_coeffs(i,2) * T))
    end do

    ! Convert delta_epsilon and delta_psi into arcseconds
    delta_epsilon = delta_epsilon / 10000.0
    delta_psi = delta_psi / 1000.0 

    U = T / 100.0
    epsilon0 = (23.0 * 360.0) + (26.0 * 60.0) + 21.448 ! epsilon0 is now in seconds
    epsilon0 = epsilon0 - (4680.93 * U)
    epsilon0 = epsilon0 - (1.55 * (U ** 2))
    epsilon0 = epsilon0 + (1999.25 * (U ** 3))
    epsilon0 = epsilon0 - (51.38 * (U ** 4))
    epsilon0 = epsilon0 - (249.67 * (U ** 5))
    epsilon0 = epsilon0 - (39.05 * (U ** 6))
    epsilon0 = epsilon0 + (7.12 * (U ** 7))
    epsilon0 = epsilon0 + (27.87 * (U ** 8))
    epsilon0 = epsilon0 + (5.79 * (U ** 9))
    epsilon0 = epsilon0 + (2.45 * (U ** 10))

    epsilon0 = epsilon0 + delta_epsilon
    epsilon0 = mod(epsilon0, 26.0) ! if the obliquity of the ecliptic is calculated more than 10,000 years from J2000.0. the numbers become silly. This keeps it within sensible boundaries.
    ! nut = (epsilon0, delta_psi)
    nut(1) = epsilon0
    nut(2) = delta_psi
  end subroutine nutation

  subroutine getsid(jday, midnight)
    ! Calculate sidereal time at Greenwich
    ! Based on Meeus, chapter 12
    real(8), intent(in) :: jday ! Julian Day in question; must end in 0.5 because we're interested in midnight
    !real(8), intent(in) :: inst ! time since midnight that we're interested in
    !real(8), intent(out), dimension(2) :: sid ! Sidereal time at midnight and at the desired moment
    real(8), intent(out) :: midnight

    real(8) :: T ! Julian centuries since J2000.0
    ! real(8) :: midnight ! Sidereal time at midnight
    ! real(8) :: alpha
    real(8) :: corr ! correction to mean sidereal time to get apparent sidereal time
    real(8), dimension(2) :: epsi ! nutation factors
    real(8) :: d2r ! convert degrees to radians

    T = (jday - 2451545.0) / 36525.0
    midnight = 100.46061837 + (36000.770053608 * T) + (0.000387933 * T * T) - ((T ** 3) / 38710000.0)
    !alpha = inst * 1.00273790935
    !theta = midnight + alpha

    d2r = 4.0 * atan(1.0) / 180.0
    call nutation(jday, epsi)
    corr = (epsi(2) * cos(d2r * epsi(1))) / 15.0
    midnight = midnight + corr
    !theta = theta + corr
    !sid = (midnight, theta)
  end subroutine getsid

  subroutine sidstant(jday, inst, answer)
    ! Calculate sidereal time for any instant a Greenwich
    ! Based on Meeus, chapter 12
    real(8), intent(in) :: jday ! Julian Day in question
    real(8), intent(in) :: inst ! time since midnight that we're interested in
    real(8), intent(out) :: answer ! sidereal time at inst

    real(8) :: midnight ! sidereal time at midnight
    real(8) :: inc ! amount to add to the time at midnight

    call getsid(jday, midnight)
    inc = inst * 1.00273790935
    answer = midnight + inc
  end subroutine sidstant

  subroutine riset(jday, lon, lat, deltat, ra2000, dec2000, distance, rv, deltara, deltadec, id, time)
    ! Calculate time of rising and setting
    ! Based on Meeus, chapter 15

    ! This ignores the effect of atmospheric refraction because it's very very very small and too unpredicatable.
    
    real(8), intent(in) :: jday ! Julian Day in question
    real(8), intent(in) :: lon ! observer's longitude, in degrees
    real(8), intent(in) :: lat ! observer's latitude, in degrees
    real(8), intent(in) :: deltat ! difference between universal time and dynamical time
    !real(8), dimension(2), intent(in) :: yesterday ! RA and Dec for the previous day. This algorithm assumes they are in radians
    !real(8), dimension(2), intent(in) :: today ! RA and Dec for day in question. This algorithm assumes they are in radians
    !real(8), dimension(2), intent(in) :: tomorrow ! RA and Dec for next day. This algorithm assumes they are in radians
    real(8), intent(in) :: ra2000 ! right ascension at J2000.0, in degrees. This has to be looked up.
    real(8), intent(in) :: dec2000 ! declination at J2000.0, in degrees. This has to be looked up.
    real(8), intent(in) :: distance ! distance from the sun, in parsecs. This has to be looked up.
    real(8), intent(in) :: rv ! radial velocity, in parsecs per year
    real(8), intent(in) :: deltara ! RA component of proper motion, in ARCSECONDS. This has to be looked up.
    real(8), intent(in) :: deltadec ! Dec component of proper motion, in ARCSECONDS. This has to be looked up.
    integer, intent(in) :: id !What is actually rising or setting?
    real(8), dimension(2), intent(out) :: time ! time of sunrise and sunset, in days and fractions of a day
    
    real(8) :: pi
    real(8) :: d2r ! convert degrees to radians
    real(8) :: r2d ! convert radians to degrees

    real(8) :: h0 ! standard altitude, in degrees    
    real(8) :: testval ! initial check
    ! real(8) :: approx ! approximate time related to sunset

    real(8) :: sid ! Sidereal time at midnight on the day in question
    real(8) :: bigh
    
    real(8) :: nr ! used in calculating a modification to the rise time
    real(8) :: ns ! used in calculating a modification to the set time
    real(8) :: a_ra ! used in calculating interpolation
    real(8) :: b_ra ! used in calculating interpolation
    real(8) :: c_ra  ! used in calculating interpolation
    real(8) :: a_dec ! used in calculating interpolation
    real(8) :: b_dec ! used in calculating interpolation
    real(8) :: c_dec ! used in calculating interpolation
    real(8) :: rai_r ! right ascension, interpolated, for sunrise
    real(8) :: deci_r ! declination, interpolated, for sunrise
    real(8) :: rai_s ! right ascension, interpolated, for sunset
    real(8) :: deci_s ! declination, interpolated, for sunset
    ! real(8), dimension(2) :: inr ! interpolated RA and dec for sunrise
    ! real(8), dimension(2) :: ins ! interpolated RA and dec for sunset
    real(8) :: alt_r ! altitude at sunrise
    real(8) :: alt_s ! altitude at sunset
    real(8) :: ha_r ! hour angle of rising sun, in degrees
    real(8) :: ha_s ! hour angle of setting sun, in degrees
    
    real(8) :: transit ! time the sun crosses the meridian
    real(8) :: theta_r ! sidereal time of the sunrise converted into degrees
    real(8) :: theta_s ! sidereal time of the sunset converted into degrees

    real(8) :: delta_r ! modification to get true rising time
    real(8) :: delta_s ! modification to get true setting time

    real(8), dimension(2) :: yesterday ! RA and dec of prev day
    real(8), dimension(2) :: today ! RA and dec of day
    real(8), dimension(2) :: tomorrow ! RA and dec of next day
    
    pi = 4.0 * atan(1.0)
    d2r = pi / 180.0
    r2d = 180.0 / pi

    !print *, jday, lon, lat, deltat, ra2000, dec2000, distance, rv, deltara, deltadec, id

    if (id == 0) then
       ! It's the sun
       h0 = -0.8333 ! this is in degrees, not radians
    else
       ! It's a star or a planet
       h0 = -0.5667 ! this is in degrees, not radians
    end if
    !print *, h0
    !print *, jday
    !print *, ra2000
    !print *, dec2000
    !print *, distance
    !print *, rv
    !print *, deltara
    !print *, deltadec
    
    call precession((jday - 1), ra2000, dec2000, distance, rv, deltara, deltadec, yesterday) ! right ascension and declination of the previous day, in degrees
    call precession(jday, ra2000, dec2000, distance, rv, deltara, deltadec, today) ! right ascension and declination of the day in question, in degrees
    call precession((jday + 1), ra2000, dec2000, distance, rv, deltara, deltadec, tomorrow) ! right ascension and declination of the next day, in degrees
    testval = (sin(d2r * h0) - (sin(d2r * lat) * sin(d2r * tomorrow(2)))) / (cos(d2r * lat) * cos(d2r * tomorrow(2)))
    !print *, yesterday
    !print *, today
    !print *, tomorrow
    !print *, testval
    !print *, yesterday
    !print *, today
    !print *, tomorrow
    !print *, sin(d2r * h0)
    !print *, sin(d2r * lat)
    !print *, sin(d2r * tomorrow(2))
    !print *, (sin(d2r * h0) - (sin(d2r * lat) * sin(d2r * tomorrow(2))))
    !print *, cos(d2r * lat)
    !print *, cos(d2r * tomorrow(2))
    !print *, cos(d2r * lat) * cos(d2r * tomorrow(2))
    !print *, testval
    if (abs(testval) > 1.0) then
       time = 25.0
    else
       bigh = acos(testval)
       call getsid(jday, sid)
       
       transit = (tomorrow(1) + lon - sid) / 360.0
       time(1) = transit - ((r2d * bigh) / 360.0) ! sunrise
       time(2) = transit + ((r2d * bigh) / 360.0) ! sunset

       do while (time(1) < 0)
          time(1) = time(1) + 1.0
       end do
       do while (time(2) < 0)
          time(2) = time(2) + 1.0
       end do

       time(1) = mod(time(1), 1.0)
       time(2) = mod(time(2), 1.0)

       !call getsid(jday, sid)

       theta_r = sid + (360.985647 * time(1))
       theta_s = sid + (360.985647 * time(2))

       nr = time(1) + (deltat / 86400.0)
       ns = time(2) + (deltat / 86400.0)

       ! interpolation

       ! RA interpolation
       a_ra = today(1) - yesterday(1)
       b_ra = tomorrow(1) - today(1)
       c_ra = b_ra - a_ra
       rai_r = (nr / 2.0) * (a_ra + b_ra + (nr * c_ra))
       rai_s = (ns / 2.0) * (a_ra + b_ra + (ns * c_ra))

       ! dec interpolation
       a_dec = today(2) - yesterday(2)
       b_dec = tomorrow(2) - today(2)
       c_dec = b_dec - a_dec
       deci_r = (nr / 2.0) * (a_dec + b_dec + (nr * c_dec))
       deci_s = (ns / 2.0) * (a_dec + b_dec + (ns * c_dec))

       ha_r = theta_r - lon - rai_r
       ha_s = theta_s - lon - rai_s       

       alt_r = r2d * asin((sin(d2r * lat) * sin(d2r * deci_r)) + (cos(d2r * lat) * cos(ha_r)))
       alt_s = r2d * asin((sin(d2r * lat) * sin(d2r * deci_r)) + (cos(d2r * lat) * cos(ha_r)))

       ! MAKE SURE alt IS IN DEGREES. IF NOT, BE SURE TO CORRECT IT.

       delta_r = (alt_r - h0) / (360.0 * cos(deci_r * d2r) * cos(lat * d2r) * sin(bigh))
       delta_s = (alt_s - h0) / (360.0 * cos(deci_s * d2r) * cos(lat * d2r) * sin(bigh))

       time(1) = time(1) + delta_r
       time(2) = time(2) + delta_s

       ! time(x) now gives the fraction of a day that has elapsed since midnight when the body rises (x == 1) or sets (x == 2)
       ! these can be compared to floats, Fractions, and Decimals in Python, but might need to be converted into other units
       ! multiply time(x) by 24 to get the time in hours (which will still have a decimal portion; for example, 06:30 would register as 6.5)
       ! multiply time(x) by 1440 to get the time in minutes, or by 86400 to get the time in seconds.
       ! time(x) can even be multiplied by 25920 to get the time in chalakim.
    end if
    !print *, time
  end subroutine riset
end module sidereal
