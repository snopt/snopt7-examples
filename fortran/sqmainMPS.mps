NAME          HS 118
ROWS
 N  COST
 G  LINCON1
 G  LINCON2
 G  LINCON3
 G  LINCON4
 G  LINCON5
 G  LINCON6
 G  LINCON7
 G  LINCON8
 G  LINCON9
 G  LINCON10
 G  LINCON11
 G  LINCON12
 G  LINCON13
 G  LINCON14
 G  LINCON15
 G  LINCON16
 G  LINCON17
COLUMNS
    X1        COST             2.3
    X1        LINCON1          -1.
    X1        LINCON13         1.0
    X2        LINCON5          -1.
    X2        LINCON13         1.0
    X2        COST             1.7
    X3        LINCON9          -1.
    X3        LINCON13         1.0
    X3        COST             2.2
    X4        LINCON1          1.0
    X4        LINCON2          -1.
    X4        LINCON14         1.0
    X4        COST             2.3
    X5        LINCON5          1.0
    X5        LINCON6          -1.
    X5        LINCON14         1.0
    X5        COST             1.7
    X6        LINCON9          1.0
    X6        LINCON10         -1.
    X6        LINCON14         1.0
    X6        COST             2.2
    X7        LINCON2          1.0
    X7        LINCON3          -1.
    X7        LINCON15         1.0
    X7        COST             2.3
    X8        LINCON6          1.0
    X8        LINCON7          -1.
    X8        LINCON15         1.0
    X8        COST             1.7
    X9        LINCON10         1.0
    X9        LINCON11         -1.
    X9        LINCON15         1.0
    X9        COST             2.2
    X10       LINCON3          1.0
    X10       LINCON4          -1.
    X10       LINCON16         1.0
    X10       COST             2.3
    X11       LINCON7          1.0
    X11       LINCON8          -1.
    X11       LINCON16         1.0
    X11       COST             1.7
    X12       LINCON11         1.0
    X12       LINCON12         -1.
    X12       LINCON16         1.0
    X12       COST             2.2
    X13       LINCON4          1.0
    X13       LINCON17         1.0
    X13       COST             2.3
    X14       LINCON8          1.0
    X14       LINCON17         1.0
    X14       COST             1.7
    X15       LINCON12         1.0
    X15       LINCON17         1.0
    X15       COST             2.2
RHS
    B         LINCON1        -7.0      LINCON2       -7.0
    B         LINCON3        -7.0      LINCON4       -7.0
    B         LINCON5        -7.0      LINCON6       -7.0
    B         LINCON7        -7.0      LINCON8       -7.0
    B         LINCON9        -7.0      LINCON10      -7.0
    B         LINCON11       -7.0      LINCON12      -7.0
    B         LINCON13       60.0      LINCON14      50.0
    B         LINCON15       70.0      LINCON16      85.0
    B         LINCON17      100.0
RANGES
    RANGEC9   LINCON1        13.0      LINCON2       13.0
    RANGEC9   LINCON3        13.0      LINCON4       13.0
    RANGEC9   LINCON5        14.0      LINCON6       14.0
    RANGEC9   LINCON7        14.0      LINCON8       14.0
    RANGEC9   LINCON9        13.0      LINCON10      13.0
    RANGEC9   LINCON11       13.0      LINCON12      13.0
BOUNDS
 LO BOUND1    X1              8.0
 UP BOUND1    X1             21.0
 LO BOUND1    X2             43.0
 UP BOUND1    X2             57.0
 LO BOUND1    X3              3.0
 UP BOUND1    X3             16.0
 UP BOUND1    X4             90.0
 UP BOUND1    X5            120.0
 UP BOUND1    X6             60.0
 UP BOUND1    X7             90.0
 UP BOUND1    X8            120.0
 UP BOUND1    X9             60.0
 UP BOUND1    X10            90.0
 UP BOUND1    X11           120.0
 UP BOUND1    X12            60.0
 UP BOUND1    X13            90.0
 UP BOUND1    X14           120.0
 UP BOUND1    X15            60.0
 FR INITIAL   X1             20.0
 FR INITIAL   X2             55.0
 FR INITIAL   X3             15.0
 FR INITIAL   X4             20.0
 FR INITIAL   X5             60.0
 FR INITIAL   X6             20.0
 FR INITIAL   X7             20.0
 FR INITIAL   X8             60.0
 FR INITIAL   X9             20.0
 FR INITIAL   X10            20.0
 FR INITIAL   X11            60.0
 FR INITIAL   X12            20.0
 FR INITIAL   X13            20.0
 FR INITIAL   X14            60.0
 FR INITIAL   X15            20.0
ENDATA
ENDRUN
