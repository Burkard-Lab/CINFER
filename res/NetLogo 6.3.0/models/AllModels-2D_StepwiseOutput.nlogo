breed [tumorcells tumorcell]

globals
[
  tumorcell-count
  Chrom1-final
  Chrom2-final
  Chrom3-final
  Chrom4-final
  Chrom5-final
  Chrom6-final
  Chrom7-final
  Chrom8-final
  Chrom9-final
  Chrom10-final
  Chrom11-final
  Chrom12-final
  Chrom13-final
  Chrom14-final
  Chrom15-final
  Chrom16-final
  Chrom17-final
  Chrom18-final
  Chrom19-final
  Chrom20-final
  Chrom21-final
  Chrom22-final
  ChromX-final
  fitness-final
  x
  CIN-rate
  total-misseg
  MKV
  tandem-losses
]

tumorcells-own
[
  MChrom1 ;--------------------------------- Maternal chromosome attribute
  MChrom2
  MChrom3
  MChrom4
  MChrom5
  MChrom6
  MChrom7
  MChrom8
  MChrom9
  MChrom10
  MChrom11
  MChrom12
  MChrom13
  MChrom14
  MChrom15
  MChrom16
  MChrom17
  MChrom18
  MChrom19
  MChrom20
  MChrom21
  MChrom22
  MChromX ;----------------------------------------------------------

  PChrom1 ;--------------------------------- Paternal chromosome attribute
  PChrom2
  PChrom3
  PChrom4
  PChrom5
  PChrom6
  PChrom7
  PChrom8
  PChrom9
  PChrom10
  PChrom11
  PChrom12
  PChrom13
  PChrom14
  PChrom15
  PChrom16
  PChrom17
  PChrom18
  PChrom19
  PChrom20
  PChrom21
  PChrom22
  PChromX;----------------------------------------------------------

  Chrom1 ;--------------------------------- Maternal + Paternal chromosome attribute
  Chrom2
  Chrom3
  Chrom4
  Chrom5
  Chrom6
  Chrom7
  Chrom8
  Chrom9
  Chrom10
  Chrom11
  Chrom12
  Chrom13
  Chrom14
  Chrom15
  Chrom16
  Chrom17
  Chrom18
  Chrom19
  Chrom20
  Chrom21
  Chrom22
  ChromX;----------------------------------------------------------

  MP1 ;------------------ Probability that a maternal chromosome is missegregated
  MP2
  MP3
  MP4
  MP5
  MP6
  MP7
  MP8
  MP9
  MP10
  MP11
  MP12
  MP13
  MP14
  MP15
  MP16
  MP17
  MP18
  MP19
  MP20
  MP21
  MP22
  MPX ;---------------------------------------------------


  PP1 ;------------------ Probability that a paternal chromosome is missegregated
  PP2
  PP3
  PP4
  PP5
  PP6
  PP7
  PP8
  PP9
  PP10
  PP11
  PP12
  PP13
  PP14
  PP15
  PP16
  PP17
  PP18
  PP19
  PP20
  PP21
  PP22
  PPX;---------------------------------------------------


  MGL1 ;-------------------------------------------------- Maternal gain/loss action attribute
  MGL2
  MGL3
  MGL4
  MGL5
  MGL6
  MGL7
  MGL8
  MGL9
  MGL10
  MGL11
  MGL12
  MGL13
  MGL14
  MGL15
  MGL16
  MGL17
  MGL18
  MGL19
  MGL20
  MGL21
  MGL22
  MGLX ;---------------------------------------------------

  PGL1 ;-------------------------------------------------- Paternal gain/loss action attribute
  PGL2
  PGL3
  PGL4
  PGL5
  PGL6
  PGL7
  PGL8
  PGL9
  PGL10
  PGL11
  PGL12
  PGL13
  PGL14
  PGL15
  PGL16
  PGL17
  PGL18
  PGL19
  PGL20
  PGL21
  PGL22
  PGLX;---------------------------------

  mmiseg1
  mmiseg2
  mmiseg3
  mmiseg4
  mmiseg5
  mmiseg6
  mmiseg7
  mmiseg8
  mmiseg9
  mmiseg10
  mmiseg11
  mmiseg12
  mmiseg13
  mmiseg14
  mmiseg15
  mmiseg16
  mmiseg17
  mmiseg18
  mmiseg19
  mmiseg20
  mmiseg21
  mmiseg22
  mmisegX

  pmiseg1
  pmiseg2
  pmiseg3
  pmiseg4
  pmiseg5
  pmiseg6
  pmiseg7
  pmiseg8
  pmiseg9
  pmiseg10
  pmiseg11
  pmiseg12
  pmiseg13
  pmiseg14
  pmiseg15
  pmiseg16
  pmiseg17
  pmiseg18
  pmiseg19
  pmiseg20
  pmiseg21
  pmiseg22
  pmisegX

  Fc1 ;--------------------------------------------------- normalized gene values for each chromosome
  Fc2
  Fc3
  Fc4
  Fc5
  Fc6
  Fc7
  Fc8
  Fc9
  Fc10
  Fc11
  Fc12
  Fc13
  Fc14
  Fc15
  Fc16
  Fc17
  Fc18
  Fc19
  Fc20
  Fc21
  Fc22
  Fcx ;-------------------------------------------end
  fitd1 ;------------------------discrete fitness for each chromosome for each step
  fitd2
  fitd3
  fitd4
  fitd5
  fitd6
  fitd7
  fitd8
  fitd9
  fitd10
  fitd11
  fitd12
  fitd13
  fitd14
  fitd15
  fitd16
  fitd17
  fitd18
  fitd19
  fitd20
  fitd21
  fitd22
  fitdx ;-------------------------------------end
  CS1 ;---------------chrom scores (Davoli et al. 13)
  CS2
  CS3
  CS4
  CS5
  CS6
  CS7
  CS8
  CS9
  CS10
  CS11
  CS12
  CS13
  CS14
  CS15
  CS16
  CS17
  CS18
  CS19
  CS20
  CS21
  CS22;--------end

  CINList

  w;--- weight for chrom scores


  totalmisseg
  chromavg
  karyotype
  CIN-level
  MitosisProbability
  fitness
]

to setup
  clear-all
  set-default-shape tumorcells "circle"
  set-periodicity
  setup-progenitor
  reset-ticks
  setup-progenitor-attributes
  translate-cin-rates
end

to translate-cin-rates
  ifelse initial-CIN = 0
  [set CIN-rate 0]
  [set CIN-rate 100 / initial-CIN]

end

to setup-progenitor
  create-tumorcells initial-number
  ask tumorcells
  [
    set MitosisProbability random 100
    set karyotype (initial-ploidy) * 23
    assign-gain-loss

    set MChrom1 initial-ploidy / 2
    set MChrom2 initial-ploidy / 2
    set MChrom3 initial-ploidy / 2
    set MChrom4 initial-ploidy / 2
    set MChrom5 initial-ploidy / 2
    set MChrom6 initial-ploidy / 2
    set MChrom7 initial-ploidy / 2
    set MChrom8 initial-ploidy / 2
    set MChrom9 initial-ploidy / 2
    set MChrom10 initial-ploidy / 2
    set MChrom11 initial-ploidy / 2
    set MChrom12 initial-ploidy / 2
    set MChrom13 initial-ploidy / 2
    set MChrom14 initial-ploidy / 2
    set MChrom15 initial-ploidy / 2
    set MChrom16 initial-ploidy / 2
    set MChrom17 initial-ploidy / 2
    set MChrom18 initial-ploidy / 2
    set MChrom19 initial-ploidy / 2
    set MChrom20 initial-ploidy / 2
    set MChrom21 initial-ploidy / 2
    set MChrom22 initial-ploidy / 2
    set MChromX initial-ploidy  / 2

    set PChrom1 initial-ploidy / 2
    set PChrom2 initial-ploidy / 2
    set PChrom3 initial-ploidy / 2
    set PChrom4 initial-ploidy / 2
    set PChrom5 initial-ploidy / 2
    set PChrom6 initial-ploidy / 2
    set PChrom7 initial-ploidy / 2
    set PChrom8 initial-ploidy / 2
    set PChrom9 initial-ploidy / 2
    set PChrom10 initial-ploidy / 2
    set PChrom11 initial-ploidy / 2
    set PChrom12 initial-ploidy / 2
    set PChrom13 initial-ploidy / 2
    set PChrom14 initial-ploidy / 2
    set PChrom15 initial-ploidy / 2
    set PChrom16 initial-ploidy / 2
    set PChrom17 initial-ploidy / 2
    set PChrom18 initial-ploidy / 2
    set PChrom19 initial-ploidy / 2
    set PChrom20 initial-ploidy / 2
    set PChrom21 initial-ploidy / 2
    set PChrom22 initial-ploidy / 2
    set PChromX initial-ploidy / 2

    set Chrom1 (MChrom1 + PChrom1)
    set Chrom2 (MChrom2 + PChrom2)
    set Chrom3 (MChrom3 + PChrom3)
    set Chrom4 (MChrom4 + PChrom4)
    set Chrom5 (MChrom5 + PChrom5)
    set Chrom6 (MChrom6 + PChrom6)
    set Chrom7 (MChrom7 + PChrom7)
    set Chrom8 (MChrom8 + PChrom8)
    set Chrom9 (MChrom9 + PChrom9)
    set Chrom10 (MChrom10 + PChrom10)
    set Chrom11 (MChrom11 + PChrom11)
    set Chrom12 (MChrom12 + PChrom12)
    set Chrom13 (MChrom13 + PChrom13)
    set Chrom14 (MChrom14 + PChrom14)
    set Chrom15 (MChrom15 + PChrom15)
    set Chrom16 (MChrom16 + PChrom16)
    set Chrom17 (MChrom17 + PChrom17)
    set Chrom18 (MChrom18 + PChrom18)
    set Chrom19 (MChrom19 + PChrom19)
    set Chrom20 (MChrom20 + PChrom20)
    set Chrom21 (MChrom21 + PChrom21)
    set Chrom22 (MChrom22 + PChrom22)
    set ChromX (MChromX + PChromX)

    set Fc1 0.095
    set Fc2 0.072
    set Fc3 0.056
    set Fc4 0.046
    set Fc5 0.048
    set Fc6 0.056
    set Fc7 0.052
    set Fc8 0.04
    set Fc9 0.042
    set Fc10 0.041
    set Fc11 0.055
    set Fc12 0.047
    set Fc13 0.026
    set Fc14 0.039
    set Fc15 0.034
    set Fc16 0.036
    set Fc17 0.046
    set Fc18 0.018
    set Fc19 0.047
    set Fc20 0.025
    set Fc21 0.015
    set Fc22 0.022
    set Fcx 0.041

    set CS1 0.0261
    set CS2 0.0661
    set CS3 0.0822
    set CS4 0.0470
    set CS5 0.0434
    set CS6 0.0277
    set CS7 0.1468
    set CS8 0.0872
    set CS9 0.0630
    set CS10 0.0016
    set CS11 0.0523
    set CS12 0.1128
    set CS13 -0.0191
    set CS14 0.0475
    set CS15 0.0382
    set CS16 0.0370
    set CS17 0.0528
    set CS18 -0.0344
    set CS19 0.0490
    set CS20 0.0888
    set CS21 -0.0062
    set CS22 -0.0100

        set MP1 999
    set MP2 999
    set MP3 999
    set MP4 999
    set MP5 999
    set MP6 999
    set MP7 999
    set MP8 999
    set MP9 999
    set MP10 999
    set MP11 999
    set MP12 999
    set MP13 999
    set MP14 999
    set MP15 999
    set MP16 999
    set MP17 999
    set MP18 999
    set MP19 999
    set MP20 999
    set MP21 999
    set MP22 999
    set MPX 999

    set PP1 999
    set PP2 999
    set PP3 999
    set PP4 999
    set PP5 999
    set PP6 999
    set PP7 999
    set PP8 999
    set PP9 999
    set PP10 999
    set PP11 999
    set PP12 999
    set PP13 999
    set PP14 999
    set PP15 999
    set PP16 999
    set PP17 999
    set PP18 999
    set PP19 999
    set PP20 999
    set PP21 999
    set PP22 999
    set PPX 999

    set CINList (list MP1 MP2 MP3 MP4 MP5 MP6 MP7 MP8 MP9 MP10 MP11 MP12 MP13 MP14 MP15 MP16 MP17 MP18 MP19 MP20 MP21 MP22 MPX PP1 PP2 PP3 PP4 PP5 PP6 PP7 PP8 PP9 PP10 PP11 PP12 PP13 PP14 PP15 PP16 PP17 PP18 PP19 PP20 PP21 PP22 PPX)
  ]
end

to setup-progenitor-attributes
  ask tumorcells [
  calculate-chromosome-scores
  calculate-fitness
  assign-gain-loss
  set-color
  ]
  assign-cin
end

to assign-gain-loss
ask tumorcells[
    set MGL1 random 2
    set MGL2 random 2
    set MGL3 random 2
    set MGL4 random 2
    set MGL5 random 2
    set MGL6 random 2
    set MGL7 random 2
    set MGL8 random 2
    set MGL9 random 2
    set MGL10 random 2
    set MGL11 random 2
    set MGL12 random 2
    set MGL13 random 2
    set MGL14 random 2
    set MGL15 random 2
    set MGL16 random 2
    set MGL17 random 2
    set MGL18 random 2
    set MGL19 random 2
    set MGL20 random 2
    set MGL21 random 2
    set MGL22 random 2
    set MGLX random 2

    set PGL1 random 2
    set PGL2 random 2
    set PGL3 random 2
    set PGL4 random 2
    set PGL5 random 2
    set PGL6 random 2
    set PGL7 random 2
    set PGL8 random 2
    set PGL9 random 2
    set PGL10 random 2
    set PGL11 random 2
    set PGL12 random 2
    set PGL13 random 2
    set PGL14 random 2
    set PGL15 random 2
    set PGL16 random 2
    set PGL17 random 2
    set PGL18 random 2
    set PGL19 random 2
    set PGL20 random 2
    set PGL21 random 2
    set PGL22 random 2
    set PGLX random 2
  ]
end

to set-color
    if fitness >= 0.9 [set color 19.9]
    if fitness < 0.9 [set color 19.9 - 0.5]
    if fitness < 0.85 [set color 19.9 - 1.0]
    if fitness < 0.8 [set color 19.9 - 1.5]
    if fitness < 0.75 [set color 19.9 - 2.5]
    if fitness < 0.7 [set color 19.9 - 2.5]
    if fitness < 0.65 [set color 19.9 - 3.0]
    if fitness < 0.6 [set color 19.9 - 3.5]
    if fitness < 0.55 [set color 19.9 - 4.0]
    if fitness < 0.5 [set color 19.9 - 4.5]
  set color lput 100 extract-rgb color

end

to go
if (ticks < 1) [ reset-timer ]
if ticks = end-ticks or count tumorcells >= end-population
  [
    stop
  ]
  calculate-chromosome-scores
  calculate-fitness
  check-death
  assign-mitosis-probability
  assign-cin
  assign-gain-loss
  reset-division-state

  ask tumorcells[
    if (100 * fitness) >= MitosisProbability
    [
      mitosis1
      mitosis2
  ]]
  count-tandem-losses
  output-chromosome-state
tick
end

to count-tandem-losses
  set tandem-losses count tumorcells with [
    Chrom1 = 0 or
    Chrom2 = 0 or
    Chrom3 = 0 or
    Chrom4 = 0 or
    Chrom5 = 0 or
    Chrom6 = 0 or
    Chrom7 = 0 or
    Chrom8 = 0 or
    Chrom9 = 0 or
    Chrom10 = 0 or
    Chrom11 = 0 or
    Chrom12 = 0 or
    Chrom13 = 0 or
    Chrom14 = 0 or
    Chrom15 = 0 or
    Chrom16 = 0 or
    Chrom17 = 0 or
    Chrom18 = 0 or
    Chrom19 = 0 or
    Chrom20 = 0 or
    Chrom21 = 0 or
    Chrom22 = 0 or
    ChromX = 0 ]
end

to assign-cin
if (ticks >= cin-period-start) [
    ask tumorcells [
      ifelse initial-CIN = 0
    [
    set MP1 999
    set MP2 999
    set MP3 999
    set MP4 999
    set MP5 999
    set MP6 999
    set MP7 999
    set MP8 999
    set MP9 999
    set MP10 999
    set MP11 999
    set MP12 999
    set MP13 999
    set MP14 999
    set MP15 999
    set MP16 999
    set MP17 999
    set MP18 999
    set MP19 999
    set MP20 999
    set MP21 999
    set MP22 999
    set MPX 999

    set PP1 999
    set PP2 999
    set PP3 999
    set PP4 999
    set PP5 999
    set PP6 999
    set PP7 999
    set PP8 999
    set PP9 999
    set PP10 999
    set PP11 999
    set PP12 999
    set PP13 999
    set PP14 999
    set PP15 999
    set PP16 999
    set PP17 999
    set PP18 999
    set PP19 999
    set PP20 999
    set PP21 999
    set PP22 999
    set PPX 999
    ]
      [
    set MP1 random initial-CIN
    set MP2 random initial-CIN
    set MP3 random initial-CIN
    set MP4 random initial-CIN
    set MP5 random initial-CIN
    set MP6 random initial-CIN
    set MP7 random initial-CIN
    set MP8 random initial-CIN
    set MP9 random initial-CIN
    set MP10 random initial-CIN
    set MP11 random initial-CIN
    set MP12 random initial-CIN
    set MP13 random initial-CIN
    set MP14 random initial-CIN
    set MP15 random initial-CIN
    set MP16 random initial-CIN
    set MP17 random initial-CIN
    set MP18 random initial-CIN
    set MP19 random initial-CIN
    set MP20 random initial-CIN
    set MP21 random initial-CIN
    set MP22 random initial-CIN
    set MPX random initial-CIN

    set PP1 random initial-CIN
    set PP2 random initial-CIN
    set PP3 random initial-CIN
    set PP4 random initial-CIN
    set PP5 random initial-CIN
    set PP6 random initial-CIN
    set PP7 random initial-CIN
    set PP8 random initial-CIN
    set PP9 random initial-CIN
    set PP10 random initial-CIN
    set PP11 random initial-CIN
    set PP12 random initial-CIN
    set PP13 random initial-CIN
    set PP14 random initial-CIN
    set PP15 random initial-CIN
    set PP16 random initial-CIN
    set PP17 random initial-CIN
    set PP18 random initial-CIN
    set PP19 random initial-CIN
    set PP20 random initial-CIN
    set PP21 random initial-CIN
    set PP22 random initial-CIN
    set PPX random initial-CIN
  ]
    ]
  ]
if (ticks > cin-period-end) [
      ask tumorcells [
        set MP1 999
    set MP2 999
    set MP3 999
    set MP4 999
    set MP5 999
    set MP6 999
    set MP7 999
    set MP8 999
    set MP9 999
    set MP10 999
    set MP11 999
    set MP12 999
    set MP13 999
    set MP14 999
    set MP15 999
    set MP16 999
    set MP17 999
    set MP18 999
    set MP19 999
    set MP20 999
    set MP21 999
    set MP22 999
    set MPX 999

    set PP1 999
    set PP2 999
    set PP3 999
    set PP4 999
    set PP5 999
    set PP6 999
    set PP7 999
    set PP8 999
    set PP9 999
    set PP10 999
    set PP11 999
    set PP12 999
    set PP13 999
    set PP14 999
    set PP15 999
    set PP16 999
    set PP17 999
    set PP18 999
    set PP19 999
    set PP20 999
    set PP21 999
    set PP22 999
    set PPX 999
    ]
  ]
end

to set-periodicity
    if (periodicity = "constant")[
    set cin-period-start 0
    set cin-period-end end-ticks
  ]
      if (periodicity = "early")[
    set cin-period-start 0
    set cin-period-end round(end-ticks * 0.2)
  ]
      if (periodicity = "middle")[
    set cin-period-start round(end-ticks * 0.4)
    set cin-period-end round(end-ticks * 0.6)
  ]
      if (periodicity = "late")[
    set cin-period-start round(end-ticks * 0.8)
    set cin-period-end end-ticks
  ]
end

to assign-mitosis-probability
  ask tumorcells [
    set MitosisProbability random mitosis-probability
  ]
end

to calculate-missegregations
set totalmisseg (abs(mchrom1 - 1) +
    abs(mchrom2 - 1) +
    abs(mchrom3 - 1) +
    abs(mchrom4 - 1) +
    abs(mchrom5 - 1) +
    abs(mchrom6 - 1) +
    abs(mchrom7 - 1) +
    abs(mchrom8 - 1) +
    abs(mchrom9 - 1) +
    abs(mchrom10 - 1) +
    abs(mchrom11 - 1) +
    abs(mchrom12 - 1) +
    abs(mchrom13 - 1) +
    abs(mchrom14 - 1) +
    abs(mchrom15 - 1) +
    abs(mchrom16 - 1) +
    abs(mchrom17 - 1) +
    abs(mchrom18 - 1) +
    abs(mchrom19 - 1) +
    abs(mchrom20 - 1) +
    abs(mchrom21 - 1) +
    abs(mchrom22 - 1) +
    abs(mchromx - 1) +
    abs(pchrom1 - 1) +
    abs(pchrom2 - 1) +
    abs(pchrom3 - 1) +
    abs(pchrom4 - 1) +
    abs(pchrom5 - 1) +
    abs(pchrom6 - 1) +
    abs(pchrom7 - 1) +
    abs(pchrom8 - 1) +
    abs(pchrom9 - 1) +
    abs(pchrom10 - 1) +
    abs(pchrom11 - 1) +
    abs(pchrom12 - 1) +
    abs(pchrom13 - 1) +
    abs(pchrom14 - 1) +
    abs(pchrom15 - 1) +
    abs(pchrom16 - 1) +
    abs(pchrom17 - 1) +
    abs(pchrom18 - 1) +
    abs(pchrom19 - 1) +
    abs(pchrom20 - 1) +
    abs(pchrom21 - 1) +
    abs(pchrom22 - 1) +
    abs(pchromx - 1))
end

to reset-division-state
  ask tumorcells [
  set mmiseg1 0
  set mmiseg2 0
  set mmiseg3 0
  set mmiseg4 0
  set mmiseg5 0
  set mmiseg6 0
  set mmiseg7 0
  set mmiseg8 0
  set mmiseg9 0
  set mmiseg10 0
  set mmiseg11 0
  set mmiseg12 0
  set mmiseg13 0
  set mmiseg14 0
  set mmiseg15 0
  set mmiseg16 0
  set mmiseg17 0
  set mmiseg18 0
  set mmiseg19 0
  set mmiseg20 0
  set mmiseg21 0
  set mmiseg22 0
  set mmisegX 0

  set pmiseg1 0
  set pmiseg2 0
  set pmiseg3 0
  set pmiseg4 0
  set pmiseg5 0
  set pmiseg6 0
  set pmiseg7 0
  set pmiseg8 0
  set pmiseg9 0
  set pmiseg10 0
  set pmiseg11 0
  set pmiseg12 0
  set pmiseg13 0
  set pmiseg14 0
  set pmiseg15 0
  set pmiseg16 0
  set pmiseg17 0
  set pmiseg18 0
  set pmiseg19 0
  set pmiseg20 0
  set pmiseg21 0
  set pmiseg22 0
  set pmisegX 0
  ]
end

to mitosis1
   if MP1 < 100 and MChrom1 > 0
   [
    ifelse MGL1 = 0
    [set MChrom1 Mchrom1 - 1
     set karyotype karyotype - 1]
    [set MChrom1 MChrom1 + 1
     set karyotype karyotype + 1]
  set mmiseg1 1 ]

   if MP2 < 100 and MChrom2 > 0
   [
    ifelse MGL2 = 0
    [set MChrom2 MChrom2 - 1
     set karyotype karyotype - 1]
    [set MChrom2 MChrom2 + 1
     set karyotype karyotype + 1]
  set mmiseg2 1 ]

   if MP3 < 100 and MChrom3 > 0
   [
    ifelse MGL3 = 0
    [set MChrom3 MChrom3 - 1
     set karyotype karyotype - 1]
    [set MChrom3 MChrom3 + 1
     set karyotype karyotype + 1]
  set mmiseg3 1 ]

   if MP4 < 100 and MChrom4 > 0
   [
    ifelse MGL4 = 0
    [set MChrom4 MChrom4 - 1
     set karyotype karyotype - 1]
    [set MChrom4 MChrom4 + 1
     set karyotype karyotype + 1]
  set mmiseg4 1 ]

   if MP5 < 100 and MChrom5 > 0
   [
    ifelse MGL5 = 0
    [set MChrom5 MChrom5 - 1
     set karyotype karyotype - 1]
    [set MChrom5 MChrom5 + 1
     set karyotype karyotype + 1]
  set mmiseg5 1 ]

   if MP6 < 100 and MChrom6 > 0
   [
    ifelse MGL6 = 0
    [set MChrom6 MChrom6 - 1
     set karyotype karyotype - 1]
    [set MChrom6 MChrom6 + 1
     set karyotype karyotype + 1]
  set mmiseg6 1 ]

   if MP7 < 100 and MChrom7 > 0
   [
    ifelse MGL7 = 0
    [set MChrom7 MChrom7 - 1
     set karyotype karyotype - 1]
    [set MChrom7 MChrom7 + 1
     set karyotype karyotype + 1]
  set mmiseg7 1 ]

   if MP8 < 100 and MChrom8 > 0
   [
    ifelse MGL8 = 0
    [set MChrom8 MChrom8 - 1
     set karyotype karyotype - 1]
    [set MChrom8 MChrom8 + 1
     set karyotype karyotype + 1]
  set mmiseg8 1 ]

   if MP9 < 100 and MChrom9 > 0
   [
    ifelse MGL9 = 0
    [set MChrom9 MChrom9 - 1
     set karyotype karyotype - 1]
    [set MChrom9 MChrom9 + 1
     set karyotype karyotype + 1]
  set mmiseg9 1 ]

   if MP10 < 100 and MChrom10 > 0
   [
    ifelse MGL10 = 0
    [set MChrom10 MChrom10 - 1
     set karyotype karyotype - 1]
    [set MChrom10 MChrom10 + 1
     set karyotype karyotype + 1]
  set mmiseg10 1 ]

   if MP11 < 100 and MChrom11 > 0
   [
    ifelse MGL11 = 0
    [set MChrom11 MChrom11 - 1
     set karyotype karyotype - 1]
    [set MChrom11 MChrom11 + 1
     set karyotype karyotype + 1]
  set mmiseg11 1 ]

   if MP12 < 100 and MChrom12 > 0
   [
    ifelse MGL12 = 0
    [set MChrom12 MChrom12 - 1
     set karyotype karyotype - 1]
    [set MChrom12 MChrom12 + 1
     set karyotype karyotype + 1]
  set mmiseg12 1 ]

   if MP13 < 100 and MChrom13 > 0
   [
    ifelse MGL13 = 0
    [set MChrom13 MChrom13 - 1
     set karyotype karyotype - 1]
    [set MChrom13 MChrom13 + 1
     set karyotype karyotype + 1]
  set mmiseg13 1 ]

   if MP14 < 100 and MChrom14 > 0
   [
    ifelse MGL14 = 0
    [set MChrom14 MChrom14 - 1
     set karyotype karyotype - 1]
    [set MChrom14 MChrom14 + 1
     set karyotype karyotype + 1]
  set mmiseg14 1 ]

   if MP15 < 100 and MChrom15 > 0
   [
    ifelse MGL15 = 0
    [set MChrom15 MChrom15 - 1
     set karyotype karyotype - 1]
    [set MChrom15 MChrom15 + 1
     set karyotype karyotype + 1]
  set mmiseg15 1 ]

   if MP16 < 100 and MChrom16 > 0
   [
    ifelse MGL16 = 0
    [set MChrom16 MChrom16 - 1
     set karyotype karyotype - 1]
    [set MChrom16 MChrom16 + 1
     set karyotype karyotype + 1]
  set mmiseg16 1 ]

   if MP17 < 100 and MChrom17 > 0
   [
    ifelse MGL17 = 0
    [set MChrom17 MChrom17 - 1
     set karyotype karyotype - 1]
    [set MChrom17 MChrom17 + 1
     set karyotype karyotype + 1]
  set mmiseg17 1 ]

   if MP18 < 100 and MChrom18 > 0
   [
    ifelse MGL18 = 0
    [set MChrom18 MChrom18 - 1
     set karyotype karyotype - 1]
    [set MChrom18 MChrom18 + 1
     set karyotype karyotype + 1]
  set mmiseg18 1 ]

   if MP19 < 100 and MChrom19 > 0
   [
    ifelse MGL19 = 0
    [set MChrom19 MChrom19 - 1
     set karyotype karyotype - 1]
    [set MChrom19 MChrom19 + 1
     set karyotype karyotype + 1]
  set mmiseg19 1 ]

   if MP20 < 100 and MChrom20 > 0
   [
    ifelse MGL20 = 0
    [set MChrom20 MChrom20 - 1
     set karyotype karyotype - 1]
    [set MChrom20 MChrom20 + 1
     set karyotype karyotype + 1]
  set mmiseg20 1 ]

   if MP21 < 100 and MChrom21 > 0
   [
    ifelse MGL21 = 0
    [set MChrom21 MChrom21 - 1
     set karyotype karyotype - 1]
    [set MChrom21 MChrom21 + 1
     set karyotype karyotype + 1]
  set mmiseg21 1 ]

   if MP22 < 100 and MChrom22 > 0
   [
    ifelse MGL22 = 0
    [set MChrom22 MChrom22 - 1
     set karyotype karyotype - 1]
    [set MChrom22 MChrom22 + 1
     set karyotype karyotype + 1]
  set mmiseg22 1 ]

   if MPX < 100 and MChromX > 0
   [
    ifelse MGLX = 0
    [set MChromX MChromX - 1
     set karyotype karyotype - 1]
    [set MChromX MChromX + 1
     set karyotype karyotype + 1]
  set mmisegX 1 ]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if PP1 < 100 and PChrom1 > 0
   [
    ifelse PGL1 = 0
    [set PChrom1 PChrom1 - 1
     set karyotype karyotype - 1]
    [set PChrom1 PChrom1 + 1
     set karyotype karyotype + 1]
  set pmiseg1 1 ]

  if PP2 < 100 and PChrom2 > 0
   [
    ifelse PGL2 = 0
    [set PChrom2 PChrom2 -  1
     set karyotype karyotype -  1]
    [set PChrom2 PChrom2 +  1
     set karyotype karyotype + 1]
  set pmiseg2 1 ]

  if PP3 < 100 and PChrom3 > 0
   [
    ifelse PGL3 = 0
    [set PChrom3 PChrom3 -  1
     set karyotype karyotype -  1]
    [set PChrom3 PChrom3 +  1
     set karyotype karyotype + 1]
  set pmiseg3 1 ]

  if PP4 < 100 and PChrom4 > 0
   [
    ifelse PGL4 = 0
    [set PChrom4 PChrom4 -  1
     set karyotype karyotype -  1]
    [set PChrom4 PChrom4 +  1
     set karyotype karyotype + 1]
  set pmiseg4 1 ]

  if PP5 < 100 and PChrom5 > 0
   [
    ifelse PGL5 = 0
    [set PChrom5 PChrom5 -  1
     set karyotype karyotype -  1]
    [set PChrom5 PChrom5 +  1
     set karyotype karyotype + 1]
  set pmiseg5 1 ]

  if PP6 < 100 and PChrom6 > 0
   [
    ifelse PGL6 = 0
    [set PChrom6 PChrom6 -  1
     set karyotype karyotype -  1]
    [set PChrom6 PChrom6 +  1
     set karyotype karyotype + 1]
  set pmiseg6 1 ]

  if PP7 < 100 and PChrom7 > 0
   [
    ifelse PGL7 = 0
    [set PChrom7 PChrom7 -  1
     set karyotype karyotype -  1]
    [set PChrom7 PChrom7 +  1
     set karyotype karyotype + 1]
  set pmiseg7 1 ]

  if PP8 < 100 and PChrom8 > 0
   [
    ifelse PGL8 = 0
    [set PChrom8 PChrom8 -  1
     set karyotype karyotype -  1]
    [set PChrom8 PChrom8 +  1
     set karyotype karyotype + 1]
  set pmiseg8 1 ]

  if PP9 < 100 and PChrom9 > 0
   [
    ifelse PGL9 = 0
    [set PChrom9 PChrom9 -  1
     set karyotype karyotype -  1]
    [set PChrom9 PChrom9 +  1
     set karyotype karyotype + 1]
  set pmiseg9 1 ]

  if PP10 < 100 and PChrom10 > 0
   [
    ifelse PGL10 = 0
    [set PChrom10 PChrom10 -  1
     set karyotype karyotype -  1]
    [set PChrom10 PChrom10 +  1
     set karyotype karyotype + 1]
  set pmiseg10 1 ]

  if PP11 < 100 and PChrom11 > 0
   [
    ifelse PGL11 = 0
    [set PChrom11 PChrom11 -  1
     set karyotype karyotype -  1]
    [set PChrom11 PChrom11 +  1
     set karyotype karyotype + 1]
  set pmiseg11 1 ]

  if PP12 < 100 and PChrom12 > 0
   [
    ifelse PGL12 = 0
    [set PChrom12 PChrom12 -  1
     set karyotype karyotype -  1]
    [set PChrom12 PChrom12 +  1
     set karyotype karyotype + 1]
  set pmiseg12 1 ]

  if PP13 < 100 and PChrom13 > 0
   [
    ifelse PGL13 = 0
    [set PChrom13 PChrom13 -  1
     set karyotype karyotype -  1]
    [set PChrom13 PChrom13 +  1
     set karyotype karyotype + 1]
  set pmiseg13 1 ]

  if PP14 < 100 and PChrom14 > 0
   [
    ifelse PGL14 = 0
    [set PChrom14 PChrom14 -  1
     set karyotype karyotype -  1]
    [set PChrom14 PChrom14 +  1
     set karyotype karyotype + 1]
  set pmiseg14 1 ]

  if PP15 < 100 and PChrom15 > 0
   [
    ifelse PGL15 = 0
    [set PChrom15 PChrom15 -  1
     set karyotype karyotype -  1]
    [set PChrom15 PChrom15 +  1
     set karyotype karyotype + 1]
  set pmiseg15 1 ]

  if PP16 < 100 and PChrom16 > 0
   [
    ifelse PGL16 = 0
    [set PChrom16 PChrom16 -  1
     set karyotype karyotype -  1]
    [set PChrom16 PChrom16 +  1
     set karyotype karyotype + 1]
  set pmiseg16 1 ]

  if PP17 < 100 and PChrom17 > 0
   [
    ifelse PGL17 = 0
    [set PChrom17 PChrom17 -  1
     set karyotype karyotype -  1]
    [set PChrom17 PChrom17 +  1
     set karyotype karyotype + 1]
  set pmiseg17 1 ]


  if PP18 < 100 and PChrom18 > 0
   [
    ifelse PGL18 = 0
    [set PChrom18 PChrom18 -  1
     set karyotype karyotype -  1]
    [set PChrom18 PChrom18 +  1
     set karyotype karyotype + 1]
  set pmiseg18 1 ]

  if PP19 < 100 and PChrom19 > 0
   [
    ifelse PGL19 = 0
    [set PChrom19 PChrom19 -  1
     set karyotype karyotype -  1]
    [set PChrom19 PChrom19 +  1
     set karyotype karyotype + 1]
  set pmiseg19 1 ]

  if PP20 < 100 and PChrom20 > 0
   [
    ifelse PGL20 = 0
    [set PChrom20 PChrom20 -  1
     set karyotype karyotype -  1]
    [set PChrom20 PChrom20 +  1
     set karyotype karyotype + 1]
  set pmiseg20 1 ]

  if PP21 < 100 and PChrom21 > 0
   [
    ifelse PGL21 = 0
    [set PChrom21 PChrom21 -  1
     set karyotype karyotype -  1]
    [set PChrom21 PChrom21 +  1
     set karyotype karyotype + 1]
  set pmiseg21 1 ]

  if PP22 < 100 and PChrom22 > 0
   [
    ifelse PGL22 = 0
    [set PChrom22 PChrom22 -  1
     set karyotype karyotype -  1]
    [set PChrom22 PChrom22 +  1
     set karyotype karyotype + 1]
  set pmiseg22 1 ]

  if PPX < 100 and PChromX > 0
   [
    ifelse PGLX = 0
    [set PChromX PChromX -  1
     set karyotype karyotype -  1]
    [set PChromX PChromX +  1
     set karyotype karyotype + 1]
  set pmisegX 1 ]


    set Chrom1 (MChrom1 + PChrom1)
    set Chrom2 (MChrom2 + PChrom2)
    set Chrom3 (MChrom3 + PChrom3)
    set Chrom4 (MChrom4 + PChrom4)
    set Chrom5 (MChrom5 + PChrom5)
    set Chrom6 (MChrom6 + PChrom6)
    set Chrom7 (MChrom7 + PChrom7)
    set Chrom8 (MChrom8 + PChrom8)
    set Chrom9 (MChrom9 + PChrom9)
    set Chrom10 (MChrom10 + PChrom10)
    set Chrom11 (MChrom11 + PChrom11)
    set Chrom12 (MChrom12 + PChrom12)
    set Chrom13 (MChrom13 + PChrom13)
    set Chrom14 (MChrom14 + PChrom14)
    set Chrom15 (MChrom15 + PChrom15)
    set Chrom16 (MChrom16 + PChrom16)
    set Chrom17 (MChrom17 + PChrom17)
    set Chrom18 (MChrom18 + PChrom18)
    set Chrom19 (MChrom19 + PChrom19)
    set Chrom20 (MChrom20 + PChrom20)
    set Chrom21 (MChrom21 + PChrom21)
    set Chrom22 (MChrom22 + PChrom22)
    set ChromX (MChromX + PChromX)
end

to mitosis2
    hatch 1
   [
     if MP1 < 100 and mmiseg1 = 1
   [
    ifelse MGL1 = 0
    [set MChrom1 Mchrom1 + 2
     set karyotype karyotype + 2]
    [set MChrom1 MChrom1 - 2
     set karyotype karyotype - 2
  ]]

   if MP2 < 100 and mmiseg2 = 1
   [
    ifelse MGL2 = 0
    [set MChrom2 MChrom2 + 2
     set karyotype karyotype + 2]
    [set MChrom2 MChrom2 - 2
     set karyotype karyotype - 2
  ]]

   if MP3 < 100 and mmiseg3 = 1
   [
    ifelse MGL3 = 0
    [set MChrom3 MChrom3 + 2
     set karyotype karyotype + 2]
    [set MChrom3 MChrom3 - 2
     set karyotype karyotype - 2
  ]]

   if MP4 < 100 and mmiseg4 = 1
   [
    ifelse MGL4 = 0
    [set MChrom4 MChrom4 + 2
     set karyotype karyotype + 2]
    [set MChrom4 MChrom4 - 2
     set karyotype karyotype - 2
  ]]

   if MP5 < 100 and mmiseg5 = 1
   [
    ifelse MGL5 = 0
    [set MChrom5 MChrom5 + 2
     set karyotype karyotype + 2]
    [set MChrom5 MChrom5 - 2
     set karyotype karyotype - 2
  ]]

   if MP6 < 100 and mmiseg6 = 1
   [
    ifelse MGL6 = 0
    [set MChrom6 MChrom6 + 2
     set karyotype karyotype + 2]
    [set MChrom6 MChrom6 - 2
     set karyotype karyotype - 2
  ]]

   if MP7 < 100 and mmiseg7 = 1
   [
    ifelse MGL7 = 0
    [set MChrom7 MChrom7 + 2
     set karyotype karyotype + 2]
    [set MChrom7 MChrom7 - 2
     set karyotype karyotype - 2
  ]]

   if MP8 < 100 and mmiseg8 = 1
   [
    ifelse MGL8 = 0
    [set MChrom8 MChrom8 + 2
     set karyotype karyotype + 2]
    [set MChrom8 MChrom8 - 2
     set karyotype karyotype - 2
  ]]

   if MP9 < 100 and mmiseg9 = 1
   [
    ifelse MGL9 = 0
    [set MChrom9 MChrom9 + 2
     set karyotype karyotype + 2]
    [set MChrom9 MChrom9 - 2
     set karyotype karyotype - 2
  ]]

   if MP10 < 100 and mmiseg10 = 1
   [
    ifelse MGL10 = 0
    [set MChrom10 MChrom10 + 2
     set karyotype karyotype + 2]
    [set MChrom10 MChrom10 - 2
     set karyotype karyotype - 2
  ]]

   if MP11 < 100 and mmiseg11 = 1
   [
    ifelse MGL11 = 0
    [set MChrom11 MChrom11 + 2
     set karyotype karyotype + 2]
    [set MChrom11 MChrom11 - 2
     set karyotype karyotype - 2
  ]]

   if MP12 < 100 and mmiseg12 = 1
   [
    ifelse MGL12 = 0
    [set MChrom12 MChrom12 + 2
     set karyotype karyotype + 2]
    [set MChrom12 MChrom12 - 2
     set karyotype karyotype - 2
  ]]

   if MP13 < 100 and mmiseg13 = 1
   [
    ifelse MGL13 = 0
    [set MChrom13 MChrom13 + 2
     set karyotype karyotype + 2]
    [set MChrom13 MChrom13 - 2
     set karyotype karyotype - 2
  ]]

   if MP14 < 100 and mmiseg14 = 1
   [
    ifelse MGL14 = 0
    [set MChrom14 MChrom14 + 2
     set karyotype karyotype + 2]
    [set MChrom14 MChrom14 - 2
     set karyotype karyotype - 2
  ]]

   if MP15 < 100 and mmiseg15 = 1
   [
    ifelse MGL15 = 0
    [set MChrom15 MChrom15 + 2
     set karyotype karyotype + 2]
    [set MChrom15 MChrom15 - 2
     set karyotype karyotype - 2
  ]]

   if MP16 < 100 and mmiseg16 = 1
   [
    ifelse MGL16 = 0
    [set MChrom16 MChrom16 + 2
     set karyotype karyotype + 2]
    [set MChrom16 MChrom16 - 2
     set karyotype karyotype - 2
  ]]

   if MP17 < 100 and mmiseg17 = 1
   [
    ifelse MGL17 = 0
    [set MChrom17 MChrom17 + 2
     set karyotype karyotype + 2]
    [set MChrom17 MChrom17 - 2
     set karyotype karyotype - 2
  ]]

   if MP18 < 100 and mmiseg18 = 1
   [
    ifelse MGL18 = 0
    [set MChrom18 MChrom18 + 2
     set karyotype karyotype + 2]
    [set MChrom18 MChrom18 - 2
     set karyotype karyotype - 2
  ]]

   if MP19 < 100 and mmiseg19 = 1
   [
    ifelse MGL19 = 0
    [set MChrom19 MChrom19 + 2
     set karyotype karyotype + 2]
    [set MChrom19 MChrom19 - 2
     set karyotype karyotype - 2
  ]]

   if MP20 < 100 and mmiseg20 = 1
   [
    ifelse MGL20 = 0
    [set MChrom20 MChrom20 + 2
     set karyotype karyotype + 2]
    [set MChrom20 MChrom20 - 2
     set karyotype karyotype - 2
  ]]

   if MP21 < 100 and mmiseg21 = 1
   [
    ifelse MGL21 = 0
    [set MChrom21 MChrom21 + 2
     set karyotype karyotype + 2]
    [set MChrom21 MChrom21 - 2
     set karyotype karyotype - 2
  ]]

   if MP22 < 100 and mmiseg22 = 1
   [
    ifelse MGL22 = 0
    [set MChrom22 MChrom22 + 2
     set karyotype karyotype + 2]
    [set MChrom22 MChrom22 - 2
     set karyotype karyotype - 2
  ]]

   if MPX < 100 and mmisegX = 1
   [
    ifelse MGLX = 0
    [set MChromX MChromX + 2
     set karyotype karyotype + 2]
    [set MChromX MChromX - 2
     set karyotype karyotype - 2
  ]]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if PP1 < 100 and pmiseg1 = 1
   [
    ifelse PGL1 = 0
    [set PChrom1 PChrom1 + 2
     set karyotype karyotype + 2]
    [set PChrom1 PChrom1 - 2
     set karyotype karyotype - 2
  ]]

  if PP2 < 100 and pmiseg2 = 1
   [
    ifelse PGL2 = 0
    [set PChrom2 PChrom2 + 2
     set karyotype karyotype + 2]
    [set PChrom2 PChrom2 - 2
     set karyotype karyotype - 2
  ]]

  if PP3 < 100 and pmiseg3 = 1
   [
    ifelse PGL3 = 0
    [set PChrom3 PChrom3 + 2
     set karyotype karyotype + 2]
    [set PChrom3 PChrom3 - 2
     set karyotype karyotype - 2
  ]]

  if PP4 < 100 and pmiseg4 = 1
   [
    ifelse PGL4 = 0
    [set PChrom4 PChrom4 + 2
     set karyotype karyotype + 2]
    [set PChrom4 PChrom4 - 2
     set karyotype karyotype - 2
  ]]

  if PP5 < 100 and pmiseg5 = 1
   [
    ifelse PGL5 = 0
    [set PChrom5 PChrom5 + 2
     set karyotype karyotype + 2]
    [set PChrom5 PChrom5 - 2
     set karyotype karyotype - 2
  ]]

  if PP6 < 100 and pmiseg6 = 1
   [
    ifelse PGL6 = 0
    [set PChrom6 PChrom6 + 2
     set karyotype karyotype + 2]
    [set PChrom6 PChrom6 - 2
     set karyotype karyotype - 2
  ]]

  if PP7 < 100 and pmiseg7 = 1
   [
    ifelse PGL7 = 0
    [set PChrom7 PChrom7 + 2
     set karyotype karyotype + 2]
    [set PChrom7 PChrom7 - 2
     set karyotype karyotype - 2
  ]]

  if PP8 < 100 and pmiseg8 = 1
   [
    ifelse PGL8 = 0
    [set PChrom8 PChrom8 + 2
     set karyotype karyotype + 2]
    [set PChrom8 PChrom8 - 2
     set karyotype karyotype - 2
  ]]

  if PP9 < 100 and pmiseg9 = 1
   [
    ifelse PGL9 = 0
    [set PChrom9 PChrom9 + 2
     set karyotype karyotype + 2]
    [set PChrom9 PChrom9 - 2
     set karyotype karyotype - 2
  ]]

  if PP10 < 100 and pmiseg10 = 1
   [
    ifelse PGL10 = 0
    [set PChrom10 PChrom10 + 2
     set karyotype karyotype + 2]
    [set PChrom10 PChrom10 - 2
     set karyotype karyotype - 2
  ]]

  if PP11 < 100 and pmiseg11 = 1
   [
    ifelse PGL11 = 0
    [set PChrom11 PChrom11 + 2
     set karyotype karyotype + 2]
    [set PChrom11 PChrom11 - 2
     set karyotype karyotype - 2
  ]]

  if PP12 < 100 and pmiseg12 = 1
   [
    ifelse PGL12 = 0
    [set PChrom12 PChrom12 + 2
     set karyotype karyotype + 2]
    [set PChrom12 PChrom12 - 2
     set karyotype karyotype - 2
  ]]

  if PP13 < 100 and pmiseg13 = 1
   [
    ifelse PGL13 = 0
    [set PChrom13 PChrom13 + 2
     set karyotype karyotype + 2]
    [set PChrom13 PChrom13 - 2
     set karyotype karyotype - 2
  ]]

  if PP14 < 100 and pmiseg14 = 1
   [
    ifelse PGL14 = 0
    [set PChrom14 PChrom14 + 2
     set karyotype karyotype + 2]
    [set PChrom14 PChrom14 - 2
     set karyotype karyotype - 2
  ]]

  if PP15 < 100 and pmiseg15 = 1
   [
    ifelse PGL15 = 0
    [set PChrom15 PChrom15 + 2
     set karyotype karyotype + 2]
    [set PChrom15 PChrom15 - 2
     set karyotype karyotype - 2
  ]]

  if PP16 < 100 and pmiseg16 = 1
   [
    ifelse PGL16 = 0
    [set PChrom16 PChrom16 + 2
     set karyotype karyotype + 2]
    [set PChrom16 PChrom16 - 2
     set karyotype karyotype - 2
  ]]

  if PP17 < 100 and pmiseg17 = 1
   [
    ifelse PGL17 = 0
    [set PChrom17 PChrom17 + 2
     set karyotype karyotype + 2]
    [set PChrom17 PChrom17 - 2
     set karyotype karyotype - 2
  ]]


  if PP18 < 100 and pmiseg18 = 1
   [
    ifelse PGL18 = 0
    [set PChrom18 PChrom18 + 2
     set karyotype karyotype + 2]
    [set PChrom18 PChrom18 - 2
     set karyotype karyotype - 2
  ]]

  if PP19 < 100 and pmiseg19 = 1
   [
    ifelse PGL19 = 0
    [set PChrom19 PChrom19 + 2
     set karyotype karyotype + 2]
    [set PChrom19 PChrom19 - 2
     set karyotype karyotype - 2
  ]]

  if PP20 < 100 and pmiseg20 = 1
   [
    ifelse PGL20 = 0
    [set PChrom20 PChrom20 + 2
     set karyotype karyotype + 2]
    [set PChrom20 PChrom20 - 2
     set karyotype karyotype - 2
  ]]

  if PP21 < 100 and pmiseg21 = 1
   [
    ifelse PGL21 = 0
    [set PChrom21 PChrom21 + 2
     set karyotype karyotype + 2]
    [set PChrom21 PChrom21 - 2
     set karyotype karyotype - 2
  ]]

  if PP22 < 100 and pmiseg22 = 1
   [
    ifelse PGL22 = 0
    [set PChrom22 PChrom22 + 2
     set karyotype karyotype + 2]
    [set PChrom22 PChrom22 - 2
     set karyotype karyotype - 2
  ]]

  if PPX < 100 and pmisegX = 1
   [
    ifelse PGLX = 0
    [set PChromX PChromX + 2
     set karyotype karyotype + 2]
    [set PChromX PChromX - 2
     set karyotype karyotype - 2
  ]]


    set Chrom1 (MChrom1 + PChrom1)
    set Chrom2 (MChrom2 + PChrom2)
    set Chrom3 (MChrom3 + PChrom3)
    set Chrom4 (MChrom4 + PChrom4)
    set Chrom5 (MChrom5 + PChrom5)
    set Chrom6 (MChrom6 + PChrom6)
    set Chrom7 (MChrom7 + PChrom7)
    set Chrom8 (MChrom8 + PChrom8)
    set Chrom9 (MChrom9 + PChrom9)
    set Chrom10 (MChrom10 + PChrom10)
    set Chrom11 (MChrom11 + PChrom11)
    set Chrom12 (MChrom12 + PChrom12)
    set Chrom13 (MChrom13 + PChrom13)
    set Chrom14 (MChrom14 + PChrom14)
    set Chrom15 (MChrom15 + PChrom15)
    set Chrom16 (MChrom16 + PChrom16)
    set Chrom17 (MChrom17 + PChrom17)
    set Chrom18 (MChrom18 + PChrom18)
    set Chrom19 (MChrom19 + PChrom19)
    set Chrom20 (MChrom20 + PChrom20)
    set Chrom21 (MChrom21 + PChrom21)
    set Chrom22 (MChrom22 + PChrom22)
    set ChromX (MChromX + PChromX)

    set-color
    set-position
    ]

end

to set-position
      rt random 90
      lt random 90
      fd 1.00
end

to calculate-chromosome-scores
  ask tumorcells [
  set karyotype sum(list
      Chrom1
      Chrom2
      Chrom3
      Chrom4
      Chrom5
      Chrom6
      Chrom7
      Chrom8
      Chrom9
      Chrom10
      Chrom11
      Chrom12
      Chrom13
      Chrom14
      Chrom15
      Chrom16
      Chrom17
      Chrom18
      Chrom19
      Chrom20
      Chrom21
      Chrom22
      ChromX
      )

  set chromavg karyotype / 23

  set fitd1 (Fc1 - (Fc1 * abs(Chrom1 - chromavg)) / chromavg)
  set fitd2 (Fc2 - (Fc2 * abs(Chrom2 - chromavg)) / chromavg)
  set fitd3 (Fc3 - (Fc3 * abs(Chrom3 - chromavg)) / chromavg)
  set fitd4 (Fc4 - (Fc4 * abs(Chrom4 - chromavg)) / chromavg)
  set fitd5 (Fc5 - (Fc5 * abs(Chrom5 - chromavg)) / chromavg)
  set fitd6 (Fc6 - (Fc6 * abs(Chrom6 - chromavg)) / chromavg)
  set fitd7 (Fc7 - (Fc7 * abs(Chrom7 - chromavg)) / chromavg)
  set fitd8 (Fc8 - (Fc8 * abs(Chrom8 - chromavg)) / chromavg)
  set fitd9 (Fc9 - (Fc9 * abs(Chrom9 - chromavg)) / chromavg)
  set fitd10 (Fc10 - (Fc10 * abs(Chrom10 - chromavg)) / chromavg)
  set fitd11 (Fc11 - (Fc11 * abs(Chrom11 - chromavg)) / chromavg)
  set fitd12 (Fc12 - (Fc12 * abs(Chrom12 - chromavg)) / chromavg)
  set fitd13 (Fc13 - (Fc13 * abs(Chrom13 - chromavg)) / chromavg)
  set fitd14 (Fc14 - (Fc14 * abs(Chrom14 - chromavg)) / chromavg)
  set fitd15 (Fc15 - (Fc15 * abs(Chrom15 - chromavg)) / chromavg)
  set fitd16 (Fc16 - (Fc16 * abs(Chrom16 - chromavg)) / chromavg)
  set fitd17 (Fc17 - (Fc17 * abs(Chrom17 - chromavg)) / chromavg)
  set fitd18 (Fc18 - (Fc18 * abs(Chrom18 - chromavg)) / chromavg)
  set fitd19 (Fc19 - (Fc19 * abs(Chrom19 - chromavg)) / chromavg)
  set fitd20 (Fc20 - (Fc20 * abs(Chrom20 - chromavg)) / chromavg)
  set fitd21 (Fc21 - (Fc21 * abs(Chrom21 - chromavg)) / chromavg)
  set fitd22 (Fc22 - (Fc22 * abs(Chrom22 - chromavg)) / chromavg)
  set fitdx (Fcx - (Fcx * abs(ChromX - chromavg)) / chromavg)
  ]

end

to calculate-fitness
if model = "Hybrid" [
    calculate-fitness-hybrid
  ]
if model = "Abundance" [
    calculate-fitness-abundance
  ]
if model = "Driver" [
    calculate-fitness-driver
  ]
end

to calculate-fitness-hybrid ;----------------------------------------fitness values
  ask tumorcells [
    set fitness (precision (((
      fitd1 +
      fitd2 +
      fitd3 +
      fitd4 +
      fitd5 +
      fitd6 +
      fitd7 +
      fitd8 +
      fitd9 +
      fitd10 +
      fitd11 +
      fitd12 +
      fitd13 +
      fitd14 +
      fitd15 +
      fitd16 +
      fitd17 +
      fitd18 +
      fitd19 +
      fitd20 +
      fitd21 +
      fitd22 +
      fitdx) + (((Chrom1 * CS1) / chromavg) +
      ((Chrom2 * CS2) / chromavg) +
      ((Chrom3 * CS3) / chromavg) +
      ((Chrom4 * CS4) / chromavg) +
      ((Chrom5 * CS5) / chromavg) +
      ((Chrom6 * CS6) / chromavg) +
      ((Chrom7 * CS7) / chromavg) +
      ((Chrom8 * CS8) / chromavg) +
      ((Chrom9 * CS9) / chromavg) +
      ((Chrom10 * CS10) / chromavg) +
      ((Chrom11 * CS11) / chromavg) +
      ((Chrom12 * CS12) / chromavg) +
      ((Chrom13 * CS13) / chromavg) +
      ((Chrom14 * CS14) / chromavg) +
      ((Chrom15 * CS15) / chromavg) +
      ((Chrom16 * CS16) / chromavg) +
      ((Chrom17 * CS17) / chromavg) +
      ((Chrom18 * CS18) / chromavg) +
      ((Chrom19 * CS19) / chromavg) +
      ((Chrom20 * CS20) / chromavg) +
      ((Chrom21 * CS21) / chromavg) +
      ((Chrom22 * CS22) / chromavg))) / 2) 3) ^ selective-pressure
  ]
end

to calculate-fitness-driver
  ask tumorcells [
    set fitness (precision ((
      ((Chrom1 * CS1) / chromavg) +
      ((Chrom2 * CS2) / chromavg) +
      ((Chrom3 * CS3) / chromavg) +
      ((Chrom4 * CS4) / chromavg) +
      ((Chrom5 * CS5) / chromavg) +
      ((Chrom6 * CS6) / chromavg) +
      ((Chrom7 * CS7) / chromavg) +
      ((Chrom8 * CS8) / chromavg) +
      ((Chrom9 * CS9) / chromavg) +
      ((Chrom10 * CS10) / chromavg) +
      ((Chrom11 * CS11) / chromavg) +
      ((Chrom12 * CS12) / chromavg) +
      ((Chrom13 * CS13) / chromavg) +
      ((Chrom14 * CS14) / chromavg) +
      ((Chrom15 * CS15) / chromavg) +
      ((Chrom16 * CS16) / chromavg) +
      ((Chrom17 * CS17) / chromavg) +
      ((Chrom18 * CS18) / chromavg) +
      ((Chrom19 * CS19) / chromavg) +
      ((Chrom20 * CS20) / chromavg) +
      ((Chrom21 * CS21) / chromavg) +
      ((Chrom22 * CS22) / chromavg))) 3) ^ selective-pressure
  ]
end

to calculate-fitness-abundance
ask tumorcells [
   set fitness (precision(fitd1 +
    fitd2 +
    fitd3 +
    fitd4 +
    fitd5 +
    fitd6 +
    fitd7 +
    fitd8 +
    fitd9 +
    fitd10 +
    fitd11 +
    fitd12 +
    fitd13 +
    fitd14 +
    fitd15 +
    fitd16 +
    fitd17 +
    fitd18 +
    fitd19 +
    fitd20 +
    fitd21 +
    fitd22 +
    fitdx) 3) ^ selective-pressure
  ]
end

to check-death

  ask tumorcells [
    if fitness <= 0 [die]

    if Chrom1 <= 0  [die]
    if Chrom2 <= 0  [die]
    if Chrom3 <= 0  [die]
    if Chrom4 <= 0  [die]
    if Chrom5 <= 0  [die]
    if Chrom6 <= 0  [die]
    if Chrom7 <= 0  [die]
    if Chrom8 <= 0  [die]
    if Chrom9 <= 0  [die]
    if Chrom10 <= 0 [die]
    if Chrom11 <= 0 [die]
    if Chrom12 <= 0 [die]
    if Chrom13 <= 0 [die]
    if Chrom14 <= 0 [die]
    if Chrom15 <= 0 [die]
    if Chrom16 <= 0 [die]
    if Chrom17 <= 0 [die]
    if Chrom18 <= 0 [die]
    if Chrom19 <= 0 [die]
    if Chrom20 <= 0 [die]
    if Chrom21 <= 0 [die]
    if Chrom22 <= 0 [die]
    if ChromX <= 0  [die]

    if Chrom1 > 6  [die]
    if Chrom2 > 6  [die]
    if Chrom3 > 6  [die]
    if Chrom4 > 6  [die]
    if Chrom5 > 6  [die]
    if Chrom6 > 6  [die]
    if Chrom7 > 6  [die]
    if Chrom8 > 6  [die]
    if Chrom9 > 6  [die]
    if Chrom10 > 6 [die]
    if Chrom11 > 6 [die]
    if Chrom12 > 6 [die]
    if Chrom13 > 6 [die]
    if Chrom14 > 6 [die]
    if Chrom15 > 6 [die]
    if Chrom16 > 6 [die]
    if Chrom17 > 6 [die]
    if Chrom18 > 6 [die]
    if Chrom19 > 6 [die]
    if Chrom20 > 6 [die]
    if Chrom21 > 6 [die]
    if Chrom22 > 6 [die]
    if ChromX > 6  [die]
  ]
  ask n-of (count tumorcells / 2) tumorcells [
    if count tumorcells >= 3000 [die]
  ]
end

to output-chromosome-state
 calculate-chromosome-scores
 calculate-fitness
 set x count tumorcells

  let selected-tumorcells ifelse-value (x >= outputcount)
  [n-of outputcount tumorcells]
  [n-of count tumorcells tumorcells]
  ask selected-tumorcells [
    set color 95
    set color [6 175 201 255]
  ]

  set Chrom1-final [Chrom1] of selected-tumorcells
  set Chrom2-final [Chrom2] of selected-tumorcells
  set Chrom3-final [Chrom3] of selected-tumorcells
  set Chrom4-final [Chrom4] of selected-tumorcells
  set Chrom5-final [Chrom5] of selected-tumorcells
  set Chrom6-final [Chrom6] of selected-tumorcells
  set Chrom7-final [Chrom7] of selected-tumorcells
  set Chrom8-final [Chrom8] of selected-tumorcells
  set Chrom9-final [Chrom9] of selected-tumorcells
  set Chrom10-final [Chrom10] of selected-tumorcells
  set Chrom11-final [Chrom11] of selected-tumorcells
  set Chrom12-final [Chrom12] of selected-tumorcells
  set Chrom13-final [Chrom13] of selected-tumorcells
  set Chrom14-final [Chrom14] of selected-tumorcells
  set Chrom15-final [Chrom15] of selected-tumorcells
  set Chrom16-final [Chrom16] of selected-tumorcells
  set Chrom17-final [Chrom17] of selected-tumorcells
  set Chrom18-final [Chrom18] of selected-tumorcells
  set Chrom19-final [Chrom19] of selected-tumorcells
  set Chrom20-final [Chrom20] of selected-tumorcells
  set Chrom21-final [Chrom21] of selected-tumorcells
  set Chrom22-final [Chrom22] of selected-tumorcells
  set ChromX-final [ChromX] of selected-tumorcells
  set fitness-final [fitness] of selected-tumorcells
  if count tumorcells > 1
  [set MKV mean(list
    variance [Chrom1] of selected-tumorcells
    variance [Chrom2] of selected-tumorcells
    variance [Chrom3] of selected-tumorcells
    variance [Chrom4] of selected-tumorcells
    variance [Chrom5] of selected-tumorcells
    variance [Chrom6] of selected-tumorcells
    variance [Chrom7] of selected-tumorcells
    variance [Chrom8] of selected-tumorcells
    variance [Chrom9] of selected-tumorcells
    variance [Chrom10] of selected-tumorcells
    variance [Chrom11] of selected-tumorcells
    variance [Chrom12] of selected-tumorcells
    variance [Chrom13] of selected-tumorcells
    variance [Chrom14] of selected-tumorcells
    variance [Chrom15] of selected-tumorcells
    variance [Chrom16] of selected-tumorcells
    variance [Chrom17] of selected-tumorcells
    variance [Chrom18] of selected-tumorcells
    variance [Chrom19] of selected-tumorcells
    variance [Chrom20] of selected-tumorcells
    variance [Chrom21] of selected-tumorcells
    variance [Chrom22] of selected-tumorcells
    variance [ChromX] of selected-tumorcells
    )
  ]
end
@#$#@#$#@
GRAPHICS-WINDOW
335
13
663
342
-1
-1
7.805
1
10
1
1
1
0
1
1
1
-20
20
-20
20
1
1
1
ticks
30.0

BUTTON
29
30
95
63
Setup
setup\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
96
30
159
63
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
161
30
224
63
Step
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
226
19
304
64
Cell Count
count turtles
0
1
11

INPUTBOX
30
65
179
125
initial-number
100.0
1
0
Number

INPUTBOX
30
365
179
425
initial-CIN
511.0
1
0
Number

SLIDER
178
306
328
339
initial-ploidy
initial-ploidy
2
6
4.0
2
1
NIL
HORIZONTAL

INPUTBOX
30
125
179
185
end-ticks
80.0
1
0
Number

INPUTBOX
30
245
179
305
outputcount
300.0
1
0
Number

INPUTBOX
179
170
328
230
cin-period-end
80.0
1
0
Number

INPUTBOX
30
185
179
245
end-population
100000.0
1
0
Number

INPUTBOX
179
110
328
170
cin-period-start
0.0
1
0
Number

CHOOSER
179
65
317
110
periodicity
periodicity
"constant" "early" "middle" "late"
0

INPUTBOX
30
305
179
365
mitosis-probability
200.0
1
0
Number

CHOOSER
179
230
328
275
model
model
"Hybrid" "Driver" "Abundance"
0

SLIDER
179
274
328
307
selective-pressure
selective-pressure
0
20
0.0
1
1
NIL
HORIZONTAL

PLOT
230
425
430
575
Average Fitness
NIL
NIL
0.0
10.0
0.0
2.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [ fitness ] of tumorcells"

PLOT
30
425
230
575
Population
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"All cells" 1.0 0 -16777216 true "" "plot count tumorcells"

PLOT
430
425
630
575
Average Ploidy
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot mean [ karyotype ] of tumorcells / 23"

PLOT
701
69
1055
271
Copy Numbers
NIL
NIL
0.0
10.0
0.0
6.0
true
true
"" ""
PENS
"1" 1.0 0 -16777216 true "" "plot mean [ Chrom1 ] of tumorcells "
"X" 1.0 0 -7500403 true "" "plot mean [ ChromX ] of tumorcells"
"2" 1.0 0 -2674135 true "" "plot mean [ Chrom2 ] of tumorcells"
"3" 1.0 0 -955883 true "" "plot mean [ Chrom3 ] of tumorcells"
"4" 1.0 0 -6459832 true "" "plot mean [ Chrom4 ] of tumorcells"
"5" 1.0 0 -1184463 true "" "plot mean [ Chrom5 ] of tumorcells"
"6" 1.0 0 -10899396 true "" "plot mean [ Chrom6 ] of tumorcells"
"7" 1.0 0 -13840069 true "" "plot mean [ Chrom7 ] of tumorcells"
"8" 1.0 0 -14835848 true "" "plot mean [ Chrom8 ] of tumorcells"
"9" 1.0 0 -11221820 true "" "plot mean [ Chrom9 ] of tumorcells"
"10" 1.0 0 -13791810 true "" "plot mean [ Chrom10 ] of tumorcells"
"11" 1.0 0 -13345367 true "" "plot mean [ Chrom11 ] of tumorcells"
"12" 1.0 0 -8630108 true "" "plot mean [ Chrom12 ] of tumorcells"
"13" 1.0 0 -5825686 true "" "plot mean [ Chrom13 ] of tumorcells"

PLOT
383
575
630
777
Tandem Losses
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot tandem-losses"

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

circle-outline
false
0
Circle -1 true false 0 0 300
Circle -16777216 false false -2 -2 302

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="fullexperiment" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>CIN-rate</metric>
    <metric>Chrom1-final</metric>
    <metric>Chrom2-final</metric>
    <metric>Chrom3-final</metric>
    <metric>Chrom4-final</metric>
    <metric>Chrom5-final</metric>
    <metric>Chrom6-final</metric>
    <metric>Chrom7-final</metric>
    <metric>Chrom8-final</metric>
    <metric>Chrom9-final</metric>
    <metric>Chrom10-final</metric>
    <metric>Chrom11-final</metric>
    <metric>Chrom12-final</metric>
    <metric>Chrom13-final</metric>
    <metric>Chrom14-final</metric>
    <metric>Chrom15-final</metric>
    <metric>Chrom16-final</metric>
    <metric>Chrom17-final</metric>
    <metric>Chrom18-final</metric>
    <metric>Chrom19-final</metric>
    <metric>Chrom20-final</metric>
    <metric>Chrom21-final</metric>
    <metric>Chrom22-final</metric>
    <metric>ChromX-final</metric>
    <metric>fitness-final</metric>
    <enumeratedValueSet variable="selective-pressure">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
      <value value="8"/>
      <value value="9"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosis-probability">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-CIN">
      <value value="200"/>
      <value value="209"/>
      <value value="219"/>
      <value value="230"/>
      <value value="242"/>
      <value value="256"/>
      <value value="271"/>
      <value value="288"/>
      <value value="307"/>
      <value value="329"/>
      <value value="354"/>
      <value value="383"/>
      <value value="418"/>
      <value value="460"/>
      <value value="511"/>
      <value value="575"/>
      <value value="657"/>
      <value value="767"/>
      <value value="920"/>
      <value value="1150"/>
      <value value="1533"/>
      <value value="2300"/>
      <value value="4600"/>
      <value value="9200"/>
      <value value="13939"/>
      <value value="18400"/>
      <value value="23000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <steppedValueSet variable="end-ticks" first="1" step="2" last="80"/>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;early&quot;"/>
      <value value="&quot;middle&quot;"/>
      <value value="&quot;late&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="PopulationGrowthCINSelection" repetitions="8" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <enumeratedValueSet variable="selective-pressure">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cin-period-end">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosis-probability">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-CIN">
      <value value="0"/>
      <value value="230"/>
      <value value="2300"/>
      <value value="4600"/>
      <value value="9200"/>
      <value value="37000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cin-period-start">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-population">
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="SamplingExperiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>mean [ karyotype ] of tumorcells</metric>
    <metric>mean [ fitness ] of tumorcells</metric>
    <metric>CIN-rate</metric>
    <metric>Chrom1-final</metric>
    <metric>Chrom2-final</metric>
    <metric>Chrom3-final</metric>
    <metric>Chrom4-final</metric>
    <metric>Chrom5-final</metric>
    <metric>Chrom6-final</metric>
    <metric>Chrom7-final</metric>
    <metric>Chrom8-final</metric>
    <metric>Chrom9-final</metric>
    <metric>Chrom10-final</metric>
    <metric>Chrom11-final</metric>
    <metric>Chrom12-final</metric>
    <metric>Chrom13-final</metric>
    <metric>Chrom14-final</metric>
    <metric>Chrom15-final</metric>
    <metric>Chrom16-final</metric>
    <metric>Chrom17-final</metric>
    <metric>Chrom18-final</metric>
    <metric>Chrom19-final</metric>
    <metric>Chrom20-final</metric>
    <metric>Chrom21-final</metric>
    <metric>Chrom22-final</metric>
    <metric>ChromX-final</metric>
    <enumeratedValueSet variable="selective-pressure">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cin-period-end">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosis-probability">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-CIN">
      <value value="920"/>
      <value value="2300"/>
      <value value="4600"/>
      <value value="9200"/>
      <value value="37000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cin-period-start">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;early&quot;"/>
      <value value="&quot;late&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-population">
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="800"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="PeriodicityDiversity" repetitions="4" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>mean [ karyotype ] of tumorcells</metric>
    <metric>mean [ fitness ] of tumorcells</metric>
    <metric>CIN-rate</metric>
    <metric>Chrom1-final</metric>
    <metric>Chrom2-final</metric>
    <metric>Chrom3-final</metric>
    <metric>Chrom4-final</metric>
    <metric>Chrom5-final</metric>
    <metric>Chrom6-final</metric>
    <metric>Chrom7-final</metric>
    <metric>Chrom8-final</metric>
    <metric>Chrom9-final</metric>
    <metric>Chrom10-final</metric>
    <metric>Chrom11-final</metric>
    <metric>Chrom12-final</metric>
    <metric>Chrom13-final</metric>
    <metric>Chrom14-final</metric>
    <metric>Chrom15-final</metric>
    <metric>Chrom16-final</metric>
    <metric>Chrom17-final</metric>
    <metric>Chrom18-final</metric>
    <metric>Chrom19-final</metric>
    <metric>Chrom20-final</metric>
    <metric>Chrom21-final</metric>
    <metric>Chrom22-final</metric>
    <metric>ChromX-final</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-CIN">
      <value value="2300"/>
      <value value="4600"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="selective-pressure">
      <value value="0"/>
      <value value="2"/>
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosis-probability">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;early&quot;"/>
      <value value="&quot;middle&quot;"/>
      <value value="&quot;late&quot;"/>
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Abundance&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Hybrid&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="exampleData2" repetitions="8" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>MKV</metric>
    <metric>CIN-rate</metric>
    <metric>Chrom1-final</metric>
    <metric>Chrom2-final</metric>
    <metric>Chrom3-final</metric>
    <metric>Chrom4-final</metric>
    <metric>Chrom5-final</metric>
    <metric>Chrom6-final</metric>
    <metric>Chrom7-final</metric>
    <metric>Chrom8-final</metric>
    <metric>Chrom9-final</metric>
    <metric>Chrom10-final</metric>
    <metric>Chrom11-final</metric>
    <metric>Chrom12-final</metric>
    <metric>Chrom13-final</metric>
    <metric>Chrom14-final</metric>
    <metric>Chrom15-final</metric>
    <metric>Chrom16-final</metric>
    <metric>Chrom17-final</metric>
    <metric>Chrom18-final</metric>
    <metric>Chrom19-final</metric>
    <metric>Chrom20-final</metric>
    <metric>Chrom21-final</metric>
    <metric>Chrom22-final</metric>
    <metric>ChromX-final</metric>
    <metric>fitness-final</metric>
    <enumeratedValueSet variable="selective-pressure">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosis-probability">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-CIN">
      <value value="200"/>
      <value value="209"/>
      <value value="219"/>
      <value value="230"/>
      <value value="242"/>
      <value value="256"/>
      <value value="271"/>
      <value value="288"/>
      <value value="307"/>
      <value value="329"/>
      <value value="354"/>
      <value value="383"/>
      <value value="418"/>
      <value value="460"/>
      <value value="511"/>
      <value value="575"/>
      <value value="657"/>
      <value value="767"/>
      <value value="920"/>
      <value value="1150"/>
      <value value="1533"/>
      <value value="2300"/>
      <value value="4600"/>
      <value value="9200"/>
      <value value="13939"/>
      <value value="18400"/>
      <value value="23000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-population">
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="exampleData3" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>MKV</metric>
    <metric>CIN-rate</metric>
    <metric>Chrom1-final</metric>
    <metric>Chrom2-final</metric>
    <metric>Chrom3-final</metric>
    <metric>Chrom4-final</metric>
    <metric>Chrom5-final</metric>
    <metric>Chrom6-final</metric>
    <metric>Chrom7-final</metric>
    <metric>Chrom8-final</metric>
    <metric>Chrom9-final</metric>
    <metric>Chrom10-final</metric>
    <metric>Chrom11-final</metric>
    <metric>Chrom12-final</metric>
    <metric>Chrom13-final</metric>
    <metric>Chrom14-final</metric>
    <metric>Chrom15-final</metric>
    <metric>Chrom16-final</metric>
    <metric>Chrom17-final</metric>
    <metric>Chrom18-final</metric>
    <metric>Chrom19-final</metric>
    <metric>Chrom20-final</metric>
    <metric>Chrom21-final</metric>
    <metric>Chrom22-final</metric>
    <metric>ChromX-final</metric>
    <metric>fitness-final</metric>
    <enumeratedValueSet variable="selective-pressure">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="5"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-number">
      <value value="400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosis-probability">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-CIN">
      <value value="200"/>
      <value value="209"/>
      <value value="219"/>
      <value value="230"/>
      <value value="242"/>
      <value value="256"/>
      <value value="271"/>
      <value value="288"/>
      <value value="307"/>
      <value value="329"/>
      <value value="354"/>
      <value value="383"/>
      <value value="418"/>
      <value value="460"/>
      <value value="511"/>
      <value value="575"/>
      <value value="657"/>
      <value value="767"/>
      <value value="920"/>
      <value value="1150"/>
      <value value="1533"/>
      <value value="2300"/>
      <value value="4600"/>
      <value value="9200"/>
      <value value="13939"/>
      <value value="18400"/>
      <value value="23000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-population">
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-01-21_FullStepwiseExperiment" repetitions="4" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>CIN-rate</metric>
    <metric>Chrom1-final</metric>
    <metric>Chrom2-final</metric>
    <metric>Chrom3-final</metric>
    <metric>Chrom4-final</metric>
    <metric>Chrom5-final</metric>
    <metric>Chrom6-final</metric>
    <metric>Chrom7-final</metric>
    <metric>Chrom8-final</metric>
    <metric>Chrom9-final</metric>
    <metric>Chrom10-final</metric>
    <metric>Chrom11-final</metric>
    <metric>Chrom12-final</metric>
    <metric>Chrom13-final</metric>
    <metric>Chrom14-final</metric>
    <metric>Chrom15-final</metric>
    <metric>Chrom16-final</metric>
    <metric>Chrom17-final</metric>
    <metric>Chrom18-final</metric>
    <metric>Chrom19-final</metric>
    <metric>Chrom20-final</metric>
    <metric>Chrom21-final</metric>
    <metric>Chrom22-final</metric>
    <metric>ChromX-final</metric>
    <metric>fitness-final</metric>
    <metric>tandem-losses</metric>
    <steppedValueSet variable="selective-pressure" first="0" step="5" last="100"/>
    <enumeratedValueSet variable="initial-number">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosis-probability">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-CIN">
      <value value="200"/>
      <value value="209"/>
      <value value="219"/>
      <value value="230"/>
      <value value="242"/>
      <value value="256"/>
      <value value="271"/>
      <value value="288"/>
      <value value="307"/>
      <value value="329"/>
      <value value="354"/>
      <value value="383"/>
      <value value="418"/>
      <value value="460"/>
      <value value="511"/>
      <value value="575"/>
      <value value="657"/>
      <value value="767"/>
      <value value="920"/>
      <value value="1150"/>
      <value value="1533"/>
      <value value="2300"/>
      <value value="4600"/>
      <value value="9200"/>
      <value value="13939"/>
      <value value="18400"/>
      <value value="23000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;early&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-01-29_FullStepwiseExperiment" repetitions="4" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>CIN-rate</metric>
    <metric>Chrom1-final</metric>
    <metric>Chrom2-final</metric>
    <metric>Chrom3-final</metric>
    <metric>Chrom4-final</metric>
    <metric>Chrom5-final</metric>
    <metric>Chrom6-final</metric>
    <metric>Chrom7-final</metric>
    <metric>Chrom8-final</metric>
    <metric>Chrom9-final</metric>
    <metric>Chrom10-final</metric>
    <metric>Chrom11-final</metric>
    <metric>Chrom12-final</metric>
    <metric>Chrom13-final</metric>
    <metric>Chrom14-final</metric>
    <metric>Chrom15-final</metric>
    <metric>Chrom16-final</metric>
    <metric>Chrom17-final</metric>
    <metric>Chrom18-final</metric>
    <metric>Chrom19-final</metric>
    <metric>Chrom20-final</metric>
    <metric>Chrom21-final</metric>
    <metric>Chrom22-final</metric>
    <metric>ChromX-final</metric>
    <metric>fitness-final</metric>
    <metric>tandem-losses</metric>
    <steppedValueSet variable="selective-pressure" first="0" step="5" last="50"/>
    <enumeratedValueSet variable="initial-number">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosis-probability">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-CIN">
      <value value="200"/>
      <value value="219"/>
      <value value="242"/>
      <value value="271"/>
      <value value="307"/>
      <value value="354"/>
      <value value="418"/>
      <value value="511"/>
      <value value="657"/>
      <value value="920"/>
      <value value="1150"/>
      <value value="1533"/>
      <value value="2300"/>
      <value value="4600"/>
      <value value="9200"/>
      <value value="13939"/>
      <value value="18400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
      <value value="&quot;early&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="200"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
