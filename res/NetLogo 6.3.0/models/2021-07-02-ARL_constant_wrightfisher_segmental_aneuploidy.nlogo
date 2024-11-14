breed [tumorcells tumorcell]

globals
[
  Chr1p
  Chr1q
  Chr2p
  Chr2q
  Chr3p
  Chr3q
  Chr4p
  Chr4q
  Chr5p
  Chr5q
  Chr6p
  Chr6q
  Chr7p
  Chr7q
  Chr8p
  Chr8q
  Chr9p
  Chr9q
  Chr10p
  Chr10q
  Chr11p
  Chr11q
  Chr12p
  Chr12q
  Chr13p
  Chr13q
  Chr14p
  Chr14q
  Chr15p
  Chr15q
  Chr16p
  Chr16q
  Chr17p
  Chr17q
  Chr18p
  Chr18q
  Chr19p
  Chr19q
  Chr20p
  Chr20q
  Chr21p
  Chr21q
  Chr22p
  Chr22q
  ChrXp
  ChrXq
  fitness_pop
  x
  tandem_losses
  driver_scores
  abundance_scores
  ploidyavg
  karyolists
  variances
  mkv
]

tumorcells-own
[
  chromosomes
  fitness
  ploidy
  missegs
]

to setup
  clear-all
  reset-ticks
  set-default-shape tumorcells "circle"
  set-periodicity
  setup-progenitor
  setup-simulation
  output-chromosome-state
end

to setup-simulation
  set abundance_scores [0.04780162326 0.04340321001 0.02733655330 0.04244053843 0.02310411791 0.02997560125 0.01238194825 0.03181795548 0.01178442796 0.03787614732 0.02557718800 0.02554399243 0.01795880430 0.03231588906 0.01591727664 0.02549419907 0.01301266411 0.02572656808 0.01122010324 0.02750253116 0.01961858288 0.03629935767 0.01425749805 0.03659811781 0.00000000000 0.02333648691 0.00001659779 0.03792594068 0.00000000000 0.03701306246 0.02383442049 0.01900446480 0.01548573420 0.03553585952 0.00627396305 0.01434048698 0.02159371940 0.02813324702 0.00896280436 0.01526996299 0.00232369002 0.01233215489 0.00013278229 0.02297133562 0.01555212535 0.02499626550]
  set driver_scores [-0.002401775  0.032443615  0.029357171  0.039432666  0.032896947  0.054167358  0.017849091  0.029013244  0.042811664  0.019499344  0.023986192  0.000116254  0.098892842 0.069333137  0.027695642  0.058614269 -0.001294077  0.047026807 -0.036421793  0.011426879  0.038186214  0.018987838  0.055155098  0.062737856  0.000000000 -0.010153924 0.000000000  0.025574386  0.000000000  0.020656601  0.043347363 -0.007144357 -0.008597452  0.043634738  0.005336974 -0.026363166  0.053714158  0.005503382  0.043510246 0.049935929  0.000000000 -0.003309223  0.000000000 -0.005158138  0.000000000  0.000000000]
end

to setup-progenitor
  create-tumorcells initial-number
  ask tumorcells
  [
    set chromosomes n-values 46 [initial-ploidy]
  ]
end

to go
if (ticks = 1) [ reset-timer ]
if ticks = end-ticks or count tumorcells >= end-population or count tumorcells = 0
  [
    show timer
    stop
  ]
  if ticks > cin-period-end [set pmisseg 0]
  reset-population
  mitosis
  tick
  calculate-fitness
  count-tandem-losses
  calculate-stats
  output-chromosome-state
  check-death
end

to calculate-stats
  set variances []
  foreach (n-values 46 [i -> i]) [i -> set variances insert-item i variances variance [item i chromosomes] of tumorcells]
  set mkv mean variances
end

to reset-population
  set karyolists []

  foreach range 46 [i ->
    set karyolists insert-item i karyolists [item i chromosomes] of tumorcells
  ]

  clear-turtles
  create-tumorcells initial-number

  ask tumorcells
  [
    let new_karyotype []
    let i -1

    foreach range 46 [
      set i i + 1
      set new_karyotype insert-item i new_karyotype (one-of item i karyolists)
    ]
    set chromosomes new_karyotype
  ]
end


to count-tandem-losses
set tandem_losses count tumorcells with [
    item 0 chromosomes = 0 or
    item 1 chromosomes = 0 or
    item 2 chromosomes = 0 or
    item 3 chromosomes = 0 or
    item 4 chromosomes = 0 or
    item 5 chromosomes = 0 or
    item 6 chromosomes = 0 or
    item 7 chromosomes = 0 or
    item 8 chromosomes = 0 or
    item 9 chromosomes = 0 or
    item 10 chromosomes = 0 or
    item 11 chromosomes = 0 or
    item 12 chromosomes = 0 or
    item 13 chromosomes = 0 or
    item 14 chromosomes = 0 or
    item 15 chromosomes = 0 or
    item 16 chromosomes = 0 or
    item 17 chromosomes = 0 or
    item 18 chromosomes = 0 or
    item 19 chromosomes = 0 or
    item 20 chromosomes = 0 or
    item 21 chromosomes = 0 or
    item 22 chromosomes = 0 ]
end


to set-periodicity
    if (periodicity = "constant")[
    set cin-period-start 0
    set cin-period-end end-ticks
  ]
      if (periodicity = "early")[
    set cin-period-start 0
    set cin-period-end round(end-ticks * 0.1)
  ]


end

to mitosis
  ask tumorcells[
      let daughterKaryotype n-values 46 [0]
      let chromosomeTracker -1
      foreach chromosomes [i ->
        set chromosomeTracker chromosomeTracker + 1
        foreach range i [
          if (precision(random-float 1)3) < pmisseg [
            let gainloss random 2
            if (remainder chromosomeTracker 2) = 0 [
              if gainloss = 1 [set daughterKaryotype replace-item chromosomeTracker daughterKaryotype ((item chromosomeTracker daughterKaryotype) + 1)]
              if gainloss = 0 [set daughterKaryotype replace-item chromosomeTracker daughterKaryotype ((item chromosomeTracker daughterKaryotype) - 1)]
              if (precision(random-float 1)3) > pbreak [
                if gainloss = 1 [set daughterKaryotype replace-item (chromosomeTracker + 1) daughterKaryotype ((item (chromosomeTracker + 1) daughterKaryotype) + 1)]
                if gainloss = 0 and (item (chromosomeTracker + 1) chromosomes) > 0 [set daughterKaryotype replace-item (chromosomeTracker + 1) daughterKaryotype ((item (chromosomeTracker + 1) daughterKaryotype) - 1)]
              ]
            ]
          ]
        ]
      ]
      set chromosomes (map [[a b] -> a + b] chromosomes daughterKaryotype)
  ]
end

to calculate-fitness
    ask tumorcells[
    set ploidy (sum chromosomes) / 46
    let abundance_fitness_scores (map [[a b] -> (a - ((a * abs(b - ploidy)) / ploidy))] abundance_scores chromosomes)
    let driver_fitness_scores (map [[a b] -> ((b * a) / ploidy)] driver_scores chromosomes)

    if model = "Abundance" [set fitness precision(sum(abundance_fitness_scores) ^ s)3]
    if model = "Driver" [set fitness precision(sum(driver_fitness_scores) ^ s)3]
    if model = "Hybrid" [set fitness precision(((sum(abundance_fitness_scores) + sum(driver_fitness_scores)) / 2) ^ s)3]
    if model = "Neutral" [set fitness 1]
  ]

end

to check-death
  ask tumorcells [
    if 1 / (fitness + 0.001) > random-float 5 [die]
    if min chromosomes <= 0 [die]
    if max chromosomes > 6 [die]
  ]

end

to output-chromosome-state
 calculate-fitness
 set x count tumorcells

  let selected-tumorcells ifelse-value (x >= outputcount)
  [n-of outputcount tumorcells]
  [n-of count tumorcells tumorcells]

  set Chr1p []
  set Chr1q []
  set Chr2p []
  set Chr2q []
  set Chr3p []
  set Chr3q []
  set Chr4p []
  set Chr4q []
  set Chr5p []
  set Chr5q []
  set Chr6p []
  set Chr6q []
  set Chr7p []
  set Chr7q []
  set Chr8p []
  set Chr8q []
  set Chr9p []
  set Chr9q []
  set Chr10p []
  set Chr10q []
  set Chr11p []
  set Chr11q []
  set Chr12p []
  set Chr12q []
  set Chr13p []
  set Chr13q []
  set Chr14p []
  set Chr14q []
  set Chr15p []
  set Chr15q []
  set Chr16p []
  set Chr16q []
  set Chr17p []
  set Chr17q []
  set Chr18p []
  set Chr18q []
  set Chr19p []
  set Chr19q []
  set Chr20p []
  set Chr20q []
  set Chr21p []
  set Chr21q []
  set Chr22p []
  set Chr22q []
  set ChrXp []
  set ChrXq []
  set fitness_pop []

  foreach (range count selected-tumorcells) [i -> set Chr1q insert-item i Chr1q [item 0 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr1p insert-item i Chr1p [item 1 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr2q insert-item i Chr2q [item 2 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr2p insert-item i Chr2p [item 3 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr3p insert-item i Chr3p [item 4 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr3q insert-item i Chr3q [item 5 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr4p insert-item i Chr4p [item 6 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr4q insert-item i Chr4q [item 7 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr5p insert-item i Chr5p [item 8 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr5q insert-item i Chr5q [item 9 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr6p insert-item i Chr6p [item 10 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr6q insert-item i Chr6q [item 11 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr7p insert-item i Chr7p [item 12 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr7q insert-item i Chr7q [item 13 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr8p insert-item i Chr8p [item 14 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr8q insert-item i Chr8q [item 15 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr9p insert-item i Chr9p [item 16 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr9q insert-item i Chr9q [item 17 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr10p insert-item i Chr10p [item 18 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr10q insert-item i Chr10q [item 19 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr11p insert-item i Chr11p [item 20 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr11q insert-item i Chr11q [item 21 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr12p insert-item i Chr12p [item 22 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr12q insert-item i Chr12q [item 23 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr13p insert-item i Chr13p [item 24 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr13q insert-item i Chr13q [item 25 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr14p insert-item i Chr14p [item 26 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr14q insert-item i Chr14q [item 27 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr15p insert-item i Chr15p [item 28 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr15q insert-item i Chr15q [item 29 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr16p insert-item i Chr16p [item 30 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr16q insert-item i Chr16q [item 31 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr17p insert-item i Chr17p [item 32 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr17q insert-item i Chr17q [item 33 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr18p insert-item i Chr18p [item 34 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr18q insert-item i Chr18q [item 35 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr19p insert-item i Chr19p [item 36 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr19q insert-item i Chr19q [item 37 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr20p insert-item i Chr20p [item 38 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr20q insert-item i Chr20q [item 39 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr21p insert-item i Chr21p [item 40 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr21q insert-item i Chr21q [item 41 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr22p insert-item i Chr22p [item 42 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set Chr22q insert-item i Chr22q [item 43 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set ChrXp insert-item i ChrXp [item 44 chromosomes] of item i (sort selected-tumorcells)]
  foreach (range count selected-tumorcells) [i -> set ChrXq insert-item i ChrXq [item 45 chromosomes] of item i (sort selected-tumorcells)]

  foreach (range count selected-tumorcells) [i -> set fitness_pop insert-item i fitness_pop [fitness] of item i (sort selected-tumorcells)]
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
120.0

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
10.0
1
0
Number

SLIDER
179
275
329
308
initial-ploidy
initial-ploidy
1
6
2.0
1
1
NIL
HORIZONTAL

INPUTBOX
30
125
179
185
end-ticks
100.0
1
0
Number

INPUTBOX
31
248
180
308
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
100.0
1
0
Number

INPUTBOX
30
187
179
247
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
"constant" "early"
0

CHOOSER
179
230
328
275
model
model
"Hybrid" "Driver" "Abundance" "Neutral"
0

PLOT
230
428
430
578
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
"default" 1.0 0 -16777216 true "plot mean [fitness] of tumorcells" "plot mean [fitness] of tumorcells"

PLOT
30
428
230
578
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
428
630
578
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
"default" 1.0 0 -16777216 true "" "plot mean [ ploidy ] of tumorcells"

INPUTBOX
180
308
329
368
s
50.0
1
0
Number

INPUTBOX
31
308
180
368
pmisseg
0.1
1
0
Number

PLOT
630
428
830
578
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
"default" 1.0 0 -16777216 true "" "plot tandem_losses"

INPUTBOX
180
368
329
428
pbreak
0.0
1
0
Number

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
NetLogo 6.2.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="2021-06-29" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>Chr1</metric>
    <metric>Chr2</metric>
    <metric>Chr3</metric>
    <metric>Chr4</metric>
    <metric>Chr5</metric>
    <metric>Chr6</metric>
    <metric>Chr7</metric>
    <metric>Chr8</metric>
    <metric>Chr9</metric>
    <metric>Chr10</metric>
    <metric>Chr11</metric>
    <metric>Chr12</metric>
    <metric>Chr13</metric>
    <metric>Chr14</metric>
    <metric>Chr15</metric>
    <metric>Chr16</metric>
    <metric>Chr17</metric>
    <metric>Chr18</metric>
    <metric>Chr19</metric>
    <metric>Chr20</metric>
    <metric>Chr21</metric>
    <metric>Chr22</metric>
    <metric>ChrX</metric>
    <metric>fitness_pop</metric>
    <metric>tandem_losses</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="100"/>
    </enumeratedValueSet>
    <steppedValueSet variable="pmisseg" first="0" step="0.001" last="0.5"/>
    <steppedValueSet variable="s" first="0" step="1" last="100"/>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mitosisProbability">
      <value value="2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-07-01_low_s_long" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>Chr1p</metric>
    <metric>Chr1q</metric>
    <metric>Chr2p</metric>
    <metric>Chr2q</metric>
    <metric>Chr3p</metric>
    <metric>Chr3q</metric>
    <metric>Chr4p</metric>
    <metric>Chr4q</metric>
    <metric>Chr5p</metric>
    <metric>Chr5q</metric>
    <metric>Chr6p</metric>
    <metric>Chr6q</metric>
    <metric>Chr7p</metric>
    <metric>Chr7q</metric>
    <metric>Chr8p</metric>
    <metric>Chr8q</metric>
    <metric>Chr9p</metric>
    <metric>Chr9q</metric>
    <metric>Chr10p</metric>
    <metric>Chr10q</metric>
    <metric>Chr11p</metric>
    <metric>Chr11q</metric>
    <metric>Chr12p</metric>
    <metric>Chr12q</metric>
    <metric>Chr13p</metric>
    <metric>Chr13q</metric>
    <metric>Chr14p</metric>
    <metric>Chr14q</metric>
    <metric>Chr15p</metric>
    <metric>Chr15q</metric>
    <metric>Chr16p</metric>
    <metric>Chr16q</metric>
    <metric>Chr17p</metric>
    <metric>Chr17q</metric>
    <metric>Chr18p</metric>
    <metric>Chr18q</metric>
    <metric>Chr19p</metric>
    <metric>Chr19q</metric>
    <metric>Chr20p</metric>
    <metric>Chr20q</metric>
    <metric>Chr21p</metric>
    <metric>Chr21q</metric>
    <metric>Chr22p</metric>
    <metric>Chr22q</metric>
    <metric>ChrXp</metric>
    <metric>ChrXq</metric>
    <metric>fitness_pop</metric>
    <metric>tandem_losses</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="4500"/>
    </enumeratedValueSet>
    <steppedValueSet variable="pmisseg" first="0" step="0.001" last="0.05"/>
    <steppedValueSet variable="s" first="0" step="2" last="100"/>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pbreak">
      <value value="0.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-07-01_low_s_long_neutral" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>Chr1p</metric>
    <metric>Chr1q</metric>
    <metric>Chr2p</metric>
    <metric>Chr2q</metric>
    <metric>Chr3p</metric>
    <metric>Chr3q</metric>
    <metric>Chr4p</metric>
    <metric>Chr4q</metric>
    <metric>Chr5p</metric>
    <metric>Chr5q</metric>
    <metric>Chr6p</metric>
    <metric>Chr6q</metric>
    <metric>Chr7p</metric>
    <metric>Chr7q</metric>
    <metric>Chr8p</metric>
    <metric>Chr8q</metric>
    <metric>Chr9p</metric>
    <metric>Chr9q</metric>
    <metric>Chr10p</metric>
    <metric>Chr10q</metric>
    <metric>Chr11p</metric>
    <metric>Chr11q</metric>
    <metric>Chr12p</metric>
    <metric>Chr12q</metric>
    <metric>Chr13p</metric>
    <metric>Chr13q</metric>
    <metric>Chr14p</metric>
    <metric>Chr14q</metric>
    <metric>Chr15p</metric>
    <metric>Chr15q</metric>
    <metric>Chr16p</metric>
    <metric>Chr16q</metric>
    <metric>Chr17p</metric>
    <metric>Chr17q</metric>
    <metric>Chr18p</metric>
    <metric>Chr18q</metric>
    <metric>Chr19p</metric>
    <metric>Chr19q</metric>
    <metric>Chr20p</metric>
    <metric>Chr20q</metric>
    <metric>Chr21p</metric>
    <metric>Chr21q</metric>
    <metric>Chr22p</metric>
    <metric>Chr22q</metric>
    <metric>ChrXp</metric>
    <metric>ChrXq</metric>
    <metric>fitness_pop</metric>
    <metric>tandem_losses</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="4500"/>
    </enumeratedValueSet>
    <steppedValueSet variable="pmisseg" first="0" step="0.001" last="0.05"/>
    <enumeratedValueSet variable="s">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Neutral&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pbreak">
      <value value="0.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-07-01_short_broad_s" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>Chr1p</metric>
    <metric>Chr1q</metric>
    <metric>Chr2p</metric>
    <metric>Chr2q</metric>
    <metric>Chr3p</metric>
    <metric>Chr3q</metric>
    <metric>Chr4p</metric>
    <metric>Chr4q</metric>
    <metric>Chr5p</metric>
    <metric>Chr5q</metric>
    <metric>Chr6p</metric>
    <metric>Chr6q</metric>
    <metric>Chr7p</metric>
    <metric>Chr7q</metric>
    <metric>Chr8p</metric>
    <metric>Chr8q</metric>
    <metric>Chr9p</metric>
    <metric>Chr9q</metric>
    <metric>Chr10p</metric>
    <metric>Chr10q</metric>
    <metric>Chr11p</metric>
    <metric>Chr11q</metric>
    <metric>Chr12p</metric>
    <metric>Chr12q</metric>
    <metric>Chr13p</metric>
    <metric>Chr13q</metric>
    <metric>Chr14p</metric>
    <metric>Chr14q</metric>
    <metric>Chr15p</metric>
    <metric>Chr15q</metric>
    <metric>Chr16p</metric>
    <metric>Chr16q</metric>
    <metric>Chr17p</metric>
    <metric>Chr17q</metric>
    <metric>Chr18p</metric>
    <metric>Chr18q</metric>
    <metric>Chr19p</metric>
    <metric>Chr19q</metric>
    <metric>Chr20p</metric>
    <metric>Chr20q</metric>
    <metric>Chr21p</metric>
    <metric>Chr21q</metric>
    <metric>Chr22p</metric>
    <metric>Chr22q</metric>
    <metric>ChrXp</metric>
    <metric>ChrXq</metric>
    <metric>fitness_pop</metric>
    <metric>tandem_losses</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="4500"/>
    </enumeratedValueSet>
    <steppedValueSet variable="pmisseg" first="0" step="0.005" last="0.5"/>
    <steppedValueSet variable="s" first="0" step="2" last="100"/>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pbreak">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-07-12_3000gensv2" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go
if member? ticks [2 25 50 75 100 1000 2000 3000] [
file-open "/smph/users/arlynch2/Revisions/2021-07-02-ARL_constant_wrightfisher_segmental_aneuploidy_3000_gensv2_chromosome_output.csv"
file-print (word behaviorspace-run-number "," ticks "," model "," s "," Chr1p "," Chr1q "," Chr2p "," Chr2q "," Chr3p "," Chr3q "," Chr4p "," Chr4q "," Chr5p "," Chr5q "," Chr6p "," Chr6q "," Chr7p "," Chr7q "," Chr8p "," Chr8q "," Chr9p "," Chr9q "," Chr10p "," Chr10q "," Chr11p "," Chr11q "," Chr12p "," Chr12q "," Chr13p "," Chr13q "," Chr14p "," Chr14q "," Chr15p "," Chr15q "," Chr16p "," Chr16q "," Chr17p "," Chr17q "," Chr18p "," Chr18q "," Chr19p "," Chr19q "," Chr20p "," Chr20q "," Chr21p "," Chr21q "," Chr22p "," Chr22q "," ChrXp "," ChrXq)
file-flush

]</go>
    <metric>count turtles</metric>
    <metric>mean [fitness] of tumorcells</metric>
    <metric>mean [ploidy] of tumorcells</metric>
    <metric>mkv</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="4500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pmisseg">
      <value value="0.003"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="s">
      <value value="1"/>
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="3000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pbreak">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Abundance&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Neutral&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-population">
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-07-12_3000gens" repetitions="5" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <final>file-open "/smph/users/arlynch2/Revisions/2021-07-02-ARL_constant_wrightfisher_segmental_aneuploidy_3000_gens_chromosome_output.csv"
file-print (word behaviorspace-run-number "," Chr1p "," Chr1q "," Chr2p "," Chr2q "," Chr3p "," Chr3q "," Chr4p "," Chr4q "," Chr5p "," Chr5q "," Chr6p "," Chr6q "," Chr7p "," Chr7q "," Chr8p "," Chr8q "," Chr9p "," Chr9q "," Chr10p "," Chr10q "," Chr11p "," Chr11q "," Chr12p "," Chr12q "," Chr13p "," Chr13q "," Chr14p "," Chr14q "," Chr15p "," Chr15q "," Chr16p "," Chr16q "," Chr17p "," Chr17q "," Chr18p "," Chr18q "," Chr19p "," Chr19q "," Chr20p "," Chr20q "," Chr21p "," Chr21q "," Chr22p "," Chr22q "," ChrXp "," ChrXq)
file-flush
file-close</final>
    <metric>count turtles</metric>
    <metric>mean [fitness] of tumorcells</metric>
    <metric>mean [ploidy] of tumorcells</metric>
    <metric>mkv</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="4500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pmisseg">
      <value value="0.003"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="s">
      <value value="1"/>
      <value value="25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="3000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pbreak">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Abundance&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Neutral&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-population">
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-10-01_s1" repetitions="3" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>Chr1p</metric>
    <metric>Chr1q</metric>
    <metric>Chr2p</metric>
    <metric>Chr2q</metric>
    <metric>Chr3p</metric>
    <metric>Chr3q</metric>
    <metric>Chr4p</metric>
    <metric>Chr4q</metric>
    <metric>Chr5p</metric>
    <metric>Chr5q</metric>
    <metric>Chr6p</metric>
    <metric>Chr6q</metric>
    <metric>Chr7p</metric>
    <metric>Chr7q</metric>
    <metric>Chr8p</metric>
    <metric>Chr8q</metric>
    <metric>Chr9p</metric>
    <metric>Chr9q</metric>
    <metric>Chr10p</metric>
    <metric>Chr10q</metric>
    <metric>Chr11p</metric>
    <metric>Chr11q</metric>
    <metric>Chr12p</metric>
    <metric>Chr12q</metric>
    <metric>Chr13p</metric>
    <metric>Chr13q</metric>
    <metric>Chr14p</metric>
    <metric>Chr14q</metric>
    <metric>Chr15p</metric>
    <metric>Chr15q</metric>
    <metric>Chr16p</metric>
    <metric>Chr16q</metric>
    <metric>Chr17p</metric>
    <metric>Chr17q</metric>
    <metric>Chr18p</metric>
    <metric>Chr18q</metric>
    <metric>Chr19p</metric>
    <metric>Chr19q</metric>
    <metric>Chr20p</metric>
    <metric>Chr20q</metric>
    <metric>Chr21p</metric>
    <metric>Chr21q</metric>
    <metric>Chr22p</metric>
    <metric>Chr22q</metric>
    <metric>ChrXp</metric>
    <metric>ChrXq</metric>
    <metric>fitness_pop</metric>
    <metric>tandem_losses</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="4500"/>
    </enumeratedValueSet>
    <steppedValueSet variable="pmisseg" first="0" step="0.001" last="0.05"/>
    <enumeratedValueSet variable="s">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
      <value value="&quot;Neutral&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pbreak">
      <value value="0.5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="2021-10-18_timecourse" repetitions="10" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count tumorcells</metric>
    <metric>Chr1p</metric>
    <metric>Chr1q</metric>
    <metric>Chr2p</metric>
    <metric>Chr2q</metric>
    <metric>Chr3p</metric>
    <metric>Chr3q</metric>
    <metric>Chr4p</metric>
    <metric>Chr4q</metric>
    <metric>Chr5p</metric>
    <metric>Chr5q</metric>
    <metric>Chr6p</metric>
    <metric>Chr6q</metric>
    <metric>Chr7p</metric>
    <metric>Chr7q</metric>
    <metric>Chr8p</metric>
    <metric>Chr8q</metric>
    <metric>Chr9p</metric>
    <metric>Chr9q</metric>
    <metric>Chr10p</metric>
    <metric>Chr10q</metric>
    <metric>Chr11p</metric>
    <metric>Chr11q</metric>
    <metric>Chr12p</metric>
    <metric>Chr12q</metric>
    <metric>Chr13p</metric>
    <metric>Chr13q</metric>
    <metric>Chr14p</metric>
    <metric>Chr14q</metric>
    <metric>Chr15p</metric>
    <metric>Chr15q</metric>
    <metric>Chr16p</metric>
    <metric>Chr16q</metric>
    <metric>Chr17p</metric>
    <metric>Chr17q</metric>
    <metric>Chr18p</metric>
    <metric>Chr18q</metric>
    <metric>Chr19p</metric>
    <metric>Chr19q</metric>
    <metric>Chr20p</metric>
    <metric>Chr20q</metric>
    <metric>Chr21p</metric>
    <metric>Chr21q</metric>
    <metric>Chr22p</metric>
    <metric>Chr22q</metric>
    <metric>ChrXp</metric>
    <metric>ChrXq</metric>
    <metric>fitness_pop</metric>
    <metric>tandem_losses</metric>
    <enumeratedValueSet variable="initial-number">
      <value value="4500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pmisseg">
      <value value="0.001"/>
      <value value="0.01"/>
      <value value="0.0042"/>
      <value value="0.0046"/>
      <value value="0.0051"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="s">
      <value value="0"/>
      <value value="1"/>
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-ploidy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="end-ticks">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="periodicity">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="model">
      <value value="&quot;Hybrid&quot;"/>
      <value value="&quot;Driver&quot;"/>
      <value value="&quot;Abundance&quot;"/>
      <value value="&quot;Neutral&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="outputcount">
      <value value="300"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pbreak">
      <value value="0.5"/>
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
