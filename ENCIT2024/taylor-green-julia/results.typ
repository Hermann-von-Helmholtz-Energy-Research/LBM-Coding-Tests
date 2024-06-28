#let ivory = rgb("#e6deca")
#set page(paper: "a4", fill: ivory)
#set text(font: "Crimson Pro", size: 11pt)
#set par(justify: true, leading: 0.65em)

#show raw: txt => text(fill: rgb("#000048"), weight: 550, txt)

= Computing Environment Description
#v(6pt)

Host `freebird` is a
`Gigabyte` computer `desktop`; with
`Intel Core i3-9100F` CPU model:
`64` bits,
`core 9` gen,
`1` CPU with
`4` cores with
`256KiB` L1 cache (`4x32 KiB` for data and `4x32 KiB` for instructions),
`1024 KiB` L2 cache, and
`6 MiB` L3 cache, running at
`3299 MHz` (avg.); with
`16 GiB` total RAM (
`15.57 GiB` available, and
`4.64 GiB` used by some `234` concurrent processes) on a
`Linux` operating system with a
`6.6.32-1-MANJARO` kernel,
`x86_64` architecture, compiled with
`gcc` version
`14.1.1` with
`60% (default)` of swappiness, with
`16 GiB` of hard disk swap of which
`0.0 MiB (0.0%)` were in use during benchmark tests.

Host `stilo` is a
`Positivo Informatica SA` computer `laptop`; with
`Intel Celeron N3010` CPU model:
`64` bits,
(unknown) gen,
`1` CPU with
`2` cores with
`112KiB` L1 cache (`2x24 KiB` for data and `2x32 KiB` for instructions),
`2 MiB` L2 cache (`2x1024 KiB`), and
no L3 cache, running at
`539 MHz` (avg.); with
`4 GiB` total RAM (
`3.69 GiB` available (est.), and
`2.65 GiB` used by some `213` concurrent processes) on a
`Linux` operating system with a
`5.15.130-1-MANJARO` kernel,
`x86_64` architecture, compiled with
`gcc` version
`13.2.1` with
`60% (default)` of swappiness, with
`4 GiB` of hard disk swap of which
`518.9 MiB (12.7%)` were in use during benchmark tests.

= Benchmarks on 2D Taylor-Green Vortex Decay with `D2Q9` LBM

== C (ISO C99) Codes

The `ISO-C99` code was compiled with `GCC`, version `14.1.1` of `2024-05-22` (Copyright #sym.copyright 2024 Free Software Foundation, Inc.) in each one of the computing environments described above. ...

== Julia Codes

The `julia` language tests were performed with the following characteristics: version `1.10.2` (commit `bd47eca2c8a` official build) reporting CPU frequency of `1.04GHz`, `libopenlibm`, and `libLLVM-15.0.7` on host `stilo`.

#align(center, table(align: center + horizon,
    columns: 7,
    inset: 4pt,
    table.header[Domain size][Time Steps][Precision][Memory][Host][Speed (Mlups)][Relative Speed],
    [32 #sym.times 32], [204'800], [`Float64`], [168.47 KiB], [`stilo`], [3.61 (avg. of 3)], [],
  )
)



